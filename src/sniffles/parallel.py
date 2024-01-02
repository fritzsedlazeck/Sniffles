#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 27.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#
import copy
import logging
from argparse import Namespace

from dataclasses import dataclass
import gc
import math
from typing import Optional

import pysam

from sniffles import leadprov
from sniffles import cluster
from sniffles import sv
from sniffles import postprocessing
from sniffles import snf
from sniffles.result import Result


@dataclass
class Task:
    """
    A task is a generic unit of work sent to a child process to be worked on in parallel. Must be pickleable.
    """
    id: int
    sv_id: int
    contig: str
    start: int
    end: int
    config: Namespace
    assigned_process_id: Optional[int] = None
    lead_provider: leadprov.LeadProvider = None
    bam: object = None
    tandem_repeats: list = None
    genotype_svs: list = None
    _logger = None

    @property
    def logger(self) -> logging.Logger:
        if self._logger is None:
            self._logger = logging.getLogger(f'sniffles.progress')

        return self._logger

    def execute(self) -> Result:
        """
        Execute this Task, returning a Result object
        """
        return Result(self, [], 0)

    def build_leadtab(self, config):
        assert (self.lead_provider is None)

        if config.input_is_cram and config.reference is not None:
            self.bam = pysam.AlignmentFile(config.input, config.input_mode, require_index=True, reference_filename=config.reference)
        else:
            self.bam = pysam.AlignmentFile(config.input, config.input_mode, require_index=True)
        self.lead_provider = leadprov.LeadProvider(config, self.id * config.task_read_id_offset_mult)
        externals = self.lead_provider.build_leadtab(self.contig, self.start, self.end, self.bam)
        return externals, self.lead_provider.read_count

    def call_candidates(self, keep_qc_fails, config):
        candidates = []
        for svtype in sv.TYPES:
            for svcluster in cluster.resolve(svtype, self.lead_provider, config, self.tandem_repeats):
                for svcall in sv.call_from(svcluster, config, keep_qc_fails, self):
                    if config.dev_trace_read is not False:
                        cluster_has_read = False
                        for ld in svcluster.leads:
                            if ld.read_qname == config.dev_trace_read:
                                cluster_has_read = True
                        if cluster_has_read:
                            import copy
                            svcall_copy = copy.deepcopy(svcall)
                            svcall_copy.postprocess = None
                            print(f"[DEV_TRACE_READ] [3/4] [Task.call_candidates] Read {config.dev_trace_read} -> Cluster {svcluster.id} -> preliminary SVCall {svcall_copy}")
                    candidates.append(svcall)

        self.coverage_average_fwd, self.coverage_average_rev = postprocessing.coverage(candidates, self.lead_provider, config)
        self.coverage_average_total = self.coverage_average_fwd + self.coverage_average_rev
        return candidates

    def finalize_candidates(self, candidates, keep_qc_fails, config):
        passed = []
        for svcall in candidates:
            svcall.qc = svcall.qc and postprocessing.qc_sv(svcall, config)
            if not keep_qc_fails and not svcall.qc:
                continue
            svcall.qc = svcall.qc and postprocessing.qc_sv_support(svcall, self.coverage_average_total, config)
            if not keep_qc_fails and not svcall.qc:
                continue

            postprocessing.annotate_sv(svcall, config)

            svcall.qc = svcall.qc and postprocessing.qc_sv_post_annotate(svcall, config)

            if config.dev_trace_read:
                cluster_has_read = False
                for ld in svcall.postprocess.cluster.leads:
                    if ld.read_qname == config.dev_trace_read:
                        cluster_has_read = True
                if cluster_has_read:
                    import copy
                    svcall_copy = copy.deepcopy(svcall)
                    svcall_copy.postprocess = None
                    print(f"[DEV_TRACE_READ] [4/4] [Task.finalize_candidates] Read {config.dev_trace_read} -> Cluster {svcall.postprocess.cluster.id} -> finalized SVCall, QC={svcall_copy.qc}: {svcall_copy}")

            if not keep_qc_fails and not svcall.qc:
                continue

            svcall.finalize()  # Remove internal information (not written to output) before sending to mainthread for VCF writing
            passed.append(svcall)
        return passed


class CallTask(Task):
    """
    """


class GenotypeTask(Task):
    ...


class CombineTask(Task):
    """
    Task to merge/combine multiple SNF files into one.
    """
    MIN_BLOCKS_PER_THREAD = 100
    emit_first_block = True

    block_indices: list[int] = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.generate_blocks()

    def generate_blocks(self):
        """
        Generate a set of blocks
        """
        self.block_indices = list(range(self.start, self.end + self.config.snf_block_size, self.config.snf_block_size))

    def __str__(self):
        if len(self.block_indices) > 0:
            return f'''Task {self.id} [{self.start} ({self.block_indices[0]}) .. {self.end} ({self.block_indices[-1]})]'''
        else:
            return f'Task {self.id} [no blocks available]'

    def clone(self, first_block: int, block_count: int, new_id: int = None, emit_first: bool = True) -> 'CombineTask':
        """
        Clone this task with a subset of the blocks to be processed.
        """
        obj = copy.copy(self)
        if new_id is not None:
            obj.id = new_id
        obj.block_indices = self.block_indices[first_block:first_block+block_count]
        obj.emit_first_block = emit_first
        obj.start = obj.block_indices[0]
        obj.end = obj.block_indices[-1] + obj.config.snf_block_size
        return obj

    def scatter(self) -> list['CombineTask']:
        """
        Scatter this task to a number of tasks to be run in parallel, on block level:
        - Distribute blocks equally to tasks, each one working on consecutive blocks
        - First task gets any blocks that can not be evenly distributed
        - Tasks other than the first one have to process their first block for kept groups, and
          should not emit calls for this block (as these will be emitted by the previous task)
        """
        if self.config.threads > 1:
            if (nBlocks := len(self.block_indices)) >= self.MIN_BLOCKS_PER_THREAD*2 and not self.config.dev_disable_interblock_threads:
                parallel_tasks = min(self.config.threads, int(nBlocks / self.MIN_BLOCKS_PER_THREAD))
                blocks_per_task = int(nBlocks / parallel_tasks)
                blocks_for_first_task = nBlocks - blocks_per_task * (parallel_tasks - 1)
                if parallel_tasks > 1:
                    return [self.clone(0, blocks_for_first_task)] + [
                        self.clone(
                            blocks_for_first_task+i*blocks_per_task-1,  # Start one block earlier...
                            blocks_per_task,
                            emit_first=False,  # ...but dont emit it
                            new_id=self.id+i+1
                        ) for i in range(parallel_tasks-1)
                    ]

        return [self]

    def execute(self):
        samples_headers_snf = {}
        for snf_info in self.config.snf_input_info:
            snf_in = snf.LazySNFile(self.config, open(snf_info["filename"], "rb"), filename=snf_info["filename"])
            snf_in.read_header()
            samples_headers_snf[snf_info["internal_id"]] = snf_in

            if self.config.combine_close_handles:
                snf_in.close()

        svcalls = []

        # block_groups_keep_threshold=5000
        # TODO: Parameterize
        bin_min_size = self.config.combine_min_size
        bin_max_candidates = max(25, int(len(self.config.snf_input_info) * 0.5))
        overlap_abs = self.config.combine_overlap_abs

        sample_internal_ids = set(samples_headers_snf.keys())

        #
        # Load candidate SVs from all samples for each block separately and cluster them based on start position
        #
        candidates_processed = 0
        groups_keep = {svtype: list() for svtype in sv.TYPES}

        for cur, block_index in enumerate(self.block_indices):  # iterate over all blocks
            self.logger.info(f'Processing block {cur+1}/{len(self.block_indices)} (active calls: {sv.SVCall._counter} groups: {sv.SVGroup._counter})')
            samples_blocks = {}
            for sample_internal_id, sample_snf in samples_headers_snf.items():
                blocks = sample_snf.read_blocks(self.contig, block_index)
                samples_blocks[sample_internal_id] = blocks
                # sample_internal_id is the number of the processed file
                # blocks is a list[dict[str, list[SVCall]]] {'INS': [...], 'DEL': [...], ...}

            for svtype in sv.TYPES:
                bins = {}
                # svcandidates=[]
                for sample_internal_id in samples_headers_snf.keys():  # fetch current block for each file
                    blocks = samples_blocks[sample_internal_id]
                    if blocks is None:
                        continue
                    for block in blocks:  # usually only 1 block
                        for cand in block[svtype]:
                            # if config.combine_pass_only and (cand.qc==False or cand.filter!="PASS"):
                            #    continue

                            cand.sample_internal_id = sample_internal_id

                            bin = int(cand.pos / bin_min_size) * bin_min_size
                            if bin not in bins:
                                bins[bin] = [cand]
                            else:
                                bins[bin].append(cand)
                        candidates_processed += len(block[svtype])

                if len(bins) == 0:
                    continue

                size = 0
                svcands = []
                keep = groups_keep[svtype]
                sorted_bins = sorted(bins)
                last_bin = sorted_bins[-1]
                for curr_bin in sorted_bins:
                    svcands.extend(bins[curr_bin])  # here SVCalls from bins are collected...
                    size += bin_min_size

                    if (not self.config.combine_exhaustive and len(svcands) >= bin_max_candidates) or curr_bin == last_bin:
                        if len(svcands) == 0:
                            size = 0
                            continue

                        svgroups = cluster.resolve_block_groups(svtype, svcands, keep, self.config)
                        groups_call = []
                        keep = []
                        for group in svgroups:
                            coverage_bin = int(
                                group.pos_mean / self.config.coverage_binsize_combine) * self.config.coverage_binsize_combine
                            for non_included_sample in sample_internal_ids - group.included_samples:
                                if samples_blocks[non_included_sample] is not None and coverage_bin in samples_blocks[non_included_sample][0]["_COVERAGE"]:
                                    coverage = samples_blocks[non_included_sample][0]["_COVERAGE"][coverage_bin]
                                else:
                                    coverage = 0
                                if non_included_sample in group.coverages_nonincluded:
                                    group.coverages_nonincluded[non_included_sample] = max(
                                        coverage,
                                        group.coverages_nonincluded[non_included_sample]
                                    )
                                else:
                                    group.coverages_nonincluded[non_included_sample] = coverage

                            if abs(group.pos_mean - curr_bin) < max(size * 0.5, overlap_abs):
                                keep.append(group)
                            else:
                                groups_call.append(group)

                        if cur > 0 or self.emit_first_block:
                            if cur == 1 and not self.emit_first_block and len(self.block_indices) > 1:
                                # If we're not emitting the first block
                                svcalls.extend(call for call in sv.call_groups(groups_call, self.config, self) if not call.pos < self.block_indices[1])
                            else:
                                svcalls.extend(sv.call_groups(groups_call, self.config, self))

                        size = 0
                        svcands = []

                groups_keep[svtype] = keep

        for svtype in groups_keep:
            svcalls.extend(sv.call_groups(groups_keep[svtype], self.config, self))

        from sniffles.result import CombineResult
        return CombineResult(self, svcalls, candidates_processed)


@dataclass
class Process:
    id: int
    process: object = None
    pipe_main: object = None
    externals: list = None


def Main(proc_id, config, pipe):
    try:
        if config.dev_profile:
            import cProfile
            cProfile.runctx("Main_Internal(proc_id,config,pipe)", globals(), locals(), sort="cumulative")
        else:
            Main_Internal(proc_id, config, pipe)
    except Exception as e:
        pipe.send(["worker_exception", ""])
        raise e


def Main_Internal(proc_id, config, pipe):
    tasks = {}
    while True:
        command, arg = pipe.recv()

        if command == "call_sample":
            task = arg

            if config.snf is not None or config.no_qc:
                qc = False
            else:
                qc = True

            _, read_count = task.build_leadtab(config)
            svcandidates = task.call_candidates(qc, config)
            svcalls = task.finalize_candidates(svcandidates, not qc, config)
            if not config.no_qc:
                svcalls = [s for s in svcalls if s.qc]

            from sniffles.result import CallResult
            result = CallResult(task, svcalls, read_count)

            if config.snf is not None:  # and len(svcandidates):
                snf_filename = f"{config.snf}.tmp_{task.id}.snf"

                with open(snf_filename, "wb") as handle:
                    snf_out = snf.SNFile(config, handle)
                    for cand in svcandidates:
                        snf_out.store(cand)
                    snf_out.annotate_block_coverages(task.lead_provider)
                    snf_out.write_and_index()
                    handle.close()
                result.snf_filename = snf_filename
                result.snf_index = snf_out.get_index()
                result.snf_total_length = snf_out.get_total_length()
                result.snf_candidate_count = len(svcandidates)
                result.has_snf = True

            result.coverage_average_total = task.coverage_average_total
            pipe.send(["return_call_sample", result])
            del task
            gc.collect()

        elif command == "genotype_vcf":
            task = arg

            qc = False
            _, read_count = task.build_leadtab(config)
            svcandidates = task.call_candidates(qc, config=config)
            svcalls = task.finalize_candidates(svcandidates, not qc, config=config)

            binsize = 5000
            binedge = int(binsize / 10)
            genotype_svs_svtypes_bins = {svtype: {} for svtype in sv.TYPES}
            for genotype_sv in task.genotype_svs:
                genotype_sv.genotype_match_sv = None
                genotype_sv.genotype_match_dist = math.inf

                if genotype_sv.svtype not in genotype_svs_svtypes_bins:
                    # TODO: Warn about unsupported SVTYPE
                    continue

                bins = [int(genotype_sv.pos / binsize) * binsize]
                if genotype_sv.pos % binsize < binedge:
                    bins.append((int(genotype_sv.pos / binsize) - 1) * binsize)
                if genotype_sv.pos % binsize > binsize - binedge:
                    bins.append((int(genotype_sv.pos / binsize) + 1) * binsize)

                for bin in bins:
                    if bin not in genotype_svs_svtypes_bins[genotype_sv.svtype]:
                        genotype_svs_svtypes_bins[genotype_sv.svtype][bin] = []
                    genotype_svs_svtypes_bins[genotype_sv.svtype][bin].append(genotype_sv)

            for cand in svcandidates:
                bin = int(cand.pos / binsize) * binsize
                if bin not in genotype_svs_svtypes_bins[cand.svtype]:
                    continue
                if cand.svtype == "BND":
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        dist = abs(genotype_sv.pos - cand.pos)
                        # if minlen>0 and dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd * 2:
                        if dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd:
                            if cand.bnd_info.mate_contig == genotype_sv.bnd_info.mate_contig:
                                genotype_sv.genotype_match_sv = cand
                                genotype_sv.genotype_match_dist = dist
                else:
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        dist = abs(genotype_sv.pos - cand.pos) + abs(abs(genotype_sv.svlen) - abs(cand.svlen))
                        minlen = float(min(abs(genotype_sv.svlen), abs(cand.svlen)))
                        if minlen > 0 and dist < genotype_sv.genotype_match_dist and dist <= config.combine_match * math.sqrt(
                                minlen) and dist <= config.combine_match_max:
                            genotype_sv.genotype_match_sv = cand
                            genotype_sv.genotype_match_dist = dist

            postprocessing.coverage(task.genotype_svs, task.lead_provider, config)

            # Determine genotypes for unmatched input SVs
            for svcall in task.genotype_svs:
                coverage_list = [svcall.coverage_start, svcall.coverage_center, svcall.coverage_end]
                coverage_list = [c for c in coverage_list if c is not None]
                if len(coverage_list) == 0:
                    return
                coverage = round(sum(coverage_list) / len(coverage_list))
                svcall.genotypes = {}
                if coverage > 0:
                    svcall.genotypes[0] = (0, 0, 0, coverage, 0, None)
                else:
                    svcall.genotypes[0] = config.genotype_none

            from sniffles.result import GenotypeResult
            result = GenotypeResult(task, task.genotype_svs, read_count)

            pipe.send(["return_genotype_vcf", result])
            del task
            gc.collect()

        elif command == "combine":
            task = arg

            result = task.execute()

            pipe.send(["return_combine", result])
            del task
            gc.collect()

        elif command == "finalize":
            return
