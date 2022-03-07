#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 27.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass
import gc
import math

import pysam

from sniffles import leadprov
from sniffles import cluster
from sniffles import sv
from sniffles import postprocessing
from sniffles import snf
from sniffles import util

@dataclass
class Task:
    id: int
    sv_id: int
    contig: str
    start: int
    end: int
    assigned_process_id: int=None
    lead_provider: object=None
    bam: object=None
    tandem_repeats: list=None
    genotype_svs: list=None

    def build_leadtab(self,config):
        assert(self.lead_provider==None)

        self.bam=pysam.AlignmentFile(config.input, config.input_mode, require_index=True)
        self.lead_provider=leadprov.LeadProvider(config,self.id*config.task_read_id_offset_mult)
        externals=self.lead_provider.build_leadtab(self.contig,self.start,self.end,self.bam)
        return externals,self.lead_provider.read_count

    def call_candidates(self,keep_qc_fails,config):
        candidates=[]
        for svtype in sv.TYPES:
            for svcluster in cluster.resolve(svtype,self.lead_provider,config,self.tandem_repeats):
                for svcall in sv.call_from(svcluster,config,keep_qc_fails,self):
                    candidates.append(svcall)

        self.coverage_average_fwd,self.coverage_average_rev=postprocessing.coverage(candidates,self.lead_provider,config)
        self.coverage_average_total=self.coverage_average_fwd+self.coverage_average_rev
        return candidates

    def finalize_candidates(self,candidates,keep_qc_fails,config):
        passed=[]
        for svcall in candidates:
            svcall.qc=svcall.qc and postprocessing.qc_sv(svcall,config)
            if not keep_qc_fails and not svcall.qc:
                continue
            svcall.qc=svcall.qc and postprocessing.qc_sv_support(svcall,self.coverage_average_total,config)
            if not keep_qc_fails and not svcall.qc:
                continue

            postprocessing.annotate_sv(svcall,config)

            svcall.qc=svcall.qc and postprocessing.qc_sv_post_annotate(svcall,config)
            if not keep_qc_fails and not svcall.qc:
                continue

            svcall.finalize() #Remove internal information (not written to output) before sending to mainthread for VCF writing
            passed.append(svcall)
        return passed

    def combine(self,config):
        samples_headers_snf={}
        for snf_info in config.snf_input_info:
            snf_in=snf.SNFile(config,open(snf_info["filename"],"rb"),filename=snf_info["filename"])
            snf_in.read_header()
            samples_headers_snf[snf_info["internal_id"]]=(snf_info["filename"],snf_in.header,snf_in)

            if not config.snf_combine_keep_open:
                snf_in.close()

        svcalls=[]

        #block_groups_keep_threshold=5000
        #TODO: Parameterize
        bin_min_size=config.combine_min_size
        bin_max_candidates=max(25,int(len(config.snf_input_info)*0.5))
        overlap_abs=config.combine_overlap_abs

        sample_internal_ids=set()
        for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
            sample_internal_ids.add(sample_internal_id)

        #
        # Load candidate SVs from all samples for each block separately and cluster them based on start position
        #
        had=False
        candidates_processed=0
        svtypes_candidates_bins={svtype: {} for svtype in sv.TYPES}
        groups_keep={svtype:list() for svtype in sv.TYPES}
        for block_index in range(self.start,self.end+config.snf_block_size,config.snf_block_size):
            samples_blocks={}
            for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
                blocks=sample_snf.read_blocks(self.contig,block_index)
                samples_blocks[sample_internal_id]=blocks

            for svtype in sv.TYPES:
                bins={}
                #svcandidates=[]
                for sample_internal_id,(sample_filename,sample_header,sample_snf) in samples_headers_snf.items():
                    blocks=samples_blocks[sample_internal_id]
                    if blocks==None:
                        continue
                    for block in blocks:
                        for cand in block[svtype]:
                            #if config.combine_pass_only and (cand.qc==False or cand.filter!="PASS"):
                            #    continue

                            cand.sample_internal_id=sample_internal_id

                            bin=int(cand.pos/bin_min_size)*bin_min_size
                            if not bin in bins:
                                bins[bin]=[cand]
                            else:
                                bins[bin].append(cand)
                        candidates_processed+=len(block[svtype])

                if len(bins)==0:
                    continue

                curr_bin=0
                size=0
                svcands=[]
                keep=groups_keep[svtype]
                sorted_bins=sorted(bins)
                last_bin=sorted_bins[-1]
                for curr_bin in sorted_bins:
                    svcands.extend(bins[curr_bin])
                    size+=bin_min_size

                    if (not config.combine_exhaustive and len(svcands) >= bin_max_candidates) or curr_bin == last_bin:
                        if len(svcands)==0:
                            size=0
                            continue
                        prevkept=len(keep)
                        svgroups=cluster.resolve_block_groups(svtype,svcands,keep,config)
                        groups_call=[]
                        keep=[]
                        for group in svgroups:
                            coverage_bin=int(group.pos_mean/config.coverage_binsize_combine)*config.coverage_binsize_combine
                            for non_included_sample in sample_internal_ids-group.included_samples:
                                if samples_blocks[non_included_sample]!=None and coverage_bin in samples_blocks[non_included_sample][0]["_COVERAGE"]:
                                    coverage=samples_blocks[non_included_sample][0]["_COVERAGE"][coverage_bin]
                                else:
                                    coverage=0
                                if non_included_sample in group.coverages_nonincluded:
                                    group.coverages_nonincluded[non_included_sample]=max(coverage,group.coverages_nonincluded[non_included_sample])
                                else:
                                    group.coverages_nonincluded[non_included_sample]=coverage

                            if abs(group.pos_mean - curr_bin) < max(size*0.5,overlap_abs):
                                keep.append(group)
                            else:
                                groups_call.append(group)
                        svcalls.extend(sv.call_groups(groups_call,config,self))
                        size=0
                        svcands=[]

                groups_keep[svtype]=keep

        for svtype in groups_keep:
            svcalls.extend(sv.call_groups(groups_keep[svtype],config,self))

        return svcalls,candidates_processed

@dataclass
class Process:
    id: int
    process: object=None
    pipe_main: object=None
    externals: list=None

def Main(proc_id,config,pipe):
    try:
        if config.dev_profile:
            import cProfile
            cProfile.runctx("Main_Internal(proc_id,config,pipe)",globals(),locals(),sort="cumulative")
        else:
            Main_Internal(proc_id,config,pipe)
    except Exception as e:
        pipe.send(["worker_exception",""])
        raise e

def Main_Internal(proc_id,config,pipe):
    tasks={}
    while True:
        command,arg=pipe.recv()

        if command=="call_sample":
            task=arg
            result={}

            if config.snf != None or config.no_qc:
                qc=False
            else:
                qc=True

            _,read_count=task.build_leadtab(config)
            svcandidates=task.call_candidates(qc,config)
            svcalls=task.finalize_candidates(svcandidates,not qc,config)
            if config.no_qc:
                result["svcalls"]=svcalls
            else:
                result["svcalls"]=[s for s in svcalls if s.qc]

            if config.snf != None: # and len(svcandidates):
                snf_filename=f"{config.snf}.tmp_{task.id}.snf"

                with open(snf_filename,"wb") as handle:
                    snf_out=snf.SNFile(config,handle)
                    for cand in svcandidates:
                        snf_out.store(cand)
                    snf_out.annotate_block_coverages(task.lead_provider)
                    snf_out.write_and_index()
                    handle.close()
                result["snf_filename"]=snf_filename
                result["snf_index"]=snf_out.get_index()
                result["snf_total_length"]=snf_out.get_total_length()
                result["snf_candidate_count"]=len(svcandidates)
                result["has_snf"]=True
            else:
                result["has_snf"]=False

            #if config.vcf != None:
            #    svcalls=task.finalize_candidates(svcandidates,config)
            #    result["svcalls"]=svcalls

            result["task_id"]=task.id
            result["processed_read_count"]=read_count
            result["coverage_average_total"]=task.coverage_average_total
            pipe.send(["return_call_sample",result])
            del task
            gc.collect()

        elif command=="genotype_vcf":
            task=arg
            result={}

            qc=False
            _,read_count=task.build_leadtab(config)
            svcandidates=task.call_candidates(qc,config=config)
            svcalls=task.finalize_candidates(svcandidates,not qc,config=config)

            binsize=5000
            binedge=int(binsize/10)
            genotype_svs_svtypes_bins={svtype:{} for svtype in sv.TYPES}
            for genotype_sv in task.genotype_svs:
                genotype_sv.genotype_match_sv=None
                genotype_sv.genotype_match_dist=math.inf

                if not genotype_sv.svtype in genotype_svs_svtypes_bins:
                    #TODO: Warn about unsupported SVTYPE
                    continue

                bins=[int(genotype_sv.pos/binsize)*binsize]
                if genotype_sv.pos%binsize < binedge:
                    bins.append((int(genotype_sv.pos/binsize)-1)*binsize)
                if genotype_sv.pos%binsize > binsize-binedge:
                    bins.append((int(genotype_sv.pos/binsize)+1)*binsize)

                for bin in bins:
                    if not bin in genotype_svs_svtypes_bins[genotype_sv.svtype]:
                         genotype_svs_svtypes_bins[genotype_sv.svtype][bin]=[]
                    genotype_svs_svtypes_bins[genotype_sv.svtype][bin].append(genotype_sv)

            for cand in svcandidates:
                bin=int(cand.pos/binsize)*binsize
                if not bin in genotype_svs_svtypes_bins[cand.svtype]:
                    continue
                if cand.svtype=="BND":
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        dist=abs(genotype_sv.pos - cand.pos)
                        #if minlen>0 and dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd * 2:
                        if dist < genotype_sv.genotype_match_dist and dist <= config.cluster_merge_bnd:
                            if cand.bnd_info.mate_contig==genotype_sv.bnd_info.mate_contig:
                                genotype_sv.genotype_match_sv=cand
                                genotype_sv.genotype_match_dist=dist
                else:
                    for genotype_sv in genotype_svs_svtypes_bins[cand.svtype][bin]:
                        dist=abs(genotype_sv.pos - cand.pos) + abs(abs(genotype_sv.svlen) - abs(cand.svlen))
                        minlen=float(min(abs(genotype_sv.svlen),abs(cand.svlen)))
                        if minlen>0 and dist < genotype_sv.genotype_match_dist and dist <= config.combine_match * math.sqrt(minlen) and dists <= config.combine_match_max:
                            genotype_sv.genotype_match_sv=cand
                            genotype_sv.genotype_match_dist=dist

            postprocessing.coverage(task.genotype_svs,task.lead_provider,config)

            #Determine genotypes for unmatched input SVs
            for svcall in task.genotype_svs:
                coverage_list=[svcall.coverage_start,svcall.coverage_center,svcall.coverage_end]
                coverage_list=[c for c in coverage_list if c!=None]
                if len(coverage_list)==0:
                    return
                coverage=round(sum(coverage_list)/len(coverage_list))
                svcall.genotypes={}
                if coverage>0:
                    svcall.genotypes[0]=(0,0,0,coverage,0)
                else:
                    svcall.genotypes[0]=config.genotype_none

            result={}
            result["task_id"]=task.id
            result["processed_read_count"]=read_count
            result["svcalls"]=task.genotype_svs
            pipe.send(["return_genotype_vcf",result])
            del task
            gc.collect()

        elif command=="combine":
            task=arg
            result={}

            result["svcalls"],candidates_processed=task.combine(config)

            result["task_id"]=task.id
            result["processed_read_count"]=candidates_processed
            pipe.send(["return_combine",result])
            del task
            gc.collect()

        elif command=="finalize":
            return
