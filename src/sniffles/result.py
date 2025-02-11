#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     05.12.2023
# Author:      Hermann Romanek
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
import logging
import os
import pickle

from sniffles.config import SnifflesConfig
from sniffles.snf import SNFile
from sniffles.sv import SVCall
from sniffles.vcf import VCF


log = logging.getLogger(__name__)


class Result:
    """
    A generic result of a sniffles task executed by a worker process.
    """
    processed_read_count: int
    task_id: int
    run_id: str
    contig: str
    processed_read_count: int
    svcalls: list[SVCall]
    svcount: int
    error: bool = False

    def __init__(self, task: 'Task', svcalls: list[SVCall], candidates_processed: int):
        self.task_id = task.id
        self.run_id = task.config.run_id
        self.contig = task.contig
        self.processed_read_count = candidates_processed
        self.svcount = len(svcalls)
        self.store_calls(svcalls)

    def store_calls(self, svcalls: list[SVCall]) -> None:
        """
        Store the given calls to this result
        """
        self.svcalls = svcalls

    def emit(self, vcf_out: VCF = None, **kwargs) -> int:
        """
        Emit this result to a file. Returns the number of records written.
        """
        if vcf_out is not None:
            calls = self.svcalls
            if calls:
                for call in calls:
                    vcf_out.write_call(call)
                log.debug(f"Wrote {len(calls)} calls from {self} to VCF.")
            else:
                log.debug(f'No calls for {self}')
            return len(calls)
        else:
            log.debug(f'No vcf output file specified.')
            return 0


class CallResult(Result):
    coverage_average_total: float
    has_snf = False
    snf_filename = None
    snf_index = None
    snf_total_length = None
    snf_candidate_count = None

    def emit(self, **kwargs) -> int:
        res = super().emit(**kwargs)
        if snf_out := kwargs.get('snf_out'):
            snf_out.add_result(self)
        return res


class GenotypeResult(Result):
    """
    Result of a genotyping run
    """

    def emit(self, vcf_out: VCF = None, snf_out: SNFile = None, **kwargs) -> int:
        if vcf_out is not None:
            genotype_lineindex_order = kwargs['genotype_lineindex_order']
            genotype_lineindex_svcalls_returned = {}

            for svcall in self.svcalls:
                genotype_lineindex_svcalls_returned[svcall.raw_vcf_line_index] = svcall

            for lineindex in genotype_lineindex_order:
                if lineindex in genotype_lineindex_svcalls_returned:
                    vcf_out.rewrite_genotype(genotype_lineindex_svcalls_returned[lineindex])

        return len(self.svcalls)


class CombineResult(Result):
    """
    Result of a combine run for one task, simple variant with calls in memory. Must be pickleable.
    """
    def store_calls(self, svcalls: list[SVCall]) -> None:
        from sniffles.config import SnifflesConfig
        if SnifflesConfig.GLOBAL.sort:
            svcalls = sorted(svcalls, key=lambda call: call.pos)

        try:
            self.svcalls.extend(svcalls)
        except AttributeError:
            self.svcalls = list(svcalls)

    def finalize(self):
        if SnifflesConfig.GLOBAL.sort:
            self.svcalls.sort(key=lambda call: call.pos)

    def emit(self, vcf_out: VCF = None, **kwargs) -> int:
        res = super().emit(vcf_out=vcf_out, **kwargs)
        from sniffles.config import SnifflesConfig
        if (config := SnifflesConfig.GLOBAL).dev_population_snf:
            from sniffles.snfp import PopulationSNF
            PopulationSNF.create(config, self.svcalls)

        return res

    def __str__(self):
        return f'CombineResult #{self.task_id}'


class CombineResultTmpFile(CombineResult):
    """
    Result of a combine run, with calls in a temporary file instead of memory
    """
    _highest_position_call: int = -1  # maximum position of last emitted call, for sorting
    unsorted: bool = False
    _initialized: bool = False

    @property
    def tmpfile_name(self) -> str:
        return os.path.join(SnifflesConfig.GLOBAL.tmp_dir, f'result-{self.run_id}-{self.task_id:04}.part.vcf')

    @property
    def tmpfile_unsorted(self) -> str:
        """
        Name of this task/results temporary file for unsorted calls
        """
        return os.path.join(SnifflesConfig.GLOBAL.tmp_dir, f'result-{self.run_id}-{self.task_id:04}-unsorted.part.vcf')

    def store_calls(self, svcalls):
        from sniffles.config import SnifflesConfig
        offset = 0

        if not self._initialized:
            if os.path.exists(self.tmpfile_name):
                self.cleanup()
            self._initialized = True

        if SnifflesConfig.GLOBAL.sort and svcalls:
            svcalls = list(sorted(svcalls, key=lambda call: call.pos))

            while offset < len(svcalls) and svcalls[offset].pos < self._highest_position_call:
                log.debug(f'Unsorted call detected: {self._highest_position_call} > {svcalls[0]}')
                offset += 1

            if offset > 0:  # unsorted calls
                self.unsorted = True
                with open(self.tmpfile_unsorted, 'a') as f:
                    vcf = VCF(SnifflesConfig.GLOBAL, f)
                    for call in svcalls[:offset]:
                        vcf.write_call(call)

            self._highest_position_call = svcalls[-1].pos

        log.debug(f'Storing {len(svcalls)} calls for task {self.task_id} to {self.tmpfile_name}')

        with open(self.tmpfile_name, 'a') as f:
            vcf = VCF(SnifflesConfig.GLOBAL, f)
            vcf.open_reference(generate_index=False)
            for call in svcalls[offset:]:
                vcf.write_call(call)

    def finalize(self):
        ...

    def _emit_sorted(self) -> int:
        # contig, position, _ = line.split('\t', 2)
        # position = int(position)
        # while position > self._unsorted_calls[-1].pos:
        #     vcf_out.write_call(self._unsorted_calls.pop())
        #     n += 1
        return 0

    def emit(self, vcf_out: VCF = None, **kwargs) -> int:
        n = 0
        with open(self.tmpfile_name, 'r') as f:
            for line in f:
                vcf_out.handle.write(line)
                n += 1

        vcf_out.call_count += n
        self.cleanup()
        return n

    def cleanup(self):
        try:
            os.unlink(self.tmpfile_name)
        except FileNotFoundError:
            ...


class CombineResultTmpFilePopulationSNF(CombineResultTmpFile):
    """
    Combine result with temporary files plus population SNF creation
    """
    @property
    def snf_filename(self) -> str:
        return os.path.join(SnifflesConfig.GLOBAL.tmp_dir, f'result-{self.run_id}-{self.task_id:04}.part.snf')

    has_snf = True
    snf_index: dict
    snf_total_length: int
    snf_candidate_count: int = 0
    svcalls: list = None

    def store_calls(self, svcalls: list[SVCall]) -> None:
        super().store_calls(svcalls)

        if self.svcalls is None:
            self.svcalls = []
        self.svcalls.extend(svcalls)

    def finalize(self):
        with open(self.snf_filename, 'wb') as handle:
            from sniffles.snfp import PopulationSNF
            snfp = PopulationSNF(SnifflesConfig.GLOBAL, handle)
            c = 0
            for call in self.svcalls:
                c += 1 if snfp.store(call) else 0
            snfp.write_and_index()

        self.snf_index = snfp.get_index()
        self.snf_total_length = snfp.get_total_length()
        self.snf_candidate_count = c

    def emit(self, vcf_out: VCF = None, **kwargs) -> int:
        res = super().emit(vcf_out=vcf_out, **kwargs)

        if psnf_out := kwargs.get('psnf_out'):
            psnf_out.add_result(self)

        return res


class ErrorResult:
    error = True

    def __init__(self, msg: str):
        self.message = msg

    def __str__(self):
        return self.message

    def emit(self, **kwargs) -> int:
        log.error(f'{self.message}')
        return 0
