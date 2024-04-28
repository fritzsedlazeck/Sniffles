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

    def store_calls(self, svcalls):
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
    def __str__(self):
        return f'CombineResult #{self.task_id}'


class CombineResultTmpFile(CombineResult):
    """
    Result of a combine run, with calls in a temporary file instead of memory
    """
    @property
    def tmpfile_name(self) -> str:
        return f'result-{self.run_id}-{self.task_id}.part'

    def store_calls(self, svcalls):
        with open(self.tmpfile_name, 'wb') as f:
            f.write(pickle.dumps(svcalls))

    @property
    def svcalls(self) -> list[SVCall]:
        with open(self.tmpfile_name, 'rb') as f:
            return pickle.loads(f.read())

    def emit(self, **kwargs) -> int:
        res = super().emit(**kwargs)
        self.cleanup()
        return res

    def cleanup(self):
        os.unlink(self.tmpfile_name)


class ErrorResult:
    error = True

    def __init__(self, ex: Exception):
        self.message = f'{ex}'

    def __str__(self):
        return self.message
