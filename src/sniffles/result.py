import os
import pickle
import sys

from sniffles.sv import SVCall
from sniffles.vcf import VCF


class Result:
    """
    A generic result of a sniffles task executed by a worker process.
    """
    processed_read_count: int
    task_id: int
    processed_read_count: int
    svcalls: list[SVCall]
    svcount: int
    error: bool = False

    def __init__(self, task: 'Task', svcalls: list[SVCall], candidates_processed: int):
        self.task_id = task.id
        self.processed_read_count = candidates_processed
        self.svcount = len(svcalls)
        self.store_calls(svcalls)

    def store_calls(self, svcalls):
        self.svcalls = svcalls

    def cleanup(self):
        """
        Optional clean up code after writing this result
        """


class CallResult(Result):
    coverage_average_total: float
    has_snf = False
    snf_filename = None
    snf_index = None
    snf_total_length = None
    snf_candidate_count = None


class GenotypeResult(Result):
    """
    Result of a genotyping run
    """


class CombineResult(Result):
    """
    Result of a combine run for one task, simple variant with calls in memory. Must be pickleable.
    """
    def emit(self, file: VCF):
        """
        Emit this result to the given file
        """
        for call in self.svcalls:
            file.write_call(call)


class CombineResultTmpFile(CombineResult):
    """
    Result of a combine run, with calls in a temporary file instead of memory
    """
    @property
    def tmpfile_name(self) -> str:
        return f'tmp{self.task_id}.part'

    def store_calls(self, svcalls):
        with open(self.tmpfile_name, 'wb') as f:
            f.write(pickle.dumps(svcalls))

    @property
    def svcalls(self) -> list[SVCall]:
        with open(self.tmpfile_name, 'rb') as f:
            return pickle.loads(f.read())

    def cleanup(self):
        os.unlink(self.tmpfile_name)


class ErrorResult:
    error = True

    def __init__(self, ex: Exception):
        self.message = f'{ex}'

    def __str__(self):
        return self.message
