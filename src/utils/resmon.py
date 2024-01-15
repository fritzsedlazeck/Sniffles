"""
Resource monitoring utilities (currently memory) for Sniffles
"""
import logging
import os
import time
from dataclasses import dataclass
from threading import Thread


@dataclass
class MemoryUsage:
    """
    A snapshot of memory usage for this process and all its workers
    """
    pid: int
    mem: float


class ResourceMonitor:
    _pid: int
    _parent: float
    _unknown: float
    _workers: dict[int, MemoryUsage]

    def __init__(self, worker_slots: int, update_interval: float = 10.0, filename: str = None):
        self._slots = worker_slots
        self._update_interval = update_interval
        self._pid = os.getpid()

        try:
            from psutil import Process
        except ImportError as ex:
            logging.getLogger('sniffles.dependencies').warning(f'Dependency psutil not available - resource monitoring is disabled.')
            self._process = None
        else:
            self._process = Process(self._pid)

            self._running = True

            if filename is None:
                filename = f'memory-{os.environ.get("SLURM_JOB_ID") or os.getpid()}.csv'
            if filename is not None:
                self._file = open(filename, "w")
                self._file.write(self._generate_header())
            else:
                self._file = None

            self._filename = filename

            self.update()
            self._update_thread = Thread(target=self._run_update, daemon=True)
            self._update_thread.start()

    @property
    def parent(self) -> float:
        return self._parent

    @property
    def workers(self) -> dict[int, float]:
        """
        dict pid -> memory
        """
        return {
            obj.pid: obj.mem for obj in self._workers.values() if obj is not None
        }

    @property
    def total(self) -> float:
        return self.parent + sum(self.workers.values())

    def stop(self):
        self._running = False

    def set_slot_process(self, slot: int, process: 'SnifflesWorker'):
        self._slots[slot] = process

    def _generate_header(self):
        return f'''total,pidp,memp,{','.join(f'pid{i},mem{i}' for i in range(self._slots))}\n'''

    def update(self):
        self._parent = self._process.memory_info().rss / 1024.0
        self._workers = {
            child.pid: MemoryUsage(child.pid, child.memory_info().rss / 1024.0) for child in self._process.children()
        }

        if self._file:
            self._file.write(f'''{','.join(self.generate())}\n''')
            self._file.flush()

        # print(f'Memory usage: {self.total / 1024:.2f} KiB')

    def _run_update(self):
        while self._running:
            self.update()
            time.sleep(self._update_interval)

        if self._file:
            self._file.close()

    def generate(self):
        yield str(self.total)
        yield str(self._pid)
        yield str(self.parent)
        for k, v in self.workers.items():
            yield str(k)
            yield str(v)

    @property
    def filename(self) -> str:
        return self._filename
