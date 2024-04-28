#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     14.01.2024
# Author:      Hermann Romanek
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
"""
Resource monitoring utilities (currently memory) for Sniffles
"""
import logging
import os
import time
from dataclasses import dataclass
from threading import Thread

from sniffles.config import SnifflesConfig


log = logging.getLogger('sniffles.ResourceMonitor')

@dataclass
class MemoryUsage:
    """
    A snapshot of memory usage for this process and all its workers
    """
    pid: int
    mem: float  # RSS in KiB


class ResourceMonitor:
    _pid: int
    _parent: float
    _unknown: float
    _workers: dict[int, MemoryUsage]
    _filename = None

    WORKER_MEMORY_LIMIT = 1024 * 1024 * 2  # in KiB

    def __init__(self, config: SnifflesConfig):

        self._slots = config.threads
        self._update_interval = config.dev_monitor_memory or 30.0
        self._pid = os.getpid()

        try:
            from psutil import Process
        except ImportError:
            logging.getLogger('sniffles.dependencies').warning(f'Dependency psutil not available - resource monitoring is disabled.')
            self._process = None
            self._workers = {}
        else:
            self._process = Process(self._pid)

            self._running = True

            if config.dev_monitor_memory:
                self._filename = filename = f'memory-{config.run_id}.csv'
                try:
                    self._file = open(filename, "w")
                    self._file.write(self._generate_header())
                except:
                    log.exception(f'Unable to write memory usage to {filename}')
                else:
                    log.info(f'Logging memory usage to {filename}')
            else:
                self._file = None

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
        """
        Stop memory monitoring thread. This may take up to update_interval to complete.
        """
        self._running = False

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

    def __call__(self, worker_id: int, worker_pid: int) -> bool:
        """
        Call the resource monitor to get a recycle hint. Returns True if the given process should recycle.
        """
        if worker_pid in self._workers:
            return self._workers[worker_pid].mem > self.WORKER_MEMORY_LIMIT

        return False
