"""
Resource monitoring utilities (currently memory) for Sniffles
"""
import os
import time
from dataclasses import dataclass
from threading import Thread

import psutil


@dataclass
class MemoryUsage:
    """
    A snapshot of memory usage for this process and all its workers
    """
    parent: float
    workers: dict[int, float]


class ResourceMonitor:
    _data: MemoryUsage = None
    _pid: int

    def __init__(self, update_interval: float = 10.0):
        self._update_interval = update_interval
        self._pid = os.getpid()
        self._process = psutil.Process(self._pid)
        self._running = True
        self.update()
        self._update_thread = Thread(target=self._run_update, daemon=True)
        self._update_thread.start()

    @property
    def parent(self) -> float:
        return self._data.parent

    @property
    def workers(self) -> dict:
        return self._data.workers

    @property
    def total(self) -> float:
        return self._data.parent + sum(self._data.workers.values())

    def stop(self):
        self._running = False

    def update(self):
        self._data = MemoryUsage(parent=self._process.memory_info().rss, workers={
            child.pid: child.memory_info().rss for child in self._process.children()
        })
        print(f'Memory usage: {self.total / 1024:.2f} KiB')

    def _run_update(self):
        while self._running:
            self.update()
            time.sleep(self._update_interval)
