#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     01.07.2024
# Author:      Hermann Romanek
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
import logging
from dataclasses import dataclass
from typing import Optional

log = logging.getLogger(__name__)


@dataclass
class Region:
    """
    A region from a bed file.

    Coordinates are 0-based, half open
    """
    contig: str
    start: int
    end: int

    @classmethod
    def from_bed_line(cls, line: str) -> Optional['Region']:
        try:
            contig, start, end = line.split('\t')[:3]
            return Region(
                contig=contig,
                start=int(start),
                end=int(end),
            )
        except ValueError as ex:
            log.warning(f'Invalid region line: {line} ({ex})')
            return None

    @classmethod
    def from_string(cls, value: str) -> Optional['Region']:
        try:
            contig, pos = value.split(':')
            start, end = map(int, pos.split('-'))
            return Region(
                contig=contig,
                start=start,
                end=end,
            )
        except ValueError as ex:
            log.warning(f'Invalid region string: {value} ({ex})')
            return None

    def __str__(self) -> str:
        return f'{self.contig}:{self.start}-{self.end}'
