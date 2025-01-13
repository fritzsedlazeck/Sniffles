#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     26.09.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
import logging
import os
import pickle
import json
import gzip
import math
from functools import cached_property
from typing import Optional

from sniffles import sv
from sniffles.config import SnifflesConfig


log = logging.getLogger(__name__)


class SNFileBase:
    """
    Base class for SNF-Files
    """
    header_length: int
    _header: Optional[dict]

    def __init__(self, config: SnifflesConfig, handle, filename=None):
        self.config = config
        self.handle = handle
        self.filename = filename
        self.blocks = {}
        self._header = None
        self._index = {}
        self.total_length = 0
        self._results = []

    @classmethod
    def open(cls, filename: str) -> 'SNFileBase':
        obj = cls(SnifflesConfig.GLOBAL, open(filename, 'rb'), filename)
        obj.read_header()
        return obj

    @property
    def index(self) -> dict:
        return self._index

    @property
    def header(self) -> dict:
        return self._header

    @cached_property
    def population(self) -> object | None:
        """
        Is this a population SNF?
        """
        return self.header.get('population', None)

    @cached_property
    def reqc(self) -> bool:
        """
        Was this file created by an old Version of sniffles we want to redo qc?
        """
        if self.config.reqc == 'auto':
            try:
                build, _, _ = self.header['config']['build'].partition('-')
            except (KeyError, AttributeError):
                log.warning(f'Unable to determine version of SNF file {self.filename} for auto-reqc')
                return True
            else:
                return build < '2.5.3'
        else:
            return self.config.reqc

    def is_open(self) -> bool:
        return self.handle is not False

    def _open(self):
        if self.handle is not False:
            self.close()
        self.handle = open(self.filename, "rb")

    def store(self, svcand):
        block_index = int(svcand.pos / self.config.snf_block_size) * self.config.snf_block_size
        if block_index not in self.blocks:
            self.blocks[block_index] = {svtype: [] for svtype in sv.TYPES}
            self.blocks[block_index]["_COVERAGE"] = {}
        if not self.config.output_rnames:
            svcand.rnames = None
        self.blocks[block_index][svcand.svtype].append(svcand)

    def serialize_block(self, block_id):
        return pickle.dumps(self.blocks[block_id])

    def unserialize_block(self, data: bytes):
        return pickle.loads(data)

    def write_and_index(self):
        if not self.is_open():
            self._open()
        offset = 0
        for block_id in sorted(self.blocks):
            data = gzip.compress(self.serialize_block(block_id))
            self.handle.write(data)
            data_len = len(data)
            self._index[block_id] = (offset, data_len)
            offset += data_len
            self.total_length += data_len
        if self.config.combine_close_handles:
            self.close()

    def read_header(self):
        if not self.is_open():
            self._open()
        try:
            header_text = self.handle.readline()
            self.header_length = len(header_text)
            self._header = json.loads(header_text.strip())
        except Exception as e:
            print(f"Error when reading SNF header from '{self.filename}': {e}. The file may not be a valid .snf file or could have been corrupted.")
            raise e
        self._index = self._header["index"]
        if self.config.combine_close_handles:
            self.close()

    def read_blocks(self, contig, block_index):
        if not self.is_open():
            self._open()
        block_index = str(block_index)
        if contig not in self.index:
            if self.config.combine_close_handles:
                self.close()
            return None

        if block_index not in self.index[contig]:
            if self.config.combine_close_handles:
                self.close()
            return None

        blocks = []
        for block_data_start, block_data_length in self.index[contig][block_index]:
            try:
                self.handle.seek(self.header_length + block_data_start)
                data = gzip.decompress(self.handle.read(block_data_length))
                blocks.append(self.unserialize_block(data))
            except Exception as e:
                print(f"Error when reading block '{contig}.{block_index}' from '{self.filename}': {e}. The file may not be a valid .snf file or could have been corrupted.")
                if self.config.combine_close_handles:
                    self.close()
                raise e
        if self.config.combine_close_handles:
            self.close()
        return blocks

    def get_index(self):
        return self.index

    def get_total_length(self):
        return self.total_length

    def add_result(self, result):
        if result.has_snf:
            self._results.append(result)

    def _calculate_contig_coverages(self, contigs: list[str]) -> dict[str, float]:
        contig_coverages: dict[str, list[float]] = {
            c: [] for c in contigs
        }

        for r in self._results:
            contig_coverages[r.contig].append(r.coverage_average_total)

        res = {}
        for contig in contig_coverages:
            res[contig] = sum(contig_coverages[contig]) / len(contig_coverages[contig]) if len(contig_coverages[contig]) > 0 else 0

        return res

    def _create_header(self, config: SnifflesConfig, main_index: dict, snf_candidate_count: int) -> dict:
        """
        Generate the header for this SNF file.
        """
        return {"config": config.__dict__, "index": main_index, "snf_candidate_count": snf_candidate_count}

    def write_results(self, config: SnifflesConfig, contigs: list[str]) -> int:
        """
        Writes all added results (regional temporary .snf files) to this file. Returns SNF candidate count
        """
        main_index = {}
        offset = 0
        snf_candidate_count = sum(r.snf_candidate_count for r in self._results)
        parts_sorted = sorted(self._results, key=lambda r: r.task_id)
        for part in parts_sorted:
            part_contig = part.contig
            if part_contig not in main_index:
                main_index[part_contig] = {}
            for block, (part_block_start, part_block_len) in part.snf_index.items():
                if block not in main_index[part_contig]:
                    main_index[part_contig][block] = []
                main_index[part_contig][block].append((part_block_start + offset, part_block_len))
            offset += part.snf_total_length

        config.contig_coverages = self._calculate_contig_coverages(contigs)
        header = self._create_header(config, main_index, snf_candidate_count)
        header_json = json.dumps(header, default=lambda obj: "<Unstored_Object>") + "\n"
        self.handle.write(header_json.encode())

        for part in parts_sorted:
            with open(part.snf_filename, "rb") as part_handle:
                part_data = part_handle.read()
                part_handle.close()
            self.handle.write(part_data)
            os.remove(part.snf_filename)

        return snf_candidate_count

    def close(self) -> None:
        if self.handle is not False:
            self.handle.close()
            self.handle = False

    def get_all_blocks(self, contig: str) -> dict[str, dict]:
        """
        Fetches all blocks for given contig.
        """
        blocks = {}
        if contig in self.index:
            for block_start in self.index[contig].keys():
                blocks[block_start] = self.read_blocks(contig, block_start)[0]
        return blocks

    def get_full_coverage(self, contig: str):
        coverage = {}
        for b in self.get_all_blocks(contig).values():
            coverage.update(b['_COVERAGE'])
        return coverage


class SNFile(SNFileBase):
    def annotate_block_coverages(self, lead_provider, resolution=500):
        config = self.config
        start_bin = lead_provider.covrtab_min_bin
        end_bin = int(lead_provider.end / config.coverage_binsize) * config.coverage_binsize
        coverage_fwd = 0
        coverage_rev = 0

        coverage_sum = 0
        bin_count = 0

        coverage_binsize_combine = self.config.coverage_binsize_combine
        snf_block_size = config.snf_block_size

        for bin in range(start_bin, end_bin + config.coverage_binsize, config.coverage_binsize):
            if bin in lead_provider.covrtab_fwd:
                coverage_fwd += lead_provider.covrtab_fwd[bin]

            if bin in lead_provider.covrtab_rev:
                coverage_rev += lead_provider.covrtab_rev[bin]

            coverage_sum += coverage_fwd + coverage_rev
            bin_count += 1

            if bin % coverage_binsize_combine == 0:
                block_index = int(bin / snf_block_size) * snf_block_size

                coverage_total_curr = math.ceil(coverage_sum / float(bin_count))
                if coverage_total_curr > 0:
                    if block_index not in self.blocks:
                        self.blocks[block_index] = {svtype: [] for svtype in sv.TYPES}
                        self.blocks[block_index]["_COVERAGE"] = {}

                    self.blocks[block_index]["_COVERAGE"][bin] = coverage_total_curr

                coverage_sum = 0
                bin_count = 0


class RemoteIndexSNFile(SNFile):
    """
    An SNFile, but its header data will be provided by a remote source.
    """

    @property
    def index(self) -> dict:
        if self._header is None:
            self.read_header()

        return super().index

    @property
    def header(self) -> dict:
        if self._header is None:
            self.read_header()

        return super().header
