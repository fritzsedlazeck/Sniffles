#!/usr/bin/env python3
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     18.06.2025
# Author:      Luis Paulin
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
from pysam import AlignmentFile as Bam
from pysam import FastaFile as Fasta
from pathlib import Path
from sniffles.sv import SVCall
from sniffles.config import SnifflesConfig
import logging
from spoa import poa


GAP = "-"
DEBUG = False
log = logging.getLogger('sniffles.local_asm')


# For SPOA
# tweaked SW:0, NW:1, semi:2
class SPOA(object):
    def __init__(self):
        self.algorithm = {"local": 0, "global": 1, "semig": 2}
        # sv size
        self.sv_small = 400
        self.sv_middl = 1200
        self.sv_large = 10000
        # scores SV small x <= 400
        self.sv_s_match = 10  # m MATCH
        self.sv_s_misma = -10  # n MISMATCH
        self.sv_s_gapop = -20  # g GAP OPEN
        self.sv_s_gapex = 0  # e GAP EXTENSION
        self.sv_s_gap2o = -10  # q GAP OPEN penalty of the second affine function
        self.sv_s_gap2e = -1  # c GAP EXTENSION penalty of the second affine function
        # scores SV mid 400 < X <= 1200
        self.sv_m_match = 13  # m MATCH
        self.sv_m_misma = -5  # n MISMATCH
        self.sv_m_gapop = -25  # g GAP OPEN
        self.sv_m_gapex = 0  # e GAP EXTENSION
        self.sv_m_gap2o = -25  # q GAP OPEN penalty of the second affine function
        self.sv_m_gap2e = -1  # c GAP EXTENSION penalty of the second affine function
        # scores SV large 1200 < X < 10000
        self.sv_l_match = 10  # m MATCH
        self.sv_l_misma = -10  # n MISMATCH
        self.sv_l_gapop = -30  # g GAP OPEN
        self.sv_l_gapex = -1  # e GAP EXTENSION
        self.sv_l_gap2o = -15  # q GAP OPEN penalty of the second affine function
        self.sv_l_gap2e = -1  # c GAP EXTENSION penalty of the second affine function
        # set scores
        self.match = 0
        self.miss = 0
        self.gap_open = 0
        self.gap_expand = 0
        self.gap_open2 = 0
        self.gap_expand2 = 0

    def set(self, svlen: int):
        # score_match, score_mismatch, score_gap_open, score_gap_extension, score_gap2o, score_gap2e
        if abs(svlen) <= self.sv_small:
            a, b, c, d, e, f = (self.sv_s_match, self.sv_s_misma, self.sv_s_gapop, self.sv_s_gapex,
                                self.sv_s_gap2o, self.sv_s_gap2e)
        elif self.sv_small < abs(svlen) <= self.sv_middl:
            a, b, c, d, e, f = (self.sv_m_match, self.sv_m_misma, self.sv_m_gapop, self.sv_m_gapex,
                                self.sv_m_gap2o, self.sv_m_gap2e)
        else:
            a, b, c, d, e, f = (self.sv_l_match, self.sv_l_misma, self.sv_l_gapop, self.sv_l_gapex,
                                self.sv_l_gap2o, self.sv_l_gap2e)
        self.match, self.miss, self.gap_open, self.gap_expand, self.gap_open2, self.gap_expand2 = a, b, c, d, e, f


class LocalAsm(object):
    def __init__(self, sv: SVCall):
        log.debug(f'Local asm in SV with filter {sv.filter}')
        self.sv = sv
        self.info_sv = ""

    def update_sv_cand(self, sv_pos: int, sv_seq: str):
        log.info(self.info_sv)
        if "MOSAIC_VAF" == self.sv.filter:
            self.sv.set_info("MOSAIC", True)
        if "MOSAIC_VAF" == self.sv.filter:
            self.sv.filter = f"GT"
        else:
            self.sv.filter = f"PASS"
        self.sv.qc = True
        self.sv.pos = sv_pos
        self.sv.set_info("LASM", True)
        if "DEL" == self.sv.svtype:
            self.sv.end = sv_pos+len(sv_seq)+1
        elif "INS" == self.sv.svtype:
            self.sv.end = sv_pos+1

    @staticmethod
    def open_bam(bam_file_in: str, ref: str) -> Bam:
        if "cram" in bam_file_in:
            return Bam(bam_file_in, mode=r"rc", require_index=True, reference_filename=ref)
        else:
            return Bam(bam_file_in, mode=r"rb", require_index=True)

    @staticmethod
    def ref_seq_by_region(fas: Fasta, region: str):
        return fas.fetch(region=region)

    @staticmethod
    def region_pos(region: str):
        [_, posit] = region.split(":")
        [start, end] = posit.split("-")
        return int(start), int(end)

    def select_padding(self, read_type: str, default_pad: int = 2000):
        svlen_pad = self.sv.svlen if "sv" == read_type else round(self.sv.svlen/2)
        max_padding = default_pad*3
        if abs(self.sv.svlen) <= 400:
            return min(max(round(svlen_pad*0.5), round(default_pad*0.70)), max_padding)
        elif 400 < abs(self.sv.svlen) <= 1200:
            return min(max(round(svlen_pad*0.5), round(default_pad*1.0)), max_padding)
        else:
            return min(max(round(svlen_pad*0.75), round(default_pad*1.2)), max_padding)

    # Analysis
    def read_seq_by_name_del(self, bam: Bam, region: str, rnames: list[str]):
        # reference and reads start and end at the same position,
        # no gaps should exist in the start/end of the alignment
        contig = ""
        # svlen = self.sv.end - self.sv.pos
        start, stop = [], []
        read_seqs = []
        read_seqs_ref = []
        extra_pad = 100  # during consensus building some bases are lost
        # sv.start  - r.reference_start
        for read in bam.fetch(region=region):
            if "" == contig:
                contig = read.reference_name
            if read.query_name in rnames:
                if read.query_sequence is not None:
                    padding = self.select_padding("sv")
                    pos = self.sv.pos - read.reference_start - padding
                    end = self.sv.end - read.reference_start + padding + 1
                    read_seq_del = read.query_sequence[pos:end]
                    read_len = len(read.query_sequence)
                    if len(read_seq_del) >= 2*padding and (pos > 0 and 0 < end < read_len):
                        read_seqs.append(read_seq_del)
                        start.append(self.sv.pos - padding - extra_pad)
                        stop.append(self.sv.end + padding + extra_pad)

        region = f'{contig}:{min(start)}-{max(stop)}' if len(read_seqs) > 0 else ""
        return read_seqs, read_seqs_ref, region

    def solve_del(self, region: str, sv_aln: str, ref_aln: str, eps: float = 0.15, max_gaps_aln: int = 3):
        svlen = abs(self.sv.svlen)
        [ref_pos, _] = self.region_pos(region)
        # zero_based = -1
        gap_size = 0
        del_pos = 0
        n_gaps = 0
        count_gap = True
        # ref has gap similar size to DEL
        ref_pos_calc = 0
        all_pos = 0
        for char in sv_aln:
            all_pos += 1
            if GAP == char:
                gap_size += 1
                if count_gap:
                    n_gaps += 1
                    count_gap = False
            else:
                count_gap = True
                if abs(gap_size - svlen)/float(svlen) <= eps and gap_size > 0:
                    del_pos = ref_pos+ref_pos_calc
                    if DEBUG or True:
                        log.debug(f'{del_pos},{ref_pos_calc},{gap_size},{abs(gap_size-svlen)/float(svlen)},{eps}')
                    break
                ref_pos_calc += 1
                gap_size = 0
        count_pos = 0
        ref_seq = ""
        for _ in ref_aln:
            count_pos += 1
            if count_pos == ref_pos_calc:
                ref_seq = ref_aln[count_pos:(count_pos+gap_size)]
                break
        return [del_pos, ref_seq, len(ref_seq) > 0 and
                (abs((gap_size - svlen)/float(svlen)) <= eps and gap_size > 0 and n_gaps <= max_gaps_aln)]

    def read_seq_by_name_ins(self, bam: Bam, region: str, rnames: list[str]):
        # reference and reads start and end at the same position,
        # no gaps should exist in the start/end of the alignment
        contig = ""
        start, stop = [], []
        read_seqs_sv = []
        read_seqs_ref = []
        # sv.start  - r.reference_start
        for read in bam.fetch(region=region):
            if "" == contig:
                contig = read.reference_name
            if read.query_name in rnames:
                if read.query_sequence is not None:
                    padding = self.select_padding("sv")
                    pos = max(self.sv.pos - read.reference_start - padding, 0)
                    end = self.sv.pos - read.reference_start + self.sv.svlen + padding
                    read_seq_ins = read.query_sequence[pos:end]
                    read_len = len(read.query_sequence)
                    if len(read_seq_ins) >= self.sv.svlen+(2*padding) and (pos > 0 and 0 < end < read_len):
                        read_seqs_sv.append(read_seq_ins)
                        start.append(self.sv.pos - padding)
                        stop.append(self.sv.pos + self.sv.svlen + padding)

        region = f'{contig}:{min(start)}-{max(stop)}' if len(read_seqs_sv) > 0 else ""
        return read_seqs_sv, read_seqs_ref, region

    def solve_ins(self, region: str, sv_aln: str, ref_aln: str, eps: float = 0.15, max_gaps_aln: int = 3):
        [ref_pos, _] = self.region_pos(region)
        # zero_based = -1
        gap_size = 0
        ins_pos = 0
        n_gaps = 0
        count_gap = True
        # ref has gap similar size to INS
        ref_pos_calc = 0
        all_pos = 0
        for char in ref_aln:
            all_pos += 1
            if GAP == char:
                gap_size += 1
                if count_gap:
                    n_gaps += 1
                    count_gap = False
            else:
                count_gap = True
                if abs((gap_size - self.sv.svlen)/self.sv.svlen) <= eps and gap_size > 0 and ref_pos > 0:
                    ins_pos = ref_pos+ref_pos_calc
                    if DEBUG:
                        log.debug(f'{ins_pos}, {ref_pos_calc}, {gap_size}, '
                                  f'{abs(float(gap_size - self.sv.svlen) / self.sv.svlen)}, {eps}')
                    break
                ref_pos_calc += 1
                gap_size = 0
        count_pos = 0
        ins_seq = ""
        for _ in sv_aln:
            count_pos += 1
            if count_pos == ref_pos_calc:
                ins_seq = sv_aln[count_pos:(count_pos+gap_size)]
                break
        return [ins_pos, ins_seq, len(ins_seq) > 0 and
                (abs((gap_size - self.sv.svlen)/self.sv.svlen) <= eps and gap_size > 0 and n_gaps <= max_gaps_aln)]

    def assembly(self, config: SnifflesConfig, min_reads: int = 5, max_reads: int = 30, fetch_read: int = 5000) -> bool:
        file_path = Path(config.input)
        if file_path.suffix not in [".bam", ".cram"]:
            log.warning(f'Input {config.input} is not BAM or CRAM, skipping local assembly')
            return False
        if config.reference is None:
            log.warning(f'Reference not provided, skipping local assembly')
            return False
        bam = self.open_bam(config.input, config.reference)
        fas = Fasta(config.reference)
        spoa = SPOA()
        region_sv = f'{self.sv.contig}:{self.sv.pos}-{self.sv.end}'
        region_sv_fetch = f'{self.sv.contig}:{self.sv.pos-fetch_read}-{self.sv.end+fetch_read}'
        svlen = abs(self.sv.svlen)
        self.info_sv = f'{self.sv.id},{region_sv},{svlen},NA,NA,{self.sv.filter},init'
        sv_pos = 0
        sv_seq = ""  # inserted for INS, deleted-ref for DEL
        spoa.set(svlen)
        if "DEL" == self.sv.svtype:
            read_seq, read_seq_ref, region_fas = self.read_seq_by_name_del(bam, region_sv_fetch, self.sv.rnames)
        elif "INS" == self.sv.svtype:
            read_seq, read_seq_ref, region_fas = self.read_seq_by_name_ins(bam, region_sv_fetch, self.sv.rnames)
        else:
            # already filtered, should not happen
            log.warning(f'[WARN] SVTYPE: {self.sv.svtype} not supported, skipping local assembly')
            return False
        # get fasta
        report_sv = False  # change to True if the sequence can be solved, i.e. SVLEN constrain is met
        if len(read_seq) >= min_reads:
            ref_seq = self.ref_seq_by_region(fas, region_fas)
            # default
            min_cov = round(len(read_seq)*0.50)
            read_seq = [read_seq[i] for i in range(0, max_reads)] if len(read_seq) > max_reads else read_seq
            consensus, _ = poa(sequences=read_seq, algorithm=spoa.algorithm["local"], genmsa=True, min_coverage=min_cov)
            # alt = self.sv.alts[0] if self.sv.alts is not None else "N"
            _, msa = poa(sequences=[consensus, ref_seq], algorithm=spoa.algorithm["local"], genmsa=True, min_coverage=1,
                         m=spoa.match, n=spoa.miss, g=spoa.gap_open, e=spoa.gap_expand,
                         q=spoa.gap_open2, c=spoa.gap_expand2)
            if len(msa) == 2:
                [sv_aln, ref_aln] = msa
            else:
                return report_sv
            if "INS" == self.sv.svtype:
                [sv_pos, sv_seq, report_sv] = self.solve_ins(region_fas, sv_aln, ref_aln)
            elif "DEL" == self.sv.svtype:
                [sv_pos, sv_seq, report_sv] = self.solve_del(region_fas, sv_aln, ref_aln)
            self.info_sv = f'{self.sv.id},{self.sv.contig}:{self.sv.pos},{svlen},{len(read_seq)},{self.sv.filter},done'
        # return info_sv, sv_pos, report_sv, sv_seq, msa
        if report_sv:
            self.update_sv_cand(sv_pos, sv_seq)
        return report_sv
