#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     28.05.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
from dataclasses import dataclass
import re
import itertools
from typing import Optional

import pysam

# for: --dev-cache
import os
import sys
import pickle
# end: for: --dev-cache

from sniffles import util
from sniffles import sv
from sniffles.region import Region


@dataclass
class Lead:
    read_id: int = None  # or tuple[int, str, str] for phased reads, with (read_id, HP, PS)
    read_qname: str = None
    contig: str = None
    ref_start: int = None
    ref_end: int = None
    qry_start: int = None
    qry_end: int = None
    strand: str = None
    mapq: int = None
    nm: float = None
    source: str = None
    svtype: str = None
    svlen: Optional[int] = None
    seq: Optional[str] = None
    svtypes_starts_lens: list = None


def CIGAR_analyze(cigar):
    buf = ""
    readspan = 0
    refspan = 0
    clip_start = None
    clip = 0
    for c in cigar:
        if c.isnumeric():
            buf += c
        else:
            oplen = int(buf)
            h = False
            if c in "MIX=":
                readspan += oplen
                h = True
            if c in "MDX=N":
                refspan += oplen
                h = True
            if not h:
                if c in "SH":
                    if clip_start == None and readspan + refspan > 0:
                        clip_start = clip
                    clip += oplen
                else:
                    raise f"Unknown CIGAR operation: '{c}'"
            buf = ""
    if clip_start == None:
        clip_start = clip
    return clip_start, clip - clip_start, refspan, readspan


def CIGAR_analyze_regex(cigar):
    # TODO: Obsolete
    opsums = {"M": 0, "I": 0, "D": 0, "=": 0, "X": 0, "N": 0, "H": 0, "S": 0}
    iter = re.split(r"(\d+)", cigar)
    for i in range(1, len(iter) - 1, 2):
        op = iter[i + 1]
        if op != "H" and op != "S":
            readstart_fwd = opsums["H"] + opsums["S"]
        opsums[iter[i + 1]] += int(iter[i])
    readstart_rev = opsums["H"] + opsums["S"] - readstart_fwd
    refspan = opsums["M"] + opsums["D"] + opsums["="] + opsums["X"] + opsums["N"]
    readspan = opsums["M"] + opsums["I"] + opsums["="] + opsums["X"]
    return readstart_fwd, readstart_rev, refspan, readspan


def CIGAR_tolist(cigar):
    # TODO: Obsolete (see CIGAR_tolist_analyze)
    """
    CIGAR string : str -> List of CIGAR operation & length tuples : [(op1:char, op1_length:int),...]
    """
    buf = ""
    ops = []
    for c in cigar:
        if c.isnumeric():
            buf += c
        else:
            ops.append((c, int(buf)))
            buf = ""
    return ops


def CIGAR_listrefspan(ops):
    # TODO: Obsolete (see CIGAR_analyze)
    # TODO(Potential): Detect&utilize minimap2 condensed supplementary alignment cigar for speed
    return sum(oplen for op, oplen in ops if op == "M" or op == "D" or op == "=" or op == "X" or op == "N")


def CIGAR_listreadspan(ops):
    # TODO: Obsolete (see CIGAR_analyze)
    return sum(oplen for op, oplen in ops if op == "M" or op == "I" or op == "=" or op == "X")


def CIGAR_listreadstart_fwd(ops):
    # TODO: Obsolete (see CIGAR_analyze)
    """
    Position in query read where CIGAR alignment starts (i.e. taking into account start clipping)
    """
    op, oplen = ops[0]
    op2, op2len = ops[1]
    if op == "H" or op == "S":
        assert (op2 != "H" and op2 != "S")
        return oplen
    else:
        return 0


def CIGAR_listreadstart_rev(ops):
    # TODO: Obsolete (see CIGAR_analyze)
    """
    Position in query read where CIGAR alignment starts (i.e. taking into account start clipping)
    """
    op, oplen = ops[-1]
    op2, op2len = ops[-2]
    if op == "H" or op == "S":
        assert (op2 != "H" and op2 != "S")
        return oplen
    else:
        return 0


OPTAB = {pysam.CMATCH: (1, 1, 0),
         pysam.CEQUAL: (1, 1, 0),
         pysam.CDIFF: (1, 1, 0),
         pysam.CINS: (1, 0, 1),
         pysam.CDEL: (0, 1, 1),
         pysam.CREF_SKIP: (0, 1, 0),
         pysam.CSOFT_CLIP: (1, 0, 1),
         pysam.CHARD_CLIP: (0, 0, 0),
         pysam.CPAD: (0, 0, 0)}
#      pysam.CBACK:      (0,0,0)}

OPLIST = [(0, 0, 0) for i in range(max(int(k) for k in OPTAB.keys()) + 1)]
for k, v in OPTAB.items():
    OPLIST[int(k)] = v


def read_iterindels(read_id, read, contig, config, use_clips, read_nm):
    minsvlen = config.minsvlen_screen
    longinslen = config.long_ins_length / 2.0
    seq_cache_maxlen = config.dev_seq_cache_maxlen
    qname = read.query_name
    mapq = read.mapping_quality
    strand = "-" if read.is_reverse else "+"
    CINS = pysam.CINS
    CDEL = pysam.CDEL
    CSOFT_CLIP = pysam.CSOFT_CLIP

    pos_read = 0
    pos_ref = read.reference_start

    for op, oplength in read.cigartuples:
        add_read, add_ref, event = OPLIST[op]
        if event and oplength >= minsvlen:
            if op == CINS:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read + oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           oplength,
                           seq=read.query_sequence[pos_read:pos_read + oplength] if oplength <= seq_cache_maxlen else None)
            elif op == CDEL:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref + oplength,
                           pos_ref,
                           pos_read,
                           pos_read,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "DEL",
                           -oplength)
            elif use_clips and op == CSOFT_CLIP and oplength >= longinslen:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read + oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           None,
                           seq=None)
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength


def get_cigar_indels(read_id, read, contig, config, use_clips, read_nm):
    minsvlen = config.minsvlen_screen
    longinslen = config.long_ins_length / 2.0
    seq_cache_maxlen = config.dev_seq_cache_maxlen
    qname = read.query_name
    mapq = read.mapping_quality
    strand = "-" if read.is_reverse else "+"
    CINS = pysam.CINS
    CDEL = pysam.CDEL
    CSOFT_CLIP = pysam.CSOFT_CLIP

    INS_SUM = 0
    DEL_SUM = 0

    pos_read = 0
    pos_ref = read.reference_start
    for op, oplength in read.cigartuples:
        add_read, add_ref, event = OPLIST[op]
        if event:
            if op == CINS:
                INS_SUM += oplength
            elif op == CDEL:
                DEL_SUM += oplength
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength

    return INS_SUM, DEL_SUM


def read_itersplits_bnd(read_id, read, contig, config, read_nm):
    assert (read.is_supplementary)
    # SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads = []
    supps = [part.split(",") for part in read.get_tag("SA").split(";") if len(part) > 0]

    if len(supps) > config.max_splits_base + config.max_splits_kb * (read.query_length / 1000.0):
        return

    if read.is_reverse:
        qry_start = read.query_length - read.query_alignment_end
    else:
        qry_start = read.query_alignment_start

    curr_lead = Lead(read_id,
                     read.query_name,
                     contig,
                     read.reference_start,
                     read.reference_start + read.reference_length,
                     qry_start,
                     qry_start + read.query_alignment_length,
                     "-" if read.is_reverse else "+",
                     read.mapping_quality,
                     read_nm,
                     "SPLIT_SUP",
                     "?")
    all_leads.append(curr_lead)

    prim_refname, prim_pos, prim_strand, prim_cigar, prim_mapq, prim_nm = supps[0]
    if prim_refname == contig:
        # Primary alignment is on this chromosome, no need to parse the supplementary
        return

    minpos_curr_chr = min(itertools.chain([read.reference_start], (int(pos) for refname, pos, strand, cigar, mapq, nm in supps if refname == contig)))
    if minpos_curr_chr < read.reference_start:
        # Only process splits once per chr (there may be multiple supplementary alignments on the same chr)
        return

    for refname, pos, strand, cigar, mapq, nm in supps:
        mapq = int(mapq)
        nm = int(nm)
        # if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev = (strand == "-")

        try:
            readstart_fwd, readstart_rev, refspan, readspan = CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero = int(pos) - 1
        split_qry_start = readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start + readspan,
                              strand,
                              mapq,
                              read_nm,
                              "SPLIT_SUP",
                              "?"))

    sv.classify_splits(read, all_leads, config, contig)

    for lead in all_leads:
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            if svtype == "BND":
                bnd = Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=config.bnd_cluster_length,
                           seq=None)
                bnd.bnd_info = arg
                # print(lead.contig,svstart,bnd.bnd_info)
                yield bnd


def read_itersplits(read_id, read, contig, config, read_nm):
    # SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads = []
    supps = [part.split(",") for part in read.get_tag("SA").split(";") if len(part) > 0]
    trace_read = config.dev_trace_read != False and config.dev_trace_read == read.query_name

    if len(supps) > config.max_splits_base + config.max_splits_kb * (read.query_length / 1000.0):
        return

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] passed max_splits check")

    # QC on: 18Aug21, HG002.ont.chr22; O.K.
    # cigarl=CIGAR_tolist(read.cigarstring)
    # if read.is_reverse:
    #    cigarl.reverse()

    # if read.is_reverse:
    #    assert(read.query_length-read.query_alignment_end == CIGAR_listreadstart(cigarl))
    # else:
    #    assert(read.query_alignment_start == CIGAR_listreadstart(cigarl))

    # assert(CIGAR_listrefspan(cigarl)==read.reference_length)
    # assert(CIGAR_listreadspan(cigarl)==read.query_alignment_length)
    # End QC

    if read.is_reverse:
        qry_start = read.query_length - read.query_alignment_end
    else:
        qry_start = read.query_alignment_start

    curr_lead = Lead(read_id,
                     read.query_name,
                     contig,
                     read.reference_start,
                     read.reference_start + read.reference_length,
                     qry_start,
                     qry_start + read.query_alignment_length,
                     "-" if read.is_reverse else "+",
                     read.mapping_quality,
                     read_nm,
                     "SPLIT_PRIM",
                     "?")
    all_leads.append(curr_lead)

    # QC on: 18Aug21; O.K.
    # assert(read.reference_length == CIGAR_listrefspan(CIGAR_tolist(read.cigarstring)))
    # assert(read.query_alignment_start == CIGAR_listreadstart(CIGAR_tolist(read.cigarstring)))
    # assert(read.query_alignment_length == CIGAR_listreadspan(CIGAR_tolist(read.cigarstring)))
    # End QC

    for refname, pos, strand, cigar, mapq, nm in supps:
        mapq = int(mapq)
        nm = int(nm)
        # if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev = (strand == "-")

        try:
            readstart_fwd, readstart_rev, refspan, readspan = CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero = int(pos) - 1
        split_qry_start = readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start + readspan,
                              strand,
                              mapq,
                              read_nm,
                              "SPLIT_SUP",
                              "?"))

        # QC on: 08Sep21; O.K.
        # cigarl=CIGAR_tolist(cigar)
        # assert(CIGAR_listrefspan(cigarl)==refspan)
        # assert(CIGAR_listreadspan(cigarl)==readspan)
        # assert(CIGAR_listreadstart_fwd(cigarl)==readstart_fwd)
        # assert(CIGAR_listreadstart_rev(cigarl)==readstart_rev)
        # End QC

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] all_leads: {all_leads}")

    all_leads = sv.classify_splits(read, all_leads, config, contig)

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] classify_splits(all_leads): {all_leads}")

    """
    if config.dev_trace_read != False:
        print(read.query_name)
        if read.query_name == config.dev_trace_read:
            for lead_i, lead in enumerate(all_leads):
                for svtype, svstart, arg in lead.svtypes_starts_lens:
                    min_mapq=min(lead.mapq,all_leads[max(0,lead_i-1)].mapq)
                    keep=True
                    if not config.dev_keep_lowqual_splits and min_mapq < config.mapq:
                        keep=False
                    print(f"[DEV_TRACE_READ] [REPORT_LEAD_SPLIT] Splits identified from read {read.read_id}")
    """

    for lead_i, lead in enumerate(all_leads):
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            min_mapq = min(lead.mapq, all_leads[max(0, lead_i - 1)].mapq)
            if not config.dev_keep_lowqual_splits and min_mapq < config.mapq:
                continue

            if svtype == "BND":
                bnd = Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=config.bnd_cluster_length,
                           seq=None)
                bnd.bnd_info = arg
                yield bnd

            elif svtype != "NOSV":
                svlen = arg

                yield Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart + svlen if svlen != None and svtype != "INS" else svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=svlen,
                           seq=lead.seq if svtype == "INS" else None)


class LeadProvider:
    def __init__(self, config, read_id_offset):
        self.config = config

        self.leadtab = {}
        self.leadcounts = {}

        for svtype in sv.TYPES:
            self.leadtab[svtype] = {}
            self.leadcounts[svtype] = 0

        self.covrtab_fwd = {}
        self.covrtab_rev = {}
        self.covrtab_min_bin = None
        # self.covrtab_read_start={}
        # self.covrtab_read_end={}

        self.read_id = read_id_offset
        self.read_count = 0

        self.contig = None
        self.start = None
        self.end = None

    def record_lead(self, ld, pos_leadtab):
        leadtab_svtype = self.leadtab[ld.svtype]
        if pos_leadtab in leadtab_svtype:
            leadtab_svtype[pos_leadtab].append(ld)
            lead_count = len(leadtab_svtype[pos_leadtab])
            if lead_count > self.config.consensus_max_reads_bin:
                ld.seq = None
        else:
            leadtab_svtype[pos_leadtab] = [ld]
            lead_count = 1
        self.leadcounts[ld.svtype] += 1

    def build_leadtab(self, regions: list[Region], bam):

        assert (self.contig is None)
        assert (self.start is None)
        assert (self.end is None)

        externals = []
        ld_binsize = self.config.cluster_binsize

        for region in regions:
            self.contig = region.contig
            self.start = region.start if self.start is None else min(region.start, self.start)
            self.end = region.end if self.end is None else max(region.start, self.end)
            self.covrtab_min_bin = int(self.start / self.config.coverage_binsize) * self.config.coverage_binsize

            for ld in self.iter_region(bam, region):
                ld_contig, ld_ref_start = ld.contig, ld.ref_start

                if region.contig == ld_contig and region.start <= ld_ref_start < region.end:
                    pos_leadtab = int(ld_ref_start / ld_binsize) * ld_binsize
                    self.record_lead(ld, pos_leadtab)
                else:
                    externals.append(ld)

        return externals

    def iter_region(self, bam, region: Region):
        leads_all = []
        binsize = self.config.cluster_binsize
        coverage_binsize = self.config.coverage_binsize
        coverage_shift_bins = self.config.coverage_shift_bins
        coverage_shift_min_aln_len = self.config.coverage_shift_bins_min_aln_length
        long_ins_threshold = self.config.long_ins_length * 0.5
        qc_nm = self.config.qc_nm_measure
        phase = self.config.phase
        advanced_tags = qc_nm or phase
        mapq_min = self.config.mapq
        alen_min = self.config.min_alignment_length
        nm_sum = 0
        nm_count = 0
        trace_read = self.config.dev_trace_read

        read_mq20 = 0
        for read in bam.fetch(region.contig, region.start, region.end, until_eof=False):
            if trace_read is not False:
                if trace_read == read.query_name:
                    print(f"[DEV_TRACE_READ] [0b/4] [LeadProvider.iter_region] [{region}] [{read.query_name}] has been fetched and is entering pre-filtering")

            # if self.read_count % 1000000 == 0:
            #    gc.collect()
            if read.reference_start < region.start or read.reference_start >= region.end:
                continue

            self.read_id += 1
            self.read_count += 1

            alen = read.query_alignment_length
            if read.mapping_quality < mapq_min or read.is_secondary or alen < alen_min:
                continue

            has_sa = read.has_tag("SA")
            use_clips = self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm = -1
            curr_read_id = self.read_id
            if advanced_tags:
                if qc_nm:
                    if read.has_tag("NM"):
                        nm_raw = read.get_tag("NM")
                        nm_ratio = read.get_tag("NM") / float(read.query_alignment_length + 1)
                        ins_sum, del_sum = get_cigar_indels(curr_read_id, read, region.contig, self.config, use_clips, read_nm=nm)
                        nm_adj = (nm_raw - (ins_sum + del_sum))
                        nm_adj_ratio = nm_adj / float(read.query_alignment_length + 1)
                        nm = nm_adj_ratio
                        nm_sum += nm
                        nm_count += 1

                if phase:
                    curr_read_id = (self.read_id, str(read.get_tag("HP")) if read.has_tag("HP") else "NULL", str(read.get_tag("PS")) if read.has_tag("PS") else "NULL")

            if trace_read is not False:
                if trace_read == read.query_name:
                    print(f"[DEV_TRACE_READ] [0b/4] [LeadProvider.iter_region] [{region}] [{read.query_name}] passed pre-filtering (whole-read), begin to extract leads")

            # Extract small indels
            for lead in read_iterindels(curr_read_id, read, region.contig, self.config, use_clips, read_nm=nm):
                if trace_read is not False:
                    if trace_read == read.query_name:
                        print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_iterindels] [{region}] [{read.query_name}] new lead: {lead}")
                yield lead

            # Extract read splits
            if has_sa:
                if read.is_supplementary:
                    if trace_read is not False:
                        if trace_read == read.query_name:
                            print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits_bnd] [{region}] [{read.query_name}] is entering read_itersplits_bnd")
                    for lead in read_itersplits_bnd(curr_read_id, read, region.contig, self.config, read_nm=nm):
                        if trace_read is not False:
                            if trace_read == read.query_name:
                                print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits_bnd] [{region}] [{read.query_name}] new lead: {lead}")
                        yield lead
                else:
                    if trace_read is not False:
                        if trace_read == read.query_name:
                            print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits] [{region}] [{read.query_name}] is entering read_itersplits")
                    for lead in read_itersplits(curr_read_id, read, region.contig, self.config, read_nm=nm):
                        if trace_read is not False:
                            if trace_read == read.query_name:
                                print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits] [{region}] [{read.query_name}] new lead: {lead}")
                        yield lead
                read_mq20 += 1 if read.mapping_quality >= 20 else 0

            # Record in coverage table
            read_end = read.reference_start + read.reference_length
            assert (read_end == read.reference_end)
            # assert(read_end>=read.reference_start)
            if read.is_reverse:
                target_tab = self.covrtab_rev
            else:
                target_tab = self.covrtab_fwd
            covr_start_bin = (int(read.reference_start / coverage_binsize) + coverage_shift_bins * (alen >= coverage_shift_min_aln_len)) * coverage_binsize
            covr_end_bin = (int(read_end / coverage_binsize) - coverage_shift_bins * (alen >= coverage_shift_min_aln_len)) * coverage_binsize

            if covr_end_bin > covr_start_bin:
                self.covrtab_min_bin = min(self.covrtab_min_bin, covr_start_bin)
                target_tab[covr_start_bin] = target_tab[covr_start_bin] + 1 if covr_start_bin in target_tab else 1

                if read_end <= self.end:
                    target_tab[covr_end_bin] = target_tab[covr_end_bin] - 1 if covr_end_bin in target_tab else -1

        self.config.average_regional_nm = nm_sum / float(max(1, nm_count))
        self.config.qc_nm_threshold = self.config.average_regional_nm
        # print(f"Contig {contig} avg. regional NM={self.config.average_regional_nm}, threshold={self.config.qc_nm_threshold}")

    def dev_leadtab_filename(self, contig, start, end):
        scriptloc = os.path.dirname(os.path.realpath(sys.argv[0]))
        if self.config.dev_cache_dir is None:
            cache_dir = f"{scriptloc}/cache"
        else:
            cache_dir = self.config.dev_cache_dir
        return f"{cache_dir}/{os.path.basename(self.config.input)}_{contig}_{start}_{end}.pickle"
