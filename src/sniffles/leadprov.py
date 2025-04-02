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
import logging
from dataclasses import dataclass
import itertools
from typing import Optional, Iterator

import numpy as np
import pysam

# for: --dev-cache
import os
import sys
# end: for: --dev-cache

from sniffles import util
from sniffles import sv
from sniffles.region import Region


log = logging.getLogger(__name__)


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
                    if clip_start is None and readspan + refspan > 0:
                        clip_start = clip
                    clip += oplen
                else:
                    raise f"Unknown CIGAR operation: '{c}'"
            buf = ""
    if clip_start is None:
        clip_start = clip
    return clip_start, clip - clip_start, refspan, readspan


# Tuples are (add_read, add_ref, event, Δ coverage)
OPTAB = {
    pysam.CMATCH:     (1, 1, 0, 1),
    pysam.CEQUAL:     (1, 1, 0, 1),
    pysam.CDIFF:      (1, 1, 0, 1),
    pysam.CINS:       (1, 0, 1, 0),
    pysam.CDEL:       (0, 1, 1, 0),
    pysam.CREF_SKIP:  (0, 1, 0, 0),
    pysam.CSOFT_CLIP: (1, 0, 1, 0),
    pysam.CHARD_CLIP: (0, 0, 0, 0),
    pysam.CPAD:       (0, 0, 0, 0),
}


OPLIST = [(0, 0, 0, 0) for i in range(max(int(k) for k in OPTAB.keys()) + 1)]
for k, v in OPTAB.items():
    OPLIST[int(k)] = v


def get_cigar_indels(read) -> tuple[int, int]:
    """
    Counts the number of inserted and deleted bases in a read.
    """
    cins = pysam.CINS
    cdel = pysam.CDEL

    ins_sum = 0
    del_sum = 0

    pos_read = 0
    pos_ref = read.reference_start
    for op, oplength in read.cigartuples:
        add_read, add_ref, event, _ = OPLIST[op]
        if event:
            if op == cins:
                ins_sum += oplength
            elif op == cdel:
                del_sum += oplength
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength

    return ins_sum, del_sum


def read_itersplits_bnd(read_id, read, contig, config, read_nm):
    """
    Iterate over supplementary alignments of a read and yield leads for BND events.
    """
    assert read.is_supplementary
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
        # Primary alignment is on this chromosome, no need to parse the supplementary, but only if mapping quality of primary is sufficient
        if int(prim_mapq) >= config.mapq:
            # Skip if primary alignment is of sufficient quality
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
                           ref_end=svstart + svlen if svlen is not None and svtype != "INS" else svstart,
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
    coverage: np.ndarray

    def __init__(self, config, read_id_offset):
        self.config = config

        self.leadtab = {}
        self.leadcounts = {}

        for svtype in sv.ALL_TYPES:
            self.leadtab[svtype] = {}
            self.leadcounts[svtype] = 0

        self.covrtab_fwd = {}
        self.covrtab_rev = {}
        self.covrtab_min_bin = None

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
        qc_nm = self.config.qc_nm_measure
        phase = self.config.phase
        advanced_tags = qc_nm or phase
        mapq_min = self.config.mapq
        alen_min = self.config.min_alignment_length
        exclude_flags = self.config.exclude_flags
        nm_sum = 0
        nm_count = 0
        trace_read = self.config.dev_trace_read
        self.coverage = np.zeros(bam.get_reference_length(region.contig), dtype=np.uint16)

        read_mq20 = 0
        for read in bam.fetch(region.contig, region.start, region.end, until_eof=False):
            if trace_read is not False:
                if trace_read == read.query_name:
                    print(f"[DEV_TRACE_READ] [0b/4] [LeadProvider.iter_region] [{region}] [{read.query_name}] has been fetched and is entering pre-filtering")

            alen = read.query_alignment_length
            if read.mapping_quality < mapq_min or read.is_secondary or alen < alen_min:
                continue

            if exclude_flags is not None:
                if read.flag & exclude_flags:
                    continue

            if read.reference_start < region.start or read.reference_start >= region.end:
                continue

            self.read_id += 1
            self.read_count += 1
            self.coverage[read.reference_start:read.reference_end] += 1

            has_sa = read.has_tag("SA")
            use_clips = self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm = -1
            curr_read_id = self.read_id
            if advanced_tags:
                if qc_nm:
                    if read.has_tag("NM"):
                        nm_raw = read.get_tag("NM")
                        ins_sum, del_sum = get_cigar_indels(read)
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

            # Extract small indels and record coverage for this read
            for lead in self.read_iterindels(curr_read_id, read, region.contig, use_clips, read_nm=nm):
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

        log.debug(f'Processed {self.read_count} reads in region {region.contig}:{region.start}-{region.end}')

        self.config.average_regional_nm = nm_sum / float(max(1, nm_count))
        self.config.qc_nm_threshold = self.config.average_regional_nm
        # print(f"Contig {contig} avg. regional NM={self.config.average_regional_nm}, threshold={self.config.qc_nm_threshold}")

    def read_iterindels(self, read_id: int, read, contig, use_clips, read_nm) -> Iterator[Lead]:
        """
        Extracts Leads from insertions/deletions.
        """
        config = self.config
        coverage = self.coverage
        minsvlen = config.minsvlen_screen
        longinslen = config.long_ins_length / 2.0
        seq_cache_maxlen = config.dev_seq_cache_maxlen
        qname = read.query_name
        mapq = read.mapping_quality
        strand = "-" if read.is_reverse else "+"
        CINS = pysam.CINS  # noqa lowercase
        CDEL = pysam.CDEL  # noqa lowercase
        CSOFT_CLIP = pysam.CSOFT_CLIP  # noqa lowercase

        pos_read = 0
        pos_ref = read.reference_start

        for op, oplength in read.cigartuples:
            add_read, add_ref, event, dcov = OPLIST[op]
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

    def dev_leadtab_filename(self, contig, start, end):
        scriptloc = os.path.dirname(os.path.realpath(sys.argv[0]))
        if self.config.dev_cache_dir is None:
            cache_dir = f"{scriptloc}/cache"
        else:
            cache_dir = self.config.dev_cache_dir
        return f"{cache_dir}/{os.path.basename(self.config.input)}_{contig}_{start}_{end}.pickle"
