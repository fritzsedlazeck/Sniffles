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
import time
from dataclasses import dataclass
from typing import Optional, Iterator

import numpy as np
import pysam
import array

# for: --dev-cache
import os
import sys
# end: for: --dev-cache

from sniffles import util
from sniffles import sv
from sniffles.region import Region
from sniffles.sv import SVCallBNDInfo

log = logging.getLogger(__name__)
UNSIGEND_SHORT = 'H'


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
    bnd_info: Optional[SVCallBNDInfo] = None
    hap: str = "0"
    phase_set: str = None
    is_sa: bool = False
    read_len: int = 0

    @classmethod
    def for_bnd(cls, read_id: int, read: pysam.AlignedSegment) -> Optional['Lead']:
        """
        Generate a Lead based on the supplementary alignment (SA tag) of the given read.
        """
        # determine side of supplementary alignment
        # only works properly for a single SA
        left = 0
        right = 0
        if read.cigartuples[0][0] in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):
            left = read.cigartuples[0][1]
        if read.cigartuples[-1][0] in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):
            right = read.cigartuples[-1][1]

        if left > right:
            ref_start = read.reference_start + 1
            is_first = False
        else:
            ref_start = read.reference_end
            is_first = True

        read_strand = "-" if read.is_reverse else "+"

        for sa in read.get_tag("SA").split(";"):
            if len(sa) > 0:
                refname, pos, strand, cigar, mapq, nm = sa.split(",")
                pos_org = pos
                pos = int(pos) - 1  # position in SA tag is 1-based, convert to 0-based
                # determine left/right on SA based on cigar string
                try:
                    left, right, refspan, readspan = CIGAR_analyze(cigar)
                except Exception as e:
                    util.error(f"Malformed CIGAR '{cigar}' with pos {pos_org} of read '{read.query_name}' ({e}). Skipping.")
                    return None
                is_reverse = right > left
                if is_reverse:
                    if is_first:
                        mate_ref_start = pos + refspan
                    else:
                        mate_ref_start = pos + refspan
                else:
                    if is_first:
                        mate_ref_start = pos + 1
                    else:
                        if strand == read_strand:
                            mate_ref_start = pos + refspan
                        else:
                            mate_ref_start = pos + 2

                return cls(
                    read_id=read_id,
                    read_qname=read.query_name,
                    contig=read.reference_name,
                    ref_start=ref_start,
                    ref_end=ref_start,
                    qry_start=read.query_alignment_start,
                    qry_end=read.query_alignment_end,
                    strand="-" if read.is_reverse else "+",
                    mapq=read.mapping_quality,
                    nm=int(nm) if read.has_tag("NM") else None,
                    source="BND_SA",
                    svtype="BND",
                    svlen=0,
                    seq=None,
                    bnd_info=SVCallBNDInfo(
                        mate_contig=refname,
                        mate_ref_start=mate_ref_start,
                        is_first=is_first,
                        is_reverse=is_reverse,
                    )
                )

        return None

    def classify(self):
        """
        Determine a potentially more descriptive classification of this BND lead (e.g. INV)
        """
        if self.svtype != 'BND' or self.bnd_info is None:
            return




def CIGAR_analyze(cigar) -> tuple[int, int, int, int]:
    """
    Analyze CIGAR string to determine clips on start/end
    :returns 4-tuple of (clipping at the start, clipping at the end, refspan, readspan)
    """
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


def get_cigar_indels(read: pysam.AlignedSegment, minoplen: int = 10) -> tuple[int, int, int, int]:
    """
    Counts the number of inserted and deleted bases in a read.
    """
    cins = pysam.CINS
    cdel = pysam.CDEL

    ins_sum = 0
    del_sum = 0
    large_ins_sum = 0
    large_del_sum = 0

    pos_read = 0
    pos_ref = read.reference_start
    for op, oplength in read.cigartuples:
        add_read, add_ref, event, _ = OPLIST[op]
        if event:
            if op == cins:
                ins_sum += oplength
                large_ins_sum += oplength if oplength > minoplen else 0
            elif op == cdel:
                del_sum += oplength
                large_del_sum += oplength if oplength > minoplen else 0
        pos_read += add_read * oplength
        pos_ref += add_ref * oplength

    return ins_sum, del_sum, large_ins_sum, large_del_sum


def read_itersplits(read_id, read, contig, config, read_nm, read_hap, read_ps) -> Iterator[Lead]:
    """
    Iterate over supplementary alignments of a primary read
    """
    # SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads = []
    supps = [part.split(",") for part in read.get_tag("SA").split(";") if len(part) > 0]
    trace_read = config.dev_trace_read and read.query_name in config.dev_trace_read

    if len(supps) > config.max_splits_base + config.max_splits_kb * (read.query_length / 1000.0):
        return

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] passed max_splits check")

    if read.is_reverse:
        qry_start = read.query_alignment_end
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
                     "?",
                     hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary)
    all_leads.append(curr_lead)

    # unterschied zu _bnd: alle supplementary alignments werden verarbeitet
    for refname, pos, strand, cigar, mapq, nm in supps:
        mapq = int(mapq)
        nm = int(nm)

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
                              "?",
                              hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary))

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] all_leads: {all_leads}")

    all_leads = sv.classify_splits(read, all_leads, config, contig)

    if trace_read:
        print(f"[DEV_TRACE_READ] [0c/4] [LeadProvider.read_itersplits] [{read.query_name}] "
              f"classify_splits(all_leads): {all_leads}")

    """
    if config.dev_trace_read != False:
        print(read.query_name)
        if read.query_name in config.dev_trace_read:
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
                           seq=None,
                           hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary)
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
                           seq=lead.seq if svtype == "INS" else None,
                           hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary)


class LeadProvider:
    coverage: np.ndarray

    def __init__(self, config, read_id_offset, contig: str):
        self.config = config

        self.leadtab = {}
        self.leadcounts = {}
        self.leadhapcount = {}

        for svtype in sv.ALL_TYPES:
            self.leadtab[svtype] = {}
            self.leadcounts[svtype] = 0
            self.leadhapcount[svtype] = {}
        self.leadhapcount["REF"] = {}  # we need to count the haps in ref too,

        self.covrtab_fwd = {}
        self.covrtab_rev = {}
        self.covrtab_min_bin = None

        self.read_id = read_id_offset
        self.read_count = 0

        self.contig = contig
        self.start = None
        self.end = None

    def record_hap_ref(self, hp_index, pos_leadtab, end_leadtab, step):
        leadtab_hapc = self.leadhapcount["REF"]
        for this_pos in range(pos_leadtab, end_leadtab, step):
            if this_pos in leadtab_hapc:
                leadtab_hapc[this_pos][hp_index] += 1
            else:
                leadtab_hapc[this_pos] = array.array(UNSIGEND_SHORT, 3*[0])
                leadtab_hapc[this_pos][hp_index] = 1

    def record_lead(self, ld, pos_leadtab):
        leadtab_svtype = self.leadtab[ld.svtype]
        leadtab_hapc = self.leadhapcount[ld.svtype]
        hp_index = int(ld.hap)
        if pos_leadtab in leadtab_svtype:
            leadtab_svtype[pos_leadtab].append(ld)
            lead_count = len(leadtab_svtype[pos_leadtab])
            if lead_count > self.config.consensus_max_reads_bin:
                ld.seq = None
            leadtab_hapc[pos_leadtab][hp_index] += 1
        else:
            leadtab_svtype[pos_leadtab] = [ld]
            lead_count = 1
            leadtab_hapc[pos_leadtab] = array.array(UNSIGEND_SHORT, 3*[0])
            leadtab_hapc[pos_leadtab][hp_index] = 1
        self.leadcounts[ld.svtype] += 1

    def _mask_N_coverage(self, regions: list[Region] = None) -> None:
        """
        Mask coverage with 0 for regions with N in the reference
        """
        if self.config.reference:
            start = time.monotonic()
            try:
                fasta = pysam.FastaFile(self.config.reference)
            except Exception as e:
                logging.warning(f'Unable to mask N regions in coverage vector, reference could not be opened: {e}', exc_info=True)
            else:
                try:
                    if regions is None:
                        mask = np.frombuffer(fasta.fetch(self.contig).encode('ascii'), dtype=np.uint8)
                    else:
                        mask = np.zeros(len(self.coverage), dtype=np.uint8)
                        for region in regions:
                            region_mask = np.frombuffer(fasta.fetch(region.contig, region.start, region.end).encode('ascii'), dtype=np.uint8)
                            mask[region.start:region.end] = region_mask
                    self.coverage[mask == 78] = 0
                except Exception as e:
                    logging.warning(f'Unable to mask N regions in coverage vector, reference could not be fetched: {e}', exc_info=True)
                else:
                    log.debug(f'Masking coverage with {self.config.reference} reference (elapsed: {time.monotonic() - start:.2f}s)')

    def build_leadtab(self, regions: list[Region], bam):
        assert (self.start is None)
        assert (self.end is None)

        externals = []
        ld_binsize = self.config.cluster_binsize
        self.coverage = np.zeros(bam.get_reference_length(self.contig), dtype=np.uint16)

        for region in regions:
            if self.contig != region.contig:
                raise Exception(f"Region contig '{region.contig}' does not match LeadProvider contig '{self.contig}'")

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

        self._mask_N_coverage(regions)

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
        ld_binsize = self.config.cluster_binsize

        read_mq20 = 0
        read: pysam.AlignedSegment
        for read in bam.fetch(region.contig, region.start, region.end, until_eof=False):
            if trace_read:
                if read.query_name in trace_read:
                    print(f"[DEV_TRACE_READ] [0b/4] [LeadProvider.iter_region] [{region}] [{read.query_name}] "
                          f"has been fetched and is entering pre-filtering")

            alen = read.query_alignment_length
            if read.mapping_quality < mapq_min or read.is_secondary or alen < alen_min:
                continue

            if exclude_flags is not None:
                if read.flag & exclude_flags:
                    continue

            if read.reference_start < region.start or read.reference_start >= region.end:
                continue

            hp = read.get_tag("HP") if read.has_tag("HP") else 0
            ps = read.get_tag("PS") if read.has_tag("PS") else "NULL"

            self.read_id += 1
            self.read_count += 1
            self.coverage[read.reference_start:read.reference_end] += 1

            has_sa = read.has_tag("SA")
            use_clips = self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm = -1
            curr_read_id = self.read_id
            if advanced_tags:
                if qc_nm or True:  # always save
                    if read.has_tag("NM"):
                        nm_raw = read.get_tag("NM")
                        ins_sum, del_sum, *large_op = get_cigar_indels(read)
                        # NOTE: previously all INS/DEL events sizes were summed up and removed from NM (nm_adj)
                        #       now only those above a given threshold (default = 10 bases)
                        nm = (nm_raw - sum(large_op)) / float(read.query_alignment_length + 1)
                        nm_sum += nm
                        nm_count += 1
                        # nm_adj = (nm_raw - (ins_sum + del_sum))
                        # nm_adj_ratio = nm_adj / float(read.query_alignment_length + 1)
                        # nm = nm_adj_ratio

            if trace_read:
                if read.query_name in trace_read:
                    print(f"[DEV_TRACE_READ] [0b/4] [LeadProvider.iter_region] [{region}] [{read.query_name}] "
                          f"passed pre-filtering (whole-read), begin to extract leads")

            # Extract small indels and record coverage for this read
            for lead in self.read_iterindels(curr_read_id, read, region.contig, use_clips, read_nm=nm, read_hap=hp,
                                             read_ps=ps):
                if trace_read:
                    if read.query_name in trace_read:
                        print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_iterindels] [{region}] [{read.query_name}] "
                              f"new lead: {lead}")
                yield lead

            # Extract read splits
            if has_sa:
                # New version for BNDs:
                if lead := Lead.for_bnd(curr_read_id, read):
                    yield lead

                # Old version for everything else:
                # TODO: may yield duplicates of events called as generic BNDs on top
                if not read.is_supplementary:
                    if trace_read:
                        if trace_read == read.query_name:
                            print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits] [{region}] [{read.query_name}] is entering read_itersplits")
                    for lead in read_itersplits(
                            curr_read_id, read, region.contig, self.config, read_nm=nm, read_hap=hp, read_ps=ps
                    ):
                        if trace_read:
                            if read.query_name in trace_read:
                                print(f"[DEV_TRACE_READ] [1/4] [leadprov.read_itersplits] [{region}] [{read.query_name}] new lead: {lead}")
                        yield lead
                read_mq20 += 1 if read.mapping_quality >= self.config.mapq else 0

            # we record counts for each lead-tab for all haplotype count appearances
            ld_ref_start, ld_ref_end = read.reference_start, read.reference_end
            pos_leadtab = int(ld_ref_start / ld_binsize) * ld_binsize
            end_leadtab = int(ld_ref_end / ld_binsize) * ld_binsize
            if region.start <= ld_ref_start < region.end:
                self.record_hap_ref(hp, pos_leadtab, end_leadtab, ld_binsize)

        log.debug(f'Processed {self.read_count} reads in region {region.contig}:{region.start}-{region.end}')

        self.config.average_regional_nm = nm_sum / float(max(1, nm_count))
        self.config.qc_nm_threshold = self.config.average_regional_nm

        # print(f"Contig {contig} avg. regional NM={self.config.average_regional_nm},
        # threshold={self.config.qc_nm_threshold}")

    def read_iterindels(self, read_id: int, read, contig, use_clips, read_nm, read_hap, read_ps) -> Iterator[Lead]:
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
        read_len = read.query_alignment_length  # or read.infer_read_length()

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
                               seq=read.query_sequence[pos_read:pos_read + oplength] if oplength <= seq_cache_maxlen else None,
                               hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary, read_len=read_len)
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
                               -oplength,
                               hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary, read_len=read_len)
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
                               seq=None,
                               hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary, read_len=read_len)
                elif op in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):
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
                               "SINGLE_LEFT" if pos_ref == read.reference_start else "SINGLE_RIGHT",
                               0,
                               seq=None,
                               hap=str(read_hap), phase_set=str(read_ps), is_sa=read.is_supplementary, read_len=read_len)

            pos_read += add_read * oplength
            pos_ref += add_ref * oplength

    def dev_leadtab_filename(self, contig, start, end):
        scriptloc = os.path.dirname(os.path.realpath(sys.argv[0]))
        if self.config.dev_cache_dir is None:
            cache_dir = f"{scriptloc}/cache"
        else:
            cache_dir = self.config.dev_cache_dir
        return f"{cache_dir}/{os.path.basename(self.config.input)}_{contig}_{start}_{end}.pickle"
