#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 28.05.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#
from dataclasses import dataclass
import re
import itertools

import pysam

#for: --dev-cache
import os
import sys
import pickle
#end: for: --dev-cache

from sniffles import util
from sniffles import sv

@dataclass
class Lead:
    read_id: int=None
    read_qname: str=None
    contig: str=None
    ref_start: int=None
    ref_end: int=None
    qry_start: int=None
    qry_end: int=None
    strand: str=None
    mapq: int=None
    nm: float=None
    source: str=None
    svtype: str=None
    svlen: int=None
    seq: str=None
    svtypes_starts_lens: list=None

def CIGAR_analyze(cigar):
    buf=""
    readspan=0
    refspan=0
    clip_start=None
    clip=0
    for c in cigar:
        if c.isnumeric():
            buf+=c
        else:
            oplen=int(buf)
            h=False
            if c in "MIX=":
                readspan+=oplen
                h=True
            if c in "MDX=N":
                refspan+=oplen
                h=True
            if not h:
                if c in "SH":
                    if clip_start==None and readspan+refspan>0:
                        clip_start=clip
                    clip+=oplen
                else:
                    raise f"Unknown CIGAR operation: '{c}'"
            buf=""
    if clip_start==None:
        clip_start=clip
    return clip_start, clip-clip_start, refspan, readspan

def CIGAR_analyze_regex(cigar):
    #TODO: Obsolete
    opsums={"M":0,"I":0,"D":0,"=":0,"X":0,"N":0,"H":0,"S":0}
    iter=re.split(r"(\d+)", cigar)
    for i in range(1,len(iter)-1,2):
        op=iter[i+1]
        if op!="H" and op!="S":
            readstart_fwd=opsums["H"]+opsums["S"]
        opsums[iter[i+1]]+=int(iter[i])
    readstart_rev=opsums["H"]+opsums["S"]-readstart_fwd
    refspan=opsums["M"]+opsums["D"]+opsums["="]+opsums["X"]+opsums["N"]
    readspan=opsums["M"]+opsums["I"]+opsums["="]+opsums["X"]
    return readstart_fwd,readstart_rev,refspan,readspan

def CIGAR_tolist(cigar):
    #TODO: Obsolete (see CIGAR_tolist_analyze)
    """
    CIGAR string : str -> List of CIGAR operation & length tuples : [(op1:char, op1_length:int),...]
    """
    buf=""
    ops=[]
    for c in cigar:
        if c.isnumeric():
            buf+=c
        else:
            ops.append((c,int(buf)))
            buf=""
    return ops

def CIGAR_listrefspan(ops):
    #TODO: Obsolete (see CIGAR_analyze)
    #TODO(Potential): Detect&utilize minimap2 condensed supplementary alignment cigar for speed
    return sum(oplen for op,oplen in ops if op=="M" or op=="D" or op=="=" or op=="X" or op=="N")

def CIGAR_listreadspan(ops):
    #TODO: Obsolete (see CIGAR_analyze)
    return sum(oplen for op,oplen in ops if op=="M" or op=="I" or op=="=" or op=="X")

def CIGAR_listreadstart_fwd(ops):
    #TODO: Obsolete (see CIGAR_analyze)
    """
    Position in query read where CIGAR alignment starts (i.e. taking into account start clipping)
    """
    op,oplen=ops[0]
    op2,op2len=ops[1]
    if op=="H" or op=="S":
        assert(op2!="H" and op2!="S")
        return oplen
    else:
        return 0

def CIGAR_listreadstart_rev(ops):
    #TODO: Obsolete (see CIGAR_analyze)
    """
    Position in query read where CIGAR alignment starts (i.e. taking into account start clipping)
    """
    op,oplen=ops[-1]
    op2,op2len=ops[-2]
    if op=="H" or op=="S":
        assert(op2!="H" and op2!="S")
        return oplen
    else:
        return 0

OPTAB={pysam.CMATCH:     (1,1,0),
       pysam.CEQUAL:     (1,1,0),
       pysam.CDIFF:      (1,1,0),
       pysam.CINS:       (1,0,1),
       pysam.CDEL:       (0,1,1),
       pysam.CREF_SKIP:  (0,1,0),
       pysam.CSOFT_CLIP: (1,0,1),
       pysam.CHARD_CLIP: (0,0,0),
       pysam.CPAD:       (0,0,0)}
#      pysam.CBACK:      (0,0,0)}

OPLIST=[(0,0,0) for i in range(max(int(k) for k in OPTAB.keys())+1)]
for k,v in OPTAB.items():
    OPLIST[int(k)]=v

def read_iterindels(read_id,read,contig,config,use_clips,read_nm):
    minsvlen=config.minsvlen_screen
    longinslen=config.long_ins_length/2.0
    seq_cache_maxlen=config.dev_seq_cache_maxlen
    qname=read.query_name
    mapq=read.mapping_quality
    strand="-" if read.is_reverse else "+"
    CINS=pysam.CINS
    CDEL=pysam.CDEL
    CSOFT_CLIP=pysam.CSOFT_CLIP

    pos_read=0
    pos_ref=read.reference_start
    for op,oplength in read.cigartuples:
        add_read,add_ref,event=OPLIST[op]
        if event and oplength >= minsvlen:
            if op==CINS:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read+oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           oplength,
                           seq=read.query_sequence[pos_read:pos_read+oplength] if oplength <= seq_cache_maxlen else None)
            elif op==CDEL:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref+oplength,
                           pos_ref,
                           pos_read,
                           pos_read,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "DEL",
                           -oplength)
            elif use_clips and op==CSOFT_CLIP and oplength >= longinslen:
                yield Lead(read_id,
                           qname,
                           contig,
                           pos_ref,
                           pos_ref,
                           pos_read,
                           pos_read+oplength,
                           strand,
                           mapq,
                           read_nm,
                           "INLINE",
                           "INS",
                           None,
                           seq=None)
        pos_read+=add_read*oplength
        pos_ref+=add_ref*oplength

def read_iterindels_unoptimized(read_id,read,contig,config,use_clips):
    minsvlen=config.minsvlen_screen
    seq_cache_maxlen=config.dev_seq_cache_maxlen
    qname=read.query_name
    mapq=read.mapping_quality
    strand="-" if read.is_reverse else "+"

    #TODO: Parse CG tag (ultra long alignments), if present
    pos_read=0
    pos_ref=read.reference_start
    for op,oplength in read.cigartuples:
        if op==pysam.CMATCH or op==pysam.CEQUAL or op==pysam.CDIFF:
            pos_read+=oplength
            pos_ref+=oplength
        elif op==pysam.CINS:
            if oplength>=minsvlen:
                #print(pos_read,pos_read+oplength)
                #print(pos_read,pos_read+oplength,read.query_sequence[pos_read:pos_read+oplength])
                if oplength <= seq_cache_maxlen:
                    seq=read.query_sequence[pos_read:pos_read+oplength]
                else:
                    seq=None
                yield Lead(read_id,qname,contig,pos_ref,pos_ref+0,pos_read,pos_read+oplength,strand,mapq,-1,"INLINE","INS",oplength,seq=seq)
            pos_read+=oplength
        elif op==pysam.CDEL:
            pos_ref+=oplength
            if oplength>=minsvlen:
                yield Lead(read_id,qname,contig,pos_ref,pos_ref+oplength,pos_read,pos_read+0,strand,mapq,-1,"INLINE","DEL",-oplength)
        elif op==pysam.CREF_SKIP:
            pos_ref+=oplength
        elif op==pysam.CSOFT_CLIP:
            if use_clips and oplength >= config.long_ins_length:
                yield Lead(read_id,qname,contig,pos_ref,pos_ref+0,pos_read,pos_read+oplength,strand,mapq,-1,"INLINE","INS",None)
            pos_read+=oplength
        elif op==pysam.CHARD_CLIP:
            #pos_ref+=oplength
            if use_clips and oplength >= config.long_ins_length:
                yield Lead(read_id,qname,contig,pos_ref,pos_ref+0,pos_read,pos_read+oplength,strand,mapq,-1,"INLINE","INS",None)
        elif op==pysam.CPAD:
            pass
        else:
            print(f"Unknown OPType {op}")
            return

def read_itersplits_bnd(read_id,read,contig,config,read_nm):
    assert(read.is_supplementary)
    #SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads=[]
    supps=[part.split(",") for part in read.get_tag("SA").split(";") if len(part)>0]

    if len(supps) > config.max_splits_base + config.max_splits_kb*(read.query_length/1000.0):
        return

    if read.is_reverse:
        qry_start=read.query_length-read.query_alignment_end
    else:
        qry_start=read.query_alignment_start

    curr_lead=Lead(read_id,
                   read.query_name,
                   contig,
                   read.reference_start,
                   read.reference_start+read.reference_length,
                   qry_start,
                   qry_start+read.query_alignment_length,
                   "-" if read.is_reverse else "+",
                   read.mapping_quality,
                   read_nm,
                   "SPLIT_SUP",
                   "?")
    all_leads.append(curr_lead)

    prim_refname,prim_pos,prim_strand,prim_cigar,prim_mapq,prim_nm=supps[0]
    if prim_refname == contig:
        #Primary alignment is on this chromosome, no need to parse the supplementary
        return

    minpos_curr_chr=min(itertools.chain([read.reference_start],(int(pos) for refname,pos,strand,cigar,mapq,nm in supps if refname==contig)))
    if minpos_curr_chr < read.reference_start:
        #Only process splits once per chr (there may be multiple supplementary alignments on the same chr)
        return

    for refname,pos,strand,cigar,mapq,nm in supps:
        mapq=int(mapq)
        nm=int(nm)
        #if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev=(strand=="-")

        try:
            readstart_fwd,readstart_rev,refspan,readspan=CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero=int(pos)-1
        split_qry_start=readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start+readspan,
                              strand,
                              mapq,
                              nm/float(readspan+1),
                              "SPLIT_SUP",
                              "?"))

    sv.classify_splits(read,all_leads,config,contig)

    for lead in all_leads:
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            if svtype=="BND":
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
                bnd.bnd_info=arg
                #print(lead.contig,svstart,bnd.bnd_info)
                yield bnd

def read_itersplits(read_id,read,contig,config,read_nm):
    #SA:refname,pos,strand,CIGAR,MAPQ,NM
    all_leads=[]
    supps=[part.split(",") for part in read.get_tag("SA").split(";") if len(part)>0]

    if len(supps) > config.max_splits_base + config.max_splits_kb*(read.query_length/1000.0):
        return

    #QC on: 18Aug21, HG002.ont.chr22; O.K.
    #cigarl=CIGAR_tolist(read.cigarstring)
    #if read.is_reverse:
    #    cigarl.reverse()

    #if read.is_reverse:
    #    assert(read.query_length-read.query_alignment_end == CIGAR_listreadstart(cigarl))
    #else:
    #    assert(read.query_alignment_start == CIGAR_listreadstart(cigarl))

    #assert(CIGAR_listrefspan(cigarl)==read.reference_length)
    #assert(CIGAR_listreadspan(cigarl)==read.query_alignment_length)
    #End QC

    if read.is_reverse:
        qry_start=read.query_length-read.query_alignment_end
    else:
        qry_start=read.query_alignment_start

    curr_lead=Lead(read_id,
                   read.query_name,
                   contig,
                   read.reference_start,
                   read.reference_start+read.reference_length,
                   qry_start,
                   qry_start+read.query_alignment_length,
                   "-" if read.is_reverse else "+",
                   read.mapping_quality,
                   read_nm,
                   "SPLIT_PRIM",
                   "?")
    all_leads.append(curr_lead)

    #QC on: 18Aug21; O.K.
    #assert(read.reference_length == CIGAR_listrefspan(CIGAR_tolist(read.cigarstring)))
    #assert(read.query_alignment_start == CIGAR_listreadstart(CIGAR_tolist(read.cigarstring)))
    #assert(read.query_alignment_length == CIGAR_listreadspan(CIGAR_tolist(read.cigarstring)))
    #End QC

    for refname,pos,strand,cigar,mapq,nm in supps:
        mapq=int(mapq)
        nm=int(nm)
        #if not config.dev_keep_lowqual_splits and mapq < config.mapq:
        #    continue

        is_rev=(strand=="-")

        try:
            readstart_fwd,readstart_rev,refspan,readspan=CIGAR_analyze(cigar)
        except Exception as e:
            util.error(f"Malformed CIGAR '{cigar}' with pos {pos} of read '{read.query_name}' ({e}). Skipping.")
            return

        pos_zero=int(pos)-1
        split_qry_start=readstart_rev if is_rev else readstart_fwd

        all_leads.append(Lead(read_id,
                              read.query_name,
                              refname,
                              pos_zero,
                              pos_zero + refspan,
                              split_qry_start,
                              split_qry_start+readspan,
                              strand,
                              mapq,
                              nm/float(readspan+1),
                              "SPLIT_SUP",
                              "?"))

        #QC on: 08Sep21; O.K.
        #cigarl=CIGAR_tolist(cigar)
        #assert(CIGAR_listrefspan(cigarl)==refspan)
        #assert(CIGAR_listreadspan(cigarl)==readspan)
        #assert(CIGAR_listreadstart_fwd(cigarl)==readstart_fwd)
        #assert(CIGAR_listreadstart_rev(cigarl)==readstart_rev)
        #End QC

    sv.classify_splits(read,all_leads,config,contig)

    for lead_i, lead in enumerate(all_leads):
        for svtype, svstart, arg in lead.svtypes_starts_lens:
            min_mapq=min(lead.mapq,all_leads[max(0,lead_i-1)].mapq)
            if not config.dev_keep_lowqual_splits and min_mapq < config.mapq:
                continue

            if svtype=="BND":
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
                bnd.bnd_info=arg
                yield bnd

            elif svtype!="NOSV":
                svlen=arg

                yield Lead(read_id=lead.read_id,
                           read_qname=lead.read_qname,
                           contig=lead.contig,
                           ref_start=svstart,
                           ref_end=svstart+svlen if svlen!=None and svtype!="INS" else svstart,
                           qry_start=lead.qry_start,
                           qry_end=lead.qry_end,
                           strand=lead.strand,
                           mapq=lead.mapq,
                           nm=lead.nm,
                           source=lead.source,
                           svtype=svtype,
                           svlen=svlen,
                           seq=lead.seq if svtype=="INS" else None)

class LeadProvider:
    def __init__(self,config,read_id_offset):
        self.config=config

        self.leadtab={}
        self.leadcounts={}

        for svtype in sv.TYPES:
            self.leadtab[svtype]={}
            self.leadcounts[svtype]=0

        self.covrtab_fwd={}
        self.covrtab_rev={}
        self.covrtab_min_bin=None
        #self.covrtab_read_start={}
        #self.covrtab_read_end={}

        self.read_id=read_id_offset
        self.read_count=0

        self.contig=None
        self.start=None
        self.end=None

    def record_lead(self,ld,pos_leadtab):
        leadtab_svtype=self.leadtab[ld.svtype]
        if pos_leadtab in leadtab_svtype:
            leadtab_svtype[pos_leadtab].append(ld)
            lead_count=len(leadtab_svtype[pos_leadtab])
            if lead_count > self.config.consensus_max_reads_bin:
                ld.seq=None
        else:
            leadtab_svtype[pos_leadtab]=[ld]
            lead_count=1
        self.leadcounts[ld.svtype]+=1

    def build_leadtab(self,contig,start,end,bam):
        if self.config.dev_cache:
            loaded_externals=self.dev_load_leadtab(contig,start,end)
            if loaded_externals!=False:
                return loaded_externals

        assert(self.contig==None)
        assert(self.start==None)
        assert(self.end==None)
        self.contig=contig
        self.start=start
        self.end=end
        self.covrtab_min_bin=int(self.start/self.config.coverage_binsize)*self.config.coverage_binsize

        externals=[]
        ld_binsize=self.config.cluster_binsize

        for ld in self.iter_region(bam,contig,start,end):
            ld_contig,ld_ref_start=ld.contig,ld.ref_start

            #TODO: Handle leads overlapping region ends (start/end)
            if contig==ld_contig and ld_ref_start >= start and ld_ref_start < end:
                pos_leadtab=int(ld_ref_start/ld_binsize)*ld_binsize
                self.record_lead(ld,pos_leadtab)
            else:
                externals.append(ld)

        if self.config.dev_cache:
            self.dev_store_leadtab(contig,start,end,externals)

        return externals

    def iter_region(self,bam,contig,start=None,end=None):
        leads_all=[]
        binsize=self.config.cluster_binsize
        coverage_binsize=self.config.coverage_binsize
        coverage_shift_bins=self.config.coverage_shift_bins
        long_ins_threshold=self.config.long_ins_length*0.5
        qc_nm=self.config.qc_nm
        phase=self.config.phase
        advanced_tags=qc_nm or phase
        mapq_min=self.config.mapq
        alen_min=self.config.min_alignment_length

        for read in bam.fetch(contig,start,end,until_eof=False):
            #if self.read_count % 1000000 == 0:
            #    gc.collect()
            if read.reference_start < start or read.reference_start >= end:
                continue

            self.read_id+=1
            self.read_count+=1

            if read.mapping_quality < mapq_min or read.is_secondary or read.query_alignment_length < alen_min:
                continue

            has_sa=read.has_tag("SA")
            use_clips=self.config.detect_large_ins and not read.is_supplementary and not has_sa

            nm=-1
            curr_read_id=self.read_id
            if advanced_tags:
                if qc_nm:
                    if read.has_tag("NM"):
                        nm=read.get_tag("NM")/float(read.query_alignment_length+1)

                if phase:
                    curr_read_id=(self.read_id,str(read.get_tag("HP")) if read.has_tag("HP") else "NULL",str(read.get_tag("PS")) if read.has_tag("PS") else "NULL")

            #Extract small indels
            for lead in read_iterindels(curr_read_id,read,contig,self.config,use_clips,read_nm=nm):
                yield lead

            #Extract read splits
            if has_sa:
                if read.is_supplementary:
                    for lead in read_itersplits_bnd(curr_read_id,read,contig,self.config,read_nm=nm):
                        yield lead
                else:
                    for lead in read_itersplits(curr_read_id,read,contig,self.config,read_nm=nm):
                        yield lead

            #Record in coverage table
            read_end=read.reference_start+read.reference_length
            assert(read_end==read.reference_end)
            #assert(read_end>=read.reference_start)
            if read.is_reverse:
                target_tab=self.covrtab_rev
            else:
                target_tab=self.covrtab_fwd
            covr_start_bin=(int(read.reference_start/coverage_binsize)+coverage_shift_bins)*coverage_binsize
            covr_end_bin=(int(read_end/coverage_binsize)-coverage_shift_bins)*coverage_binsize


            if covr_end_bin > covr_start_bin:
                self.covrtab_min_bin=min(self.covrtab_min_bin,covr_start_bin)
                target_tab[covr_start_bin]=target_tab[covr_start_bin]+1 if covr_start_bin in target_tab else 1

                if read_end <= self.end:
                    target_tab[covr_end_bin]=target_tab[covr_end_bin]-1 if covr_end_bin in target_tab else -1


    def dev_leadtab_filename(self,contig,start,end):
        scriptloc=os.path.dirname(os.path.realpath(sys.argv[0]))
        if self.config.dev_cache_dir==None:
            cache_dir=f"{scriptloc}/cache"
        else:
            cache_dir=self.config.dev_cache_dir
        return f"{cache_dir}/{os.path.basename(self.config.input)}_{contig}_{start}_{end}.pickle"

    def dev_store_leadtab(self,contig,start,end,externals):
        data={"externals":externals, "self": self}
        filename=self.dev_leadtab_filename(contig,start,end)
        with open(filename,"wb") as h:
            pickle.dump(data,h)
        print(f"(DEV/Cache) Dumped leadtab to {filename}")

    def dev_load_leadtab(self,contig,start,end):
        filename=self.dev_leadtab_filename(contig,start,end)

        if not os.path.exists(filename):
            return False

        with open(filename,"rb") as h:
            data=pickle.load(h)
        for item in data["self"].__dict__:
            self.__dict__[item]=data["self"].__dict__[item]
        print(f"(DEV/Cache) Loaded leadtab from {filename}")
        return data["externals"]

    def dev_debug_graph(self,title):
        import matplotlib.pyplot as plt
        import seaborn as sns
        print(title)
        sns.set()
        data=[]
        for k,v in self.leadtab.items():
            data.append(len(v))
            if len(data)>50000:
                break

        plt.hist(data,bins=[i for i in range(0,20)])
        #plt.savefig(filename)
        plt.title(title)
        plt.savefig(f"debug/{title}.png")
        plt.close()
