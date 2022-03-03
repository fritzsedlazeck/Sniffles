#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 15.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass

from sniffles import util

TYPES=["INS","DEL","DUP","INV","BND"]

@dataclass
class SVCallBNDInfo:
    mate_contig: str
    mate_ref_start: int
    is_first: bool
    is_reverse:bool

@dataclass
class SVCallPostprocessingInfo:
    cluster: list

@dataclass
class SVCall:
    contig: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: int
    filter: str
    info: dict

    svtype: str
    svlen: int
    end: int
    genotypes: dict

    precise: bool
    support: int
    rnames: list

    qc: bool
    nm: float
    postprocess: SVCallPostprocessingInfo

    fwd: int = None
    rev: int = None

    coverage_upstream: int = None
    coverage_downstream: int = None
    coverage_start: int = None
    coverage_center: int = None
    coverage_end: int = None

    def set_info(self,k,v):
        self.info[k]=v

    def get_info(self,k):
        return self.info[k] if k in self.info else None

    def has_info(self,k):
        return k in self.info

    def finalize(self):
        self.postprocess=None

@dataclass
class SVGroup:
    candidates: list
    pos_mean: float
    len_mean: float
    included_samples: set
    coverages_nonincluded: dict

def calculate_bounds(svtype,ref_start_mode,svlen_mode):
    if svtype=="INS":
        svstart=ref_start_mode
        svend=ref_start_mode
    elif svtype=="DEL":
        svstart=ref_start_mode+svlen_mode
        svend=ref_start_mode
    else:
        svstart=ref_start_mode
        svend=svstart+abs(svlen_mode)
    return svstart,svend

def call_from(cluster,config,keep_qc_fails,task):
    leads=cluster.leads

    svtype=cluster.svtype
    stdev_pos,stdev_len=None,None
    qc=True

    svlen=util.center(v.svlen for v in leads)

    if abs(svlen) < config.minsvlen:
        return

    #Count inline events only once per read, but split events as individual alignments, as in coverage calculation
    #inline_qnames=set(k.read_qname for k in leads if k.source=="INLINE")
    #support=len(inline_qnames)+sum(1 for k in leads if k.source!="INLINE")
    if svtype=="INS" and svlen>=config.long_ins_length:
        support_long_set=set(lead.read_qname for lead in cluster.leads_long)
        support=len(set(k.read_qname for k in leads) | support_long_set)
        support_long=len(support_long_set)
    else:
        support=len(set(k.read_qname for k in leads))
        support_long=0
    ref_start=util.center(v.ref_start for v in leads)
    stdev_pos=util.stdev(util.trim((v.ref_start for v in leads)))
    if svtype!="BND":
        stdev_len=util.stdev(util.trim((v.svlen for v in leads)))
        precise=(stdev_pos+stdev_len < config.precise)
    else:
        precise=(stdev_pos < config.precise)
    svstart,svend=calculate_bounds(svtype,ref_start,svlen)
    qual=int(util.mean(v.mapq for v in leads))

    support_fwd=sum(lead.strand == "+" for lead in leads)
    support_rev=len(leads) - support_fwd

    if config.qc_strand and (support_fwd==0 or support_rev==0):
        qc=False

    if config.qc_nm:
        nm_mean=util.mean(v.nm for v in leads)
        if nm_mean > config.qc_nm_max:
            qc=False
    else:
        nm_mean=-1

    if not keep_qc_fails and not qc:
        return

    svpi=SVCallPostprocessingInfo(cluster=cluster)

    if config.output_rnames or config.snf!=None:
        rnames=list(set(k.read_qname for k in leads))
    else:
        rnames=None

    svcall=SVCall(contig=cluster.contig,
                  pos=svstart,
                  id=f"{svtype}.{task.sv_id:X}S{task.id:X}",
                  ref="N",
                  alt=f"<{svtype}>",
                  qual=qual,
                  filter="PASS",
                  info=dict(),
                  svtype=svtype,
                  svlen=svlen,
                  end=svend,
                  genotypes=dict(),
                  precise=precise,
                  support=support,
                  rnames=rnames,
                  postprocess=svpi,
                  qc=qc,
                  nm=nm_mean,
                  fwd=support_fwd,
                  rev=support_rev)

    if svtype=="BND":
        resolve_bnd(svcall,cluster,config)
    elif svtype=="INS":
        svcall.set_info("SUPPORT_LONG",support_long)

    if stdev_pos!=None:
        svcall.set_info("STDEV_POS",stdev_pos)
    if stdev_len!=None:
        svcall.set_info("STDEV_LEN",stdev_len)

    task.sv_id+=1

    yield svcall

def merge_inner_bounds(leads,config):
    svlength_binsize=25
    read_svlengths={}
    read_starts={}
    for ld in leads:
        if not ld.read_id in read_svlengths:
            read_svlengths[ld.read_id]=0
            read_starts[ld.read_id]=ld.ref_start
        read_svlengths[ld.read_id]+=ld.svlen
        read_starts[ld.read_id]=min(ld.ref_start,read_starts[ld.read_id])
    svlen=util.center(v for k,v in read_svlengths.items())
    pos=util.center(v for k,v in read_starts.items())

    svlengths_reads={}
    for read_id, length in read_svlengths.items():
        bin=int(length/svlength_binsize)*svlength_binsize
        if not bin in svlengths_reads:
            svlengths_reads[bin]=[(read_id,length)]
        else:
            svlengths_reads[bin].append((read_id,length))

    return pos,svlen,util.stdev(util.trim((v for k,v in read_starts.items()))),util.stdev(util.trim((v for k,v in read_svlengths.items())))

def resolve_bnd(svcall,cluster,config):
    mate_contig=util.most_common_top([lead.bnd_info.mate_contig for lead in cluster.leads])
    selected=[lead for lead in cluster.leads if lead.bnd_info.mate_contig==mate_contig]
    mate_ref_start=util.center([lead.bnd_info.mate_ref_start for lead in selected])
    is_first=util.most_common_top([lead.bnd_info.is_first for lead in selected])
    is_reverse=util.most_common_top([lead.bnd_info.is_reverse for lead in selected])
    svcall.alt=(("N" if is_first else "") +
                ("]" if is_reverse else "[" ) +
                f"{mate_contig}:{mate_ref_start}" +
                ("]" if is_reverse else "[" ) +
                ("N" if not is_first else ""))
    svcall.support=len(set(k.read_qname for k in selected))
    cluster.leads=selected
    svcall.bnd_info=SVCallBNDInfo(mate_contig=mate_contig, mate_ref_start=mate_ref_start, is_first=is_first, is_reverse=is_reverse)
    svcall.set_info("CHR2",mate_contig)

def call_groups(svgroups,config,task):
    for group in svgroups:
        svcall=call_group(group,config,task)
        if svcall!=None:
            yield svcall

def call_group(svgroup,config,task):
    """For combining grouped SV calls from multiple .snf samples into one SV call"""
    first_cand=svgroup.candidates[0]

    #Filtering
    samples_count=float(len(config.snf_input_info))
    sample_internal_ids=set(sample["internal_id"] for sample in config.snf_input_info)
    total_count=len(svgroup.included_samples)
    pass_count=sum(cand.qc==True for cand in svgroup.candidates)
    qc=(pass_count > 0 and pass_count/samples_count >= config.combine_high_confidence) or (total_count/samples_count >= config.combine_low_confidence and total_count >= config.combine_low_confidence_abs)

    if not qc:
        return None

    if (not config.combine_output_filtered) and not any(cand.qc and cand.filter=="PASS" for cand in svgroup.candidates):
        return None

    if config.output_rnames:
        rnames=[]
    else:
        rnames=None

    genotypes={}
    genotyped_count=0
    for cand in svgroup.candidates:
        if config.output_rnames and cand.rnames!=None:
            rnames.extend(cand.rnames)
        if not 0 in cand.genotypes:
            cand.genotypes[0]=(".",".",0,0,cand.support)
        if cand.sample_internal_id in genotypes:
            #Intra-sample merging
            a,b,gt_qual,dr,dv=cand.genotypes[0]
            curr_a,curr_b,curr_gt_qual,curr_dr,curr_dv,curr_id=genotypes[cand.sample_internal_id]
            new_id=curr_id+","+config.id_prefix+cand.id
            if (curr_a==".") or (a != "." and (a,b) >= (curr_a,curr_b)):
                genotypes[cand.sample_internal_id]=(a,b,gt_qual,dr,dv,new_id)
            else:
                genotypes[cand.sample_internal_id]=(curr_a,curr_b,curr_gt_qual,curr_dr,curr_dv,new_id)
        else:
            a,b,gt_qual,dr,dv=cand.genotypes[0]
            genotypes[cand.sample_internal_id]=(a,b,gt_qual,dr,dv,config.id_prefix+cand.id)
        genotyped_count+=1

    for sample_internal_id in sample_internal_ids:
        if sample_internal_id in genotypes:
            continue
        coverage=svgroup.coverages_nonincluded[sample_internal_id]
        if coverage >= config.combine_null_min_coverage:
            genotypes[sample_internal_id]=(0,0,0,coverage,0,"NULL")
        else:
            genotypes[sample_internal_id]=(".",".",0,coverage,0,"NULL")

    if config.combine_consensus:
        genotypes_consensus={}
        for a,b,gt_qual,dr,dv in genotypes.values():
            if not (a,b) in genotypes_consensus:
                genotypes_consensus[(a,b)]={"count":0,"qual":list(),"dr":list(),"dv":list()}
            genotypes_consensus[(a,b)]["count"]+=1
            genotypes_consensus[(a,b)]["qual"].append(gt_qual)
            genotypes_consensus[(a,b)]["dr"].append(dr)
            genotypes_consensus[(a,b)]["dv"].append(dv)
        most_common_count=genotypes_consensus[sorted(genotypes_consensus,key=lambda k: genotypes_consensus[k]["count"],reverse=True)[0]]["count"]
        most_common_gt=[gt for gt in genotypes_consensus if genotypes_consensus[gt]["count"]==most_common_count]
        cons_a,cons_b=max(most_common_gt)
        consensus_info=genotypes_consensus[(cons_a,cons_b)]
        genotypes={0:(cons_a,cons_b,int(sum(consensus_info["qual"])/consensus_info["count"]),sum(consensus_info["dr"]),sum(consensus_info["dv"]))}
        if cons_a!=1 and cons_b!=1:
            return None

    svcall_pos=int(util.median(cand.pos for cand in svgroup.candidates))
    svcall_svlen=int(util.median(cand.svlen for cand in svgroup.candidates))
    svcall_end=svcall_pos+abs(svcall_svlen)
    svcall_alt=first_cand.alt
    svcall_alt_mindist=abs(len(svcall_alt)-svcall_svlen)
    if first_cand.svtype=="INS":
        for cand in svgroup.candidates:
            dist=abs(len(cand.alt)-svcall_svlen)
            if dist < svcall_alt_mindist:
                svcall_alt_mindist=dist
                svcall_alt=cand.alt

    svcall=SVCall(contig=first_cand.contig,
                  pos=svcall_pos,
                  id=f"{first_cand.svtype}.{task.sv_id:X}M{task.id:X}",
                  ref="N",
                  alt=svcall_alt,
                  qual=round(util.mean(int(cand.qual) for cand in svgroup.candidates)),
                  filter="PASS",
                  info=dict(),
                  svtype=first_cand.svtype,
                  svlen=svcall_svlen,
                  end=svcall_end,
                  genotypes=genotypes,
                  precise=sum(int(cand.precise) for cand in svgroup.candidates)/float(len(svgroup.candidates)) > 0.5,
                  support=round(util.mean(cand.support for cand in svgroup.candidates)),
                  rnames=rnames,
                  postprocess=None,
                  qc=True,
                  nm=-1,
                  fwd=sum(cand.fwd for cand in svgroup.candidates),
                  rev=sum(cand.rev for cand in svgroup.candidates),
                  coverage_upstream=util.mean_or_none_round(cand.coverage_upstream for cand in svgroup.candidates if cand.coverage_upstream!=None),
                  coverage_start=util.mean_or_none_round(cand.coverage_start for cand in svgroup.candidates if cand.coverage_start!=None),
                  coverage_center=util.mean_or_none_round(cand.coverage_center for cand in svgroup.candidates if cand.coverage_center!=None),
                  coverage_end=util.mean_or_none_round(cand.coverage_end for cand in svgroup.candidates if cand.coverage_end!=None),
                  coverage_downstream=util.mean_or_none_round(cand.coverage_downstream for cand in svgroup.candidates if cand.coverage_downstream!=None) )

    svcall.set_info("STDEV_POS",util.stdev(cand.pos for cand in svgroup.candidates))
    svcall.set_info("STDEV_LEN",util.stdev(cand.svlen for cand in svgroup.candidates))

    if abs(svcall.svlen) < config.minsvlen:
        return None

    task.sv_id+=1

    return svcall

def classify_splits(read,leads,config,main_contig):
    minsvlen_screen=config.minsvlen_screen
    maxsvlen_other=minsvlen_screen*5

    leads.sort(key=lambda ld: ld.qry_start)
    last=leads[0]
    last.svtypes_starts_lens=[]

    if last.qry_start >= config.long_ins_length*0.5:
        last.svtypes_starts_lens.append(("INS",last.ref_start,None))

    for i in range(1,len(leads)):
        curr=leads[i]
        curr.svtypes_starts_lens=[]
        qry_dist_abs=abs(curr.qry_start - last.qry_end)

        if curr.contig == last.contig:
            rev=(curr.strand == "-")
            fwd=not rev
            if curr.strand == last.strand:
                #
                #INS, DEL, DUP
                #
                if (fwd and (curr.qry_start - last.qry_end) >= minsvlen_screen
                    and (curr.ref_start - last.ref_end) < maxsvlen_other
                    and (curr.qry_start - last.qry_end) - (curr.ref_start - last.ref_end) >= minsvlen_screen):
                    #INS, FWD
                    svstart=curr.ref_start
                    svlen=(curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq=read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq=None
                    curr.svtypes_starts_lens.append(("INS",svstart,svlen))

                elif (rev and (curr.qry_start - last.qry_end) >= minsvlen_screen
                      and (last.ref_start - curr.ref_end) < maxsvlen_other
                      and (curr.qry_start - last.qry_end) - (last.ref_start - curr.ref_end) >= minsvlen_screen):
                    #INS, REV
                    svstart=last.ref_start
                    svlen=(curr.qry_start - last.qry_end)
                    if svlen <= config.dev_seq_cache_maxlen:
                        curr.seq=read.query_sequence[last.qry_end:curr.qry_start]
                    else:
                        curr.seq=None
                    curr.svtypes_starts_lens.append(("INS",svstart,svlen))

                elif (fwd and (curr.ref_start - last.ref_end) >= minsvlen_screen
                      and qry_dist_abs < maxsvlen_other
                      and (curr.ref_start - last.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                        #DEL, FWD
                        svstart=curr.ref_start
                        svlen=(curr.ref_start - last.ref_end)
                        curr.svtypes_starts_lens.append(("DEL",svstart,-svlen))

                elif (rev and (last.ref_start - curr.ref_end) >= minsvlen_screen
                      and qry_dist_abs < maxsvlen_other
                      and (last.ref_start - curr.ref_end) - (curr.qry_start - last.qry_end) >= minsvlen_screen):
                        #DEL, REV
                        svstart=last.ref_start
                        svlen=(last.ref_start - curr.ref_end)
                        curr.svtypes_starts_lens.append(("DEL",svstart,-svlen))

                elif fwd and curr.ref_start <= last.ref_end:
                    if qry_dist_abs < maxsvlen_other:
                        #DUP, FWD
                        svstart=curr.ref_start
                        svlen=(last.ref_end - curr.ref_start)
                        if svlen >= minsvlen_screen:
                            curr.svtypes_starts_lens.append(("DUP",svstart,svlen))

                elif rev and last.ref_start <= curr.ref_end:
                    if qry_dist_abs < maxsvlen_other:
                        #DUP, REV
                        svstart=last.ref_start
                        svlen=(curr.ref_end - last.ref_start)
                        if svlen >= minsvlen_screen:
                            curr.svtypes_starts_lens.append(("DUP",svstart,svlen))

            elif qry_dist_abs < maxsvlen_other:
                #
                #INV
                #
                if fwd and curr.ref_start <= last.ref_start:
                    #CASE B
                    svstart=curr.ref_start
                    svlen=last.ref_start-curr.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif fwd and curr.ref_start > last.ref_start:
                    #CASE C
                    svstart=last.ref_start
                    svlen=curr.ref_start-last.ref_start
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif rev and curr.ref_end >= last.ref_end:
                    #CASE A
                    svstart=last.ref_end
                    svlen=curr.ref_end-last.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

                elif rev and curr.ref_end < last.ref_end:
                    #CASE D
                    svstart=curr.ref_end
                    svlen=last.ref_end-curr.ref_end
                    if svlen >= minsvlen_screen:
                        curr.svtypes_starts_lens.append(("INV",svstart,svlen))

        elif qry_dist_abs < maxsvlen_other:
            #
            #BND
            #
            if curr.contig == main_contig:
                a,b=curr,last
            else:
                a,b=last,curr

            if a.contig == main_contig:
                is_first=a.qry_start < b.qry_start
                if is_first:
                    if a.strand=="+":
                        svstart=a.ref_end
                    else:
                        svstart=a.ref_start
                else:
                    if a.strand=="+":
                        svstart=a.ref_start
                    else:
                        svstart=a.ref_end
                a.svtypes_starts_lens.append(("BND",
                                            svstart,
                                            SVCallBNDInfo(b.contig,
                                                          b.ref_start,
                                                          is_first,
                                                          a.strand!=b.strand)))
        last=curr
