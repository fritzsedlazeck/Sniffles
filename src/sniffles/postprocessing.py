#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 27.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass

import collections

from sniffles import sv
from sniffles import util
from sniffles import consensus

import math

def annotate_sv(svcall,config):
    genotype_sv(svcall,config)

    if config.phase:
        phase_sv(svcall,config)

    if svcall.svtype=="INS" and not config.symbolic:
        merged_leads=[l for l in svcall.postprocess.cluster.leads if l.seq!=None]

        if len(merged_leads):
            best_lead=merged_leads[0]
            best_index=0
            best_diff=abs(len(best_lead.seq) - svcall.svlen) + abs(best_lead.ref_start - svcall.pos)*1.5
            for i,ld in enumerate(merged_leads):
                if i==0:
                    continue
                curr_diff=abs(len(ld.seq) - svcall.svlen) + abs(ld.ref_start - svcall.pos)*1.5
                if curr_diff < best_diff:
                    best_lead=ld
                    best_index=i
                    best_diff=curr_diff

            merged_leads.pop(best_index)
            #merged_leads_new=list()

            #for lead in merged_leads:
            #    curr_lendiff=abs(len(ld.seq) - len(best_lead.seq)) + abs(ld.ref_start - best_lead.ref_start)
            #    #if curr_lendiff < 25:
            #    merged_leads_new.append(lead)
            #merged_leads=merged_leads_new


            if len(merged_leads) >= config.consensus_min_reads and not config.no_consensus:
                klen=config.consensus_kmer_len
                skip=config.consensus_kmer_skip_base+int(len(best_lead.seq)*config.consensus_kmer_skip_seqlen_mult)
                skip_repetitive=skip

                svcall.alt=consensus.novel_from_reads(best_lead,merged_leads,klen=klen,skip=skip,skip_repetitive=skip_repetitive)
            else:
                svcall.alt=best_lead.seq

def add_request(svcall,field,pos,requests_for_coverage,config):
    bin=int(pos/config.coverage_binsize)*config.coverage_binsize
    if not bin in requests_for_coverage:
        requests_for_coverage[bin]=[]
    requests_for_coverage[bin].append((svcall,field))

def coverage(calls,lead_provider,config):
    requests_for_coverage=coverage_build_requests(calls,lead_provider,config)
    return coverage_fulfill(requests_for_coverage,calls,lead_provider,config)

def coverage_build_requests(calls,lead_provider,config):
    requests_for_coverage={}
    for svcall in calls:
        start=svcall.pos
        if svcall.svtype=="INS":
            end=start+1
        else:
            end=svcall.pos+abs(svcall.svlen)
        #if start > end:
        #    assert(svcall.svtype=="DEL")
        #    start,end=end,start
        #add_request(svcall,"coverage_start",start+config.coverage_binsize,requests_for_coverage,config)
        add_request(svcall,"coverage_start",start-config.coverage_binsize,requests_for_coverage,config)
        add_request(svcall,"coverage_center",int((start+end)/2),requests_for_coverage,config)
        #add_request(svcall,"coverage_end",end-config.coverage_binsize,requests_for_coverage,config)
        add_request(svcall,"coverage_end",end+config.coverage_binsize,requests_for_coverage,config)
        add_request(svcall,"coverage_upstream",start-config.coverage_binsize*config.coverage_updown_bins,requests_for_coverage,config)
        add_request(svcall,"coverage_downstream",end+config.coverage_binsize*config.coverage_updown_bins,requests_for_coverage,config)
    return requests_for_coverage

def coverage_fulfill(requests_for_coverage,calls,lead_provider,config):
    if len(requests_for_coverage)==0:
        return -1,-1

    start_bin=lead_provider.covrtab_min_bin
    end_bin=int(lead_provider.end/config.coverage_binsize)*config.coverage_binsize
    coverage_fwd=0
    coverage_rev=0
    coverage_fwd_total=0
    coverage_rev_total=0
    n=0

    for bin in range(start_bin, end_bin+config.coverage_binsize,config.coverage_binsize):
        n+=1

        if bin in lead_provider.covrtab_fwd:
            coverage_fwd+=lead_provider.covrtab_fwd[bin]

        if bin in lead_provider.covrtab_rev:
            coverage_rev+=lead_provider.covrtab_rev[bin]

        if bin in requests_for_coverage:
            coverage_total_curr=coverage_fwd+coverage_rev
            for svcall, field in requests_for_coverage[bin]:
                setattr(svcall,field,coverage_total_curr)

        coverage_fwd_total+=coverage_fwd
        coverage_rev_total+=coverage_rev

    average_coverage_fwd=coverage_fwd_total/float(n) if n>0 else 0
    average_coverage_rev=coverage_rev_total/float(n) if n>0 else 0
    return average_coverage_fwd,average_coverage_rev

def qc_sv_support(svcall,coverage_global,config):
    if config.minsupport == "auto":
        if not qc_support_auto(svcall,coverage_global,config):
            return False
    else:
        if not qc_support_const(svcall,config):
            return False
    return True

def rescale_support(svcall,config):
    if svcall.svtype!="INS" or svcall.svlen < config.long_ins_length:
        return svcall.support
    else:
        base=svcall.support + svcall.get_info("SUPPORT_LONG")
        scale_factor=config.long_ins_rescale_mult*(float(svcall.svlen)/config.long_ins_length)
        return round(base*(config.long_ins_rescale_base+scale_factor))

def qc_support_auto(svcall,coverage_global,config):
    support=rescale_support(svcall,config)
    #if svcall.svtype=="INS":
    #    coverage_list=[svcall.coverage_center]
    #else:
    coverage_list=[svcall.coverage_upstream,svcall.coverage_downstream]
    coverage_list=[c for c in coverage_list if c!=None]
    if len(coverage_list)==0:
        coverage_regional=coverage_global
    else:
        coverage_regional=round(sum(coverage_list)/len(coverage_list))
        if coverage_regional==0:
            coverage_regional=coverage_global
    coverage=(coverage_regional*config.minsupport_auto_regional_coverage_weight+coverage_global*(1.0-config.minsupport_auto_regional_coverage_weight))
    min_support=round(config.minsupport_auto_base+config.minsupport_auto_mult*coverage)
    return support >= min_support
    #return True

def qc_support_const(svcall,config):
    #svcall.set_info("MINSUPPORT",config.minsupport)
    return svcall.support >= config.minsupport

def qc_sv(svcall,config):
    if config.qc_stdev:
        stdev_pos=svcall.get_info("STDEV_POS")
        if stdev_pos > config.qc_stdev_abs_max:
            return False
        if svcall.svtype!="BND" and stdev_pos / abs(svcall.svlen) > 2.0:
            return False

        stdev_len = svcall.get_info("STDEV_LEN")
        if stdev_len != None:
            if svcall.svtype != "BND" and stdev_len / abs(svcall.svlen) > 1.0:
                return False
            if stdev_len > stdev_len > config.qc_stdev_abs_max:
                return False

    #if (svcall.coverage_upstream != None and svcall.coverage_upstream < config.qc_coverage) or (svcall.coverage_downstream != None and svcall.coverage_downstream < config.qc_coverage):
    if svcall.svtype != "DEL" and svcall.svtype != "INS" and (svcall.coverage_center != None and svcall.coverage_center < config.qc_coverage):
        return False

    if svcall.svtype == "DEL" and abs(svcall.svlen) >= config.long_del_length and not config.non_germline:
        if svcall.coverage_center != None and svcall.coverage_upstream != None and svcall.coverage_downstream != None and svcall.coverage_center > (svcall.coverage_upstream+svcall.coverage_downstream)/2.0 * config.long_del_coverage:
            return False
    elif svcall.svtype=="INS" and ( (svcall.coverage_upstream != None and svcall.coverage_upstream < config.qc_coverage) or (svcall.coverage_downstream != None and svcall.coverage_downstream < config.qc_coverage)):
        return False
    elif svcall.svtype == "DUP" and abs(svcall.svlen) >= config.long_dup_length and not config.non_germline:
        if svcall.coverage_center != None and svcall.coverage_upstream != None and svcall.coverage_downstream != None and svcall.coverage_center < (svcall.coverage_upstream+svcall.coverage_downstream)/2.0 * config.long_dup_coverage:
            return False

    return True

def qc_sv_post_annotate(svcall,config):
    if (len(svcall.genotypes)==0 or (svcall.genotypes[0][0]!="." and svcall.genotypes[0][0]+svcall.genotypes[0][1]<2)) and (svcall.coverage_center != None and svcall.coverage_center < config.qc_coverage):
        return False

    return True

def binomial_coef(n,k):
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def binomial_probability(k,n,p):
    try:
        #Binomial coef cancels out for likelihood ratios
        #return binomial_coef(n,k) * (p**k) * ((1.0-p)**(n-k))
        return (p**k) * ((1.0-p)**(n-k))
    except OverflowError:
        return 1.0

def likelihood_ratio(q1,q2):
    if q1/q2>0:
        try:
            return math.log(q1/q2,10)
        except ValueError:
            return 0
    else:
        return 0

def genotype_sv(svcall,config):
    normalization_target=250
    hom_ref_p=config.genotype_error
    het_p=(1.0/config.genotype_ploidy) # - config.genotype_error
    hom_var_p=1.0 - config.genotype_error
    coverage=0

    #Count inline events only once per read, but split events as individual alignments, as in coverage calculation
    leads=svcall.postprocess.cluster.leads
    support=rescale_support(svcall,config)

    if svcall.svtype=="INS":
        coverage_list=[svcall.coverage_center]
    else:
        if svcall.svtype=="DUP":
            if False and svcall.coverage_start!=None and svcall.coverage_end!=None:
                if svcall.coverage_start>svcall.coverage_end:
                    coverage_list=[svcall.coverage_end]
                else:
                    coverage_list=[svcall.coverage_start]
            else:
                coverage_list=[svcall.coverage_start,svcall.coverage_end]
            coverage+=round(support*0.75)
        elif svcall.svtype=="INV":
            coverage_list=[svcall.coverage_upstream,svcall.coverage_downstream]
            coverage+=round(support*0.5)
        else:
            coverage_list=[svcall.coverage_start,svcall.coverage_center,svcall.coverage_end]

    coverage_list=[c for c in coverage_list if c!=None]
    if len(coverage_list)==0:
        return
    coverage+=round(sum(coverage_list)/len(coverage_list))

    if support > coverage:
        coverage=support

    af=support / float(coverage)

    genotype_p=[((0,0),hom_ref_p),
                ((0,1),het_p),
                ((1,1),hom_var_p)]

    max_lead=max(support,coverage)
    if max_lead>normalization_target:
        norm=normalization_target/float(max_lead)
        normalized_support=round(support*norm)
        normalized_coverage=round(coverage*norm)
    else:
        normalized_support=support
        normalized_coverage=coverage

    genotype_likelihoods=[]
    for gt, p in genotype_p:
        q=binomial_probability(normalized_support,normalized_coverage,p)
        genotype_likelihoods.append((gt,q))
    genotype_likelihoods.sort(key=lambda k: k[1], reverse=True)

    sum_likelihoods=sum(q for gt,q in genotype_likelihoods)
    normalized_likelihoods=[ (gt,(q/sum_likelihoods)) for gt,q in genotype_likelihoods]

    gt1,q1=normalized_likelihoods[0]
    gt2,q2=normalized_likelihoods[1]
    qz=[q for gt,q in normalized_likelihoods if gt==(0,0)][0]
    genotype_z_score = min(60,int((-10) * likelihood_ratio(qz,q1)))
    genotype_quality = min(60,int((-10) * likelihood_ratio(q2,q1)))

    is_long_ins=(svcall.svtype=="INS" and svcall.svlen >= config.long_ins_length)
    if genotype_z_score < config.genotype_min_z_score and not config.non_germline and not is_long_ins:
        svcall.filter="GT"

    if is_long_ins and gt1==(0,0):
        a,b=".","."
    else:
        a,b=gt1
    svcall.genotypes[0]=(a,b,genotype_quality,coverage-support,support) #f"{a}/{b}:{genotype_quality}:{coverage-support}:{support}"
    svcall.set_info("AF",af)

def phase_sv(svcall,config):
    reads_phases={lead.read_id[0]: (lead.read_id[1],lead.read_id[2]) for lead in svcall.postprocess.cluster.leads}
    hp_list=util.most_common(hp for hp,ps in reads_phases.values())
    ps_list=util.most_common(ps for hp,ps in reads_phases.values())

    hp_support,hp=hp_list[0]
    ps_support,ps=ps_list[0]
    if hp==None:
        hp="NULL"
    if ps==None:
        ps="NULL"

    other_hp_support=sum(other_supp for other_supp, other_hp in hp_list if other_hp != hp and other_hp != "NULL")
    other_ps_support=sum(other_supp for other_supp, other_ps in ps_list if other_ps != ps and other_ps != "NULL")

    hp_filter="FAIL"
    if hp != "NULL" and hp_support > 0 and float(other_hp_support)/(hp_support+other_hp_support) < config.phase_conflict_threshold:
        hp_filter="PASS"

    ps_filter="FAIL"
    if hp != "NULL" and ps_support > 0 and float(other_ps_support)/(ps_support+other_ps_support) < config.phase_conflict_threshold:
        ps_filter="PASS"

    svcall.set_info("PHASE",f"{hp},{ps},{hp_support},{ps_support},{hp_filter},{ps_filter}")
