#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 29.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass
import statistics
import math

from sniffles import sv
from sniffles import leadprov

@dataclass
class Cluster:
    id: str
    svtype: str
    contig: str
    start: int
    end: int
    seed: int
    leads: list
    repeat: bool
    leads_long: list

    def compute_metrics(self,max_n=100):
        self.span=self.end-self.start

        n=min(len(self.leads),max_n)
        if n==0:
            self.mean_svlen=0
            self.stdev_start=0
            return

        step=int(len(self.leads)/n)
        if n>1:
            self.mean_svlen=sum(self.leads[i].svlen for i in range(0,len(self.leads),step))/float(n)
            self.stdev_start=statistics.stdev(self.leads[i].ref_start for i in range(0,len(self.leads),step))
        else:
            self.mean_svlen=self.leads[0].svlen
            self.stdev_start=0

def merge_inner(cluster,threshold):
    read_seq={}
    for ld in cluster.leads:
        if not ld.read_qname in read_seq:
            read_seq[ld.read_qname]=[]
        read_seq[ld.read_qname].append(ld)

    cluster.leads=[]
    for qname in read_seq:
        read_seq[qname].sort(key=lambda k: k.ref_start)
        to_merge=read_seq[qname][0]

        curr_lead=to_merge

        last_ref_end=to_merge.ref_end
        last_qry_end=to_merge.qry_end
        last_ref_start=to_merge.ref_start
        last_qry_start=to_merge.qry_start

        for to_merge in read_seq[qname][1:]:
            merge=(threshold==-1) or ((abs(to_merge.ref_start-last_ref_end) < threshold or abs(to_merge.ref_start-last_ref_start) < threshold) and (abs(to_merge.qry_start-last_qry_end) < threshold or abs(to_merge.qry_start-last_qry_start) < threshold))
            if merge:
                curr_lead.svlen+=to_merge.svlen
                if to_merge.seq==None or curr_lead.seq==None:
                    curr_lead.seq=None
                else:
                    curr_lead.seq+=to_merge.seq
            else:
                cluster.leads.append(curr_lead)
                curr_lead=to_merge
            last_ref_end=to_merge.ref_end
            last_qry_end=to_merge.qry_end
            last_ref_start=to_merge.ref_start
            last_qry_start=to_merge.qry_start

        cluster.leads.append(curr_lead)
    return cluster

def resplit(cluster,prop,binsize,merge_threshold_min,merge_threshold_frac):
    bins_leads={}
    for lead in cluster.leads:
        bin=int(abs(prop(lead))/binsize)*binsize
        if not bin in bins_leads:
            bins_leads[bin]=[lead]
        else:
            bins_leads[bin].append(lead)

    new_clusters=list(sorted(bins_leads.keys()))
    i=1
    while len(new_clusters) > 1 and i < len(new_clusters):
        last_cluster=new_clusters[i-1]
        curr_cluster=new_clusters[i]
        merge_threshold=max(merge_threshold_min,min(curr_cluster,last_cluster)*merge_threshold_frac)
        merge=abs(curr_cluster - last_cluster) <= merge_threshold
        if merge:
            bins_leads[new_clusters[i]].extend(bins_leads[new_clusters[i-1]])
            new_clusters.pop(i-1)
            i=max(0,i-2)
        else:
            i+=1

    for cluster_index in new_clusters:
        new_cluster=Cluster(id=cluster.id+f".{cluster_index}",
                            svtype=cluster.svtype,
                            contig=cluster.contig,
                            start=cluster.start,
                            end=cluster.end,
                            seed=cluster.seed,
                            leads=bins_leads[cluster_index],
                            repeat=cluster.repeat,
                            leads_long=cluster.leads_long)
        yield new_cluster

def resplit_bnd(cluster,merge_threshold):
    contigs_leads={}

    #Split by mate contig and mate ref start
    for lead in cluster.leads:
        if not lead.bnd_info.mate_contig in contigs_leads:
            contigs_leads[lead.bnd_info.mate_contig]={}

        bin=int(lead.bnd_info.mate_ref_start/merge_threshold)*merge_threshold if merge_threshold>0 else 0
        if not bin in contigs_leads[lead.bnd_info.mate_contig]:
            contigs_leads[lead.bnd_info.mate_contig][bin]=[lead]
        else:
            contigs_leads[lead.bnd_info.mate_contig][bin].append(lead)

    for contig in contigs_leads:
        bins=sorted(contigs_leads[contig])
        curr_leads=[]+contigs_leads[contig][bins[0]]
        last_bin=bins[0]
        for bin in bins[1:]:
            if bin - last_bin <= merge_threshold:
                curr_leads.extend(contigs_leads[contig][bin])
            else:
                if len(curr_leads):
                    new_cluster=Cluster(id=cluster.id+f".CHR2.{contig}.POS2.{bin}",
                                        svtype=cluster.svtype,
                                        contig=cluster.contig,
                                        start=cluster.start,
                                        end=cluster.end,
                                        seed=cluster.seed,
                                        leads=[k for k in curr_leads],
                                        repeat=cluster.repeat,
                                        leads_long=None)
                    yield new_cluster
                curr_leads=[]+contigs_leads[contig][bin]
            last_bin=bin
        if len(curr_leads):
            new_cluster=Cluster(id=cluster.id+f".CHR2.{contig}.POS2.{bin}",
                                svtype=cluster.svtype,
                                contig=cluster.contig,
                                start=cluster.start,
                                end=cluster.end,
                                seed=cluster.seed,
                                leads=[k for k in curr_leads],
                                repeat=cluster.repeat,
                                leads_long=None)
            curr_leads=[]
            yield new_cluster

def resolve(svtype,leadtab_provider,config,tr):
    leadtab=leadtab_provider.leadtab[svtype]
    seeds=sorted(leadtab_provider.leadtab[svtype])

    if len(seeds)==0:
        return []

    #Initialize tandem repeat region handling
    if tr != None:
        tr_index=0
        if len(tr)==0:
            tr=None
        else:
            tr_start,tr_end=tr[tr_index]

    clusters=[]
    for seed_index,seed in enumerate(seeds):
        if config.dev_call_region != None:
            if seed<config.dev_call_region["start"] or seed>config.dev_call_region["end"]:
                continue

        within_tr=False
        if tr!=None and tr_index < len(tr):
            while tr_end < seed and tr_index+1 < len(tr):
                tr_index+=1
                tr_start,tr_end=tr[tr_index]
            if seed > tr_start and seed < tr_end:
                within_tr=True
            #print(within_tr,seed,tr_index,len(tr))
        #within_tr=True

        if svtype=="INS":
            leads=[lead for lead in leadtab[seed] if lead.svlen!=None]
            leads_long=[lead for lead in leadtab[seed] if lead.svlen==None]
        else:
            leads=leadtab[seed]
            leads_long=None

        cluster=Cluster(id=f"CL.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{seed_index}",
                        svtype=svtype,
                        contig=leadtab_provider.contig,
                        start=seed,
                        end=seed+config.cluster_binsize,
                        seed=seed,
                        leads=leads,
                        repeat=within_tr or config.repeat,
                        leads_long=leads_long)

        cluster.compute_metrics()
        clusters.append(cluster)


    #Merge clusters
    cluster_count_initial=len(clusters)
    i=0
    while i < len(clusters) - 1:
        curr_cluster=clusters[i]
        next_cluster=clusters[i+1]

        #assert((next_cluster.end - curr_cluster.start) >= 0)
        #assert((next_cluster.start - curr_cluster.end) >= 0)

        #Fast Adaptive Scanning Clustering
        #  Inter-Cluster distance < (StdDev-Start-Within-A + StdDev-Stat-Within-B) * r   (r=1.0 default)
        #  Merged size (OR: Inter-cluster distance) < (MeanSVlen between both clusters) * h (h=0.5 default)
        #  Max. overall cluster size criterion?

        inner_dist=(next_cluster.start - curr_cluster.end)
        outer_dist=(next_cluster.end - curr_cluster.start)
        merge = inner_dist <= min(curr_cluster.stdev_start, next_cluster.stdev_start) * config.cluster_r
        merge = merge or ( (config.repeat or curr_cluster.repeat or next_cluster.repeat) and outer_dist <= min(config.cluster_repeat_h_max, (abs(curr_cluster.mean_svlen)+abs(next_cluster.mean_svlen)) * config.cluster_repeat_h) )
        merge = merge or (svtype=="BND" and inner_dist <= config.cluster_merge_bnd)

        if merge:
            clusters.pop(i+1)
            curr_cluster.leads+=next_cluster.leads
            if svtype=="INS":
                curr_cluster.leads_long+=next_cluster.leads_long
            curr_cluster.end=next_cluster.end
            curr_cluster.repeat=curr_cluster.repeat or next_cluster.repeat
            curr_cluster.compute_metrics()
            i=max(0,i-2)
        i+=1

    if config.dev_dump_clusters:
        filename=f"{config.input}.clusters.{svtype}.{leadtab_provider.contig}.{leadtab_provider.start}.{leadtab_provider.end}.bed"
        print(f"Dumping clusters to {filename}")
        with open(filename,"w") as h:
            for c in clusters:
                info=f"ID={c.id}, #LEADS={len(c.leads)}; "
                for ld in c.leads:
                    info+=f"(ref_start={ld.ref_start},svlen={ld.svlen},source={ld.source}); "
                h.write(f"{c.contig}\t{c.start}\t{c.end}\t\"{info}\"\n")

    for cluster in clusters:
        if len(cluster.leads)==0:
            continue

        if svtype == "BND":
            if config.dev_no_resplit:
                yield cluster
            else:
                for new_cluster in resplit_bnd(cluster,merge_threshold=config.bnd_cluster_resplit):
                    yield new_cluster
        else:
            if svtype=="INS" or svtype=="DEL":
                if cluster.repeat:
                    merge_inner_threshold=-1
                else:
                    merge_inner_threshold=config.cluster_merge_pos

                merge_inner(cluster,merge_inner_threshold)

            if not cluster.repeat and not config.dev_no_resplit:
                for new_cluster in resplit(cluster,
                                           prop=lambda lead: lead.svlen,
                                           binsize=config.resplit_binsize,
                                           merge_threshold_min=config.minsvlen,
                                           merge_threshold_frac=config.cluster_merge_len):
                    yield new_cluster
            else:
                yield cluster

def resolve_block_groups(svtype,svcands,groups_initial,config):
    """For clustering groups of SVs for combining .snfs (multi-call)"""

    #TODO: Remove sorting
    groups=groups_initial
    for svcand in sorted(svcands,key=lambda cand: cand.support,reverse=True):
        best_group=None
        best_dist=math.inf

        if svtype=="BND":
            mate_contig,mate_ref_start=svcand.bnd_info.mate_contig,svcand.bnd_info.mate_ref_start
            for group in groups:
                #TODO: Favor bigger groups in placement
                dist=abs(group.pos_mean - svcand.pos) + abs(group.bnd_mate_ref_start_mean - mate_ref_start)
                if dist < best_dist and dist <= config.cluster_merge_bnd*2 and group.bnd_mate_contig==mate_contig:
                    if not config.combine_separate_intra or not svcand.sample_internal_id in group.included_samples:
                        best_group=group
                        best_dist=dist
        else:
            for group in groups:
                #TODO: Favor bigger groups in placement
                dist=abs(group.pos_mean - svcand.pos) + abs(abs(group.len_mean) - abs(svcand.svlen))
                minlen=float(min(abs(group.len_mean),abs(svcand.svlen)))
                if minlen>0 and dist < best_dist and dist <= config.combine_match * math.sqrt(minlen):
                    if not config.combine_separate_intra or not svcand.sample_internal_id in group.included_samples:
                        best_group=group
                        best_dist=dist

        if best_group==None:
            group=sv.SVGroup(candidates=[svcand],
                             pos_mean=float(svcand.pos),
                             len_mean=float(abs(svcand.svlen)),
                             included_samples=set([svcand.sample_internal_id]),
                             coverages_nonincluded=dict())
            if svtype=="BND":
                group.bnd_mate_contig=mate_contig
                group.bnd_mate_ref_start_mean=mate_ref_start
            groups.append(group)
        else:
            group_size=len(best_group.candidates)
            best_group.pos_mean*=group_size
            best_group.len_mean*=group_size
            best_group.pos_mean+=svcand.pos
            best_group.len_mean+=abs(svcand.svlen)
            if svtype=="BND":
                best_group.bnd_mate_ref_start_mean*=group_size
                best_group.bnd_mate_ref_start_mean+=mate_ref_start

            best_group.candidates.append(svcand)
            group_size+=1
            best_group.pos_mean/=group_size
            best_group.len_mean/=group_size
            best_group.included_samples.add(svcand.sample_internal_id)

            if svtype=="BND":
                best_group.bnd_mate_ref_start_mean/=group_size

    return groups
