#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 09.09.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

from dataclasses import dataclass
import collections

from sniffles import util

@dataclass
class Group:
    """
    A sequence to which other sequences are pseudoaligned (mappings) through exact kmer matches (anchors)
    Can be reduced to a consensus sequence using the flatten() method
    """
    anchors: dict
    mappings: dict
    consensus_min: int
    consensus_max: int
    start_len: int
    origin: object
    support: list

    def flatten(self,frameshift_origin=False,return_debug=False):
        override_origin_min=3
        override_origin_indel_ratio=0.6

        consensus=""
        origin_space_consensus=""
        ins_consensus=""
        del_consensus=""
        mappings_consensus=""
        skip_until=self.consensus_min-1
        for ci in range(self.consensus_min,self.consensus_max+1,1):
            if not ci in self.mappings:
                origin_space_consensus+="-"
                ins_consensus+="-"
                del_consensus+="-"
                mappings_consensus+="-"
                continue
            type_fragments={"M":[], "I":[], "D":[]}
            type_fragments_counts={"M":0, "I":0, "D":0}
            for maptype, maplen, map_last_i, map_i, map_j, map_lead in self.mappings[ci]:
                type_fragments[maptype].append((maplen,map_lead.seq[map_last_i:map_i]))
                type_fragments_counts[maptype]+=1

            type_consensus=sorted(type_fragments, key=lambda t: len(type_fragments[t]))[-1]
            assert(len(type_fragments[type_consensus])>0)

            #Favor match over I/D when ratio is below 2

            if type_consensus=="I" and type_fragments_counts["M"]>0 and type_fragments_counts["I"] / float(type_fragments_counts["M"]+type_fragments_counts["I"]) < override_origin_indel_ratio:
                type_consensus="N"

            if type_consensus=="D" and type_fragments_counts["M"]>0 and type_fragments_counts["D"] / float(type_fragments_counts["M"]+type_fragments_counts["D"]) < override_origin_indel_ratio:
                type_consensus="N"

            if type_consensus=="M":
                #print(type_fragments_seq)
                if len(type_fragments[type_consensus]) < 3:
                    #Favor original sequence when only two mappings present
                    maplen,seq=type_fragments[type_consensus][0]
                    common_base,count=seq[0],1
                else:
                    common_base,count=collections.Counter(seq[0] for maplen,seq in type_fragments[type_consensus]).most_common(1)[0]

                if ci > skip_until:
                    consensus+=common_base
                    origin_space_consensus+=common_base
                else:
                    origin_space_consensus+="-"
                ins_consensus+="-"
                del_consensus+="-"
                mappings_consensus+=str(min(9,count))

            elif type_consensus=="N":
                #Favor original sequence when only two mappings present
                maplen,seq=type_fragments["M"][0]
                common_base,count=seq[0],1

                if ci > skip_until:
                    consensus+=common_base
                    origin_space_consensus+=common_base
                else:
                    origin_space_consensus+="-"
                ins_consensus+="-"
                del_consensus+="-"
                mappings_consensus+=str(min(9,count))

            elif frameshift_origin and type_consensus=="I":
                common_maplen,count=collections.Counter(maplen for maplen,seq in type_fragments[type_consensus]).most_common(1)[0]
                common_seqlen,count=collections.Counter(len(seq) for maplen,seq in type_fragments[type_consensus] if maplen==common_maplen).most_common(1)[0]
                origin_space_consensus+="^"
                ins_consensus+=str(min(9,maplen))
                del_consensus+="-"
                mappings_consensus+=str(min(9,count))

                selected=[seq for maplen,seq in type_fragments[type_consensus] if common_maplen==maplen and len(seq)==common_seqlen]

                for i in range(common_seqlen):
                    common_ins_base,_=collections.Counter(seq[i] for seq in selected).most_common(1)[0]
                    consensus+=common_ins_base

                skip_until=ci+common_seqlen+common_maplen-1

            elif frameshift_origin and type_consensus=="D":
                del_consensus=del_consensus[:-1]+str(min(9,maplen))

        if not return_debug:
            return consensus
        else:
            return origin_space_consensus,ins_consensus,del_consensus,mappings_consensus,consensus

def make_group(lead,klen,skip,repetitive_kmers):
    anchors={}
    for i, kmer in iter_kmers(lead.seq,klen,skip):
        if kmer in repetitive_kmers:
            continue
        if not kmer in anchors:
            anchors[kmer]=[i]
        else:
            anchors[kmer].append(i)

    return Group(anchors=anchors,
                 mappings={i:[("M", 1, i, i+1, i, lead)] for i in range(len(lead.seq))},
                 consensus_min=0,
                 consensus_max=len(lead.seq)-1,
                 start_len=len(lead.seq),
                 origin=lead,
                 support=[lead])

def iter_kmers(seq,klen,skip):
    for i in range(0,len(seq)-klen,skip):
        yield (i,seq[i:i+klen])

def from_leads(leads, initial_lead, max_grp=3, klen=15, skip=1, skip_repetitive=1, extend=False):
    #TODO: Handle DELs in flattening
    #TODO: Allow mapping across putative INDELs at sequence ends, allowing correcting them
    #TODO: Dynamically adapt kmer len threshold so short+long origs can be assembled at the same time, same for repetitive_limit
    #TODO: Ensure that kmers before and after last mapping are added as (secondary) anchors (hints 4 indels)

    repetitive_limit=len(leads)+2

    #Generate repetitive kmer list
    kmer_counts={}
    #for lead in leads:
    for i, kmer in iter_kmers(initial_lead.seq,klen=klen,skip=skip_repetitive):
        if not kmer in kmer_counts:
            kmer_counts[kmer]=1
        else:
            kmer_counts[kmer]+=1
    repetitive_kmers=set(kmer for kmer,count in kmer_counts.items() if count>repetitive_limit)

    #Create initial group
    leads=sorted(leads,key=lambda ld: len(ld.seq))
    groups=[make_group(initial_lead,klen,skip,repetitive_kmers)]

    for lead in leads:
        if lead.seq==None:
            continue

        assigned=False

        #Try to place this lead's INS seq in the group
        for group in groups[:max_grp]:
            mapped=False
            last_i,last_j,last_was_repeat=None,None,False
            first_i,first_j=None,None

            #Attempt to anchor kmers from seq to current group
            for i, kmer in iter_kmers(lead.seq,klen=klen,skip=skip):
                if not kmer in group.anchors: #and len(set(group.anchors[kmer]))==1:
                    continue

                if len(set(group.anchors[kmer])) > 1:
                    if last_j==None:
                        #Do not allow the first kmer to have multiple potential positions
                        continue

                    #Kmer occurs at multiple points in group sequence. Attempt to choose anchor most suitable for linear mapping
                    j=group.anchors[kmer][0]
                    for anchor_j_curr in group.anchors[kmer][1:]:
                        if abs(anchor_j_curr-(last_j+1)) < abs(j-(last_j+1)):
                            j=anchor_j_curr
                    last_was_repeat=True
                else:
                    #Kmer can be unambiguously placed on group sequence
                    j=group.anchors[kmer][0]

                if first_i==None and first_j==None:
                    first_i=i
                    first_j=j

                if last_i != None and last_j != None:
                    #Determine mapping type by position of current and last kmer on the lead sequence and group sequence
                    dist_lead=i-last_i
                    dist_group=j-last_j

                    if not last_was_repeat and dist_lead > dist_group:
                        maptype="D"
                        maplen=dist_group-dist_lead
                        if not last_j in group.mappings:
                            group.mappings[last_j]=[]
                        group.mappings[last_j].append((maptype, maplen, last_i, i, j, lead))

                    elif not last_was_repeat and dist_lead < dist_group:
                        maptype="I"
                        maplen=dist_group-dist_lead
                        if not last_j in group.mappings:
                            group.mappings[last_j]=[]
                        group.mappings[last_j].append((maptype, maplen, last_i, i, j, lead))

                    else:
                        maptype="M"
                        maplen=last_i-i
                        act_j=last_j
                        for act_i in range(last_i,i):
                            if not act_j in group.mappings:
                                group.mappings[act_j]=[]
                            group.mappings[act_j].append((maptype, maplen, act_i, act_i+1, last_j, lead))
                            act_j+=1

                mapped=True
                last_i=i
                last_j=j
                last_was_repeat=False

            if mapped:
                group.support.append(lead)
                assigned=True
                max_i=len(lead.seq)

                #Extend consensus group leftwards
                if extend and first_i - (first_j-group.consensus_min) > 0: # and first_j < 1.5*klen:
                    for i in range(first_i,-1,-1):
                        curr_j=first_j-(first_i-i)
                        if curr_j <= group.consensus_min: #+klen/1:
                            assert(curr_j==group.consensus_min)
                            kmer=lead.seq[i:i+klen]
                            if not kmer in repetitive_kmers:
                                group.anchors[kmer]=[curr_j]
                            group.consensus_min-=1
                            if not curr_j in group.mappings:
                                group.mappings[curr_j]=[]
                            group.mappings[curr_j].append(("M", 1, i, i+1, first_j, lead))

                #Extend consensus group rightwards
                if extend and last_j+(max_i-last_i) > group.consensus_max:
                    for i in range(last_i,max_i,1):
                        curr_j=last_j+(i-last_i)
                        if curr_j+klen >= group.consensus_max: #-klen/1:
                            assert(curr_j+klen==group.consensus_max)
                            kmer=lead.seq[i:i+klen]
                            if not kmer in repetitive_kmers:
                                group.anchors[kmer]=[curr_j]
                            group.consensus_max+=1
                            if not curr_j in group.mappings:
                                group.mappings[curr_j]=[]
                            group.mappings[curr_j].append(("M", 1, i, i+1, last_j, lead))
                break

        if not assigned:
            groups.append(make_group(lead,klen,skip,repetitive_kmers))
            #groups.sort(key=lambda g: len(g.mappings), reverse=True)

    return groups


def novel_from_reads(best_lead,other_leads,klen,skip,skip_repetitive,debug=False):
    consensus_min=2
    maxshift=klen
    minspan=0.2
    minalns=0.25
    minident=0.5
    minident_abs=5
    minbestdiff=3

    alignments=[]
    anchors={}
    taboo=set()
    for i, kmer in iter_kmers(best_lead.seq,klen=klen,skip=skip_repetitive):
        if kmer in taboo:
            continue
        if kmer in anchors:
            del anchors[kmer]
            taboo.add(kmer)
            continue
        anchors[kmer]=i

    for leadi,lead in enumerate(other_leads):
        last_i=None
        last_j=None
        conseq=""
        span=0
        for j, kmer in iter_kmers(lead.seq,klen=klen,skip=skip):
            if not kmer in anchors:
                continue
            i=anchors[kmer]
            if abs(i-j) > maxshift:
                continue
            if last_i != None and i <= last_i:
                continue

            if last_i == None:
                if j>0:
                    conseq="-"*i
            else:
                fwd_i=i-last_i
                fwd_j=j-last_j
                if len(conseq)+fwd_j > len(best_lead.seq):
                    fwd_j=len(best_lead.seq)-len(conseq)

                if fwd_i == fwd_j and fwd_j > 0:
                    span+=(j-last_j)
                    m=0
                    for l in range(1,(j-last_j)+1):
                        if lead.seq[last_j+l] == best_lead.seq[last_i+l]:
                            m+=1
                    ident=m/float((j-last_j))
                    if ident >= minident:
                        conseq+=lead.seq[last_j:j][:fwd_j]
                    else:
                        conseq+="-"*(fwd_j)
                else:
                    conseq+="-"*(fwd_j)
            last_i=i
            last_j=j

        if len(conseq) < len(best_lead.seq):
            conseq+="-"*(len(best_lead.seq)-len(conseq))

        conseq_new=[]
        h=0
        while h < len(best_lead.seq):
            if conseq[h]=="-":
                conseq_new.append("-")
                h+=1
            else:
                buffer=[]
                ident=0
                while h < len(best_lead.seq) and conseq[h]!="-":
                    ident+=(best_lead.seq[h]==conseq[h])
                    buffer.append(conseq[h])
                    h+=1
                if ident/float(len(buffer)) > minident and ident>minident_abs:
                    conseq_new.append("".join(buffer))
                else:
                    conseq_new.append("-"*len(buffer))
        conseq="".join(conseq_new)

        if span/float(len(best_lead.seq)) > minspan:
            alignments.append(conseq)

    maxal=1
    for i in range(len(best_lead.seq)):
        maxal=max(maxal,len([best_lead.seq[i]]+[a[i] for a in alignments if not a[i] in "^_"]))
    maxal=float(maxal)

    flattened=""
    for i in range(len(best_lead.seq)):
        al=[a[i] for a in alignments if not a[i]=="-"]
        if len(al) < consensus_min or len(al)/maxal < minalns:
            flattened+=best_lead.seq[i]
        else:
            top=util.most_common([best_lead.seq[i]]+al)
            if len(top)>1 and top[0][0]-top[1][0] >= minbestdiff:
                flattened+=top[0][1]
            else:
                flattened+=best_lead.seq[i]

            #if top[0][0]/float(len(al)+1) < 0.75:
            #    flattened+=best_lead.seq[i]
            #else:
            #    flattened+=top[0][1]

    #print("FLT",flattened)
    #if debug:
    #    print("B",best_lead.seq)
    #    for o in other_leads:
    #        print("O",o.seq)
    #    print("F",flattened)
    #    print("=====")
    return flattened
