#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 27.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

import statistics
import math
import sys
import time

class Sniffles2Exit(Exception):
    pass

def stdev(nums):
    nums=list(nums)
    return statistics.stdev(nums) if len(nums) > 1 else 0

def median(nums):
    return int(statistics.median(nums))

def median_or_mode(nums):
    nums=list(nums)
    top=most_common(nums)
    if len(top)>1 and (top[0][0]-top[1][0]<2):
    #if len(nums)>3:
        return median_noavg(nums)
    else:
        return median_modes(nums)

def median_noavg(nums):
    nums=sorted(list(nums))
    mid=int(len(nums)/2)
    return nums[mid]

def median_modes(nums):
    max_count=0
    counts={}
    for n in nums:
        if not n in counts:
            counts[n]=1
        else:
            counts[n]+=1
        max_count=max(max_count,counts[n])
    return median_noavg(k for k,n in counts.items() if max_count-n<3 )

def mean(nums):
    nums=list(nums)
    return sum(nums)/len(nums)

def mean_or_none(nums):
    nums=list(nums)
    if len(nums)==0:
        return None
    else:
        return sum(nums)/len(nums)

def mean_or_none_round(nums):
    r=mean_or_none(nums)
    if r==None:
        return r
    else:
        return round(r)

def trim(nums,pct=25):
    nums=sorted(list(nums))
    trim_n=int(len(nums)/float(100.0)*pct)
    if trim_n > 0:
        return nums[trim_n:-trim_n]
    else:
        return nums

def most_common(nums):
    counts={}
    for n in nums:
        if not n in counts:
            counts[n]=1
        else:
            counts[n]+=1
    return sorted(((counts[n],n) for n in counts), reverse=True)

def most_common_top(nums):
    result=most_common(nums)
    return sorted(item for count,item in result if count==result[0][0])[0]

def error(msg):
    sys.stderr.write("Sniffles2 Error: "+msg+"\n")
    sys.stderr.flush()

def fatal_error(msg):
    error(msg+" (Fatal error, exiting.)")
    exit(1)

def fatal_error_main(msg):
    error(msg+" (Fatal error, exiting.)")
    raise Sniffles2Exit

def load_tandem_repeats(filename,padding):
    contigs_tr={}
    unsorted=False
    with open(filename,"r") as handle:
        for line in handle:
            parts=line.split("\t")
            if len(parts) >= 3:
                contig,start,end=parts[:3]
                start=int(start)
                end=int(end)
                if not contig in contigs_tr:
                    contigs_tr[contig]=[]
                if len(contigs_tr[contig])>0:
                    last_start,last_end=contigs_tr[contig][-1]
                    if start < last_start:
                        unsorted=True
                contigs_tr[contig].append((max(0,int(start)-padding),int(end)+padding))

    if unsorted:
        print("Info: The tandem repeat annotations file was not sorted. Sorting it in-memory after loading... (please sort the .bed file once before to save time when running multiple samples)")
        sort_start=time.time()
        for contig in contigs_tr:
            contigs_tr[contig].sort()
        print(f"Info: Optional sorting of input tandem repeat annotations took {time.time()-sort_start:.2f}s.")

    return contigs_tr

center=median_modes
