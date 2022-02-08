#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 15.08.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

import pysam
import os

from sniffles import sv
from sniffles import util

def format_info(k,v):
    if isinstance(v,float):
        return f"{k}={v:.3f}"
    elif isinstance(v,list):
        return f"{k}={','.join(v)}"
    else:
        return f"{k}={v}"

def format_genotype(gt):
    if len(gt)==5:
        a,b,qual,dr,dv=gt
        return f"{a}/{b}:{qual}:{dr}:{dv}"
    else:
        a,b,qual,dr,dv,id=gt
        return f"{a}/{b}:{qual}:{dr}:{dv}:{id}"

class VCF:
    def __init__(self,config,handle):
        self.config=config
        self.handle=handle
        self.call_count=0
        self.info_order=["SVTYPE","SVLEN","END","SUPPORT","RNAMES","COVERAGE","STRAND"]
        if config.qc_nm:
            self.info_order.append("NM")

        self.default_genotype=config.genotype_none

        if config.mode=="combine":
            self.genotype_format=config.genotype_format+":ID"
            self.default_genotype+=tuple(["NULL"])
        else:
            self.genotype_format=config.genotype_format

        self.reference_handle=None

    def open_reference(self):
        if self.config.reference == None:
            return
        if not os.path.exists(self.config.reference+".fai") and not os.path.exists(self.config.reference+".gzi"):
            print(f"Info: Fasta index for {self.config.reference} not found. Generating with pysam.faidx (this may take a while)")
            pysam.faidx(self.config.reference)
        self.reference_handle=pysam.FastaFile(self.config.reference)

    def write_header(self,contigs_lengths):
        self.write_header_line("fileformat=VCFv4.2")
        self.write_header_line(f"source={self.config.version}_{self.config.build}")
        self.write_header_line('command="'+self.config.command+'"')
        self.write_header_line('fileDate="'+self.config.start_date+'"')
        for contig,contig_len in contigs_lengths:
            self.write_header_line(f"contig=<ID={contig},length={contig_len}>")

        self.write_header_line('ALT=<ID=INS,Description="Insertion">')
        self.write_header_line('ALT=<ID=DEL,Description="Deletion">')
        self.write_header_line('ALT=<ID=DUP,Description="Duplication">')
        self.write_header_line('ALT=<ID=INV,Description="Inversion">')
        self.write_header_line('ALT=<ID=BND,Description="Breakend; Translocation">')

        self.write_header_line('FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.write_header_line('FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        self.write_header_line('FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        self.write_header_line('FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')
        self.write_header_line('FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">')

        self.write_header_line('FILTER=<ID=PASS,Description="All filters passed">')
        self.write_header_line('FILTER=<ID=GT,Description="Genotype filter">')
        self.write_header_line('FILTER=<ID=CONS,Description="Low support for insertion consensus sequence">')
        self.write_header_line('INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">')
        self.write_header_line('INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">')
        self.write_header_line('INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">')
        self.write_header_line('INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">')
        self.write_header_line('INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">')
        self.write_header_line('INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">')
        self.write_header_line('INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only)">')
        self.write_header_line('INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">')
        self.write_header_line('INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">')
        self.write_header_line('INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">')
        self.write_header_line('INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">')
        self.write_header_line('INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">')
        self.write_header_line('INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">')
        self.write_header_line('INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">')
        self.write_header_line('INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">')
        #self.write_header_line('INFO=<ID=COVERAGE_REGIONAL,Number=A,Type=Float,Description="Coverage for local region">')
        #self.write_header_line('INFO=<ID=COVERAGE_GLOBAL,Number=A,Type=Float,Description="Coverage for global region">')
        #self.write_header_line('INFO=<ID=MINSUPPORT,Number=A,Type=Integer,Description="Determined minimum support threshold for this region">')
        #self.write_header_line('INFO=<ID=COVERAGE_UPSTREAM,Number=A,Type=Float,Description="Coverage upstream of structural variation start">')
        #self.write_header_line('INFO=<ID=COVERAGE_DOWNSTREAM,Number=A,Type=Float,Description="Coverage downstream of structural variation end">')
        self.write_header_line('INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">')
        self.write_header_line('INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">')
        self.write_header_line('INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">')
        self.write_header_line('INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">')
        self.write_header_line('INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">')
        #self.write_header_line('INFO=<ID=DBG_CLUSTER_SPAN,Number=A,Type=Float,Description="Debug info">')

        samples_header="\t".join(sample_id for internal_id,sample_id in self.config.sample_ids_vcf)
        self.write_raw(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_header}")

    def write_raw(self,text,endl="\n"):
        if self.config.vcf_output_bgz:
            self.handle.write(text.encode())
            self.handle.write(endl.encode())
        else:
            self.handle.write(text)
            self.handle.write(endl)

    def write_header_line(self,text):
        self.write_raw("##"+text)

    def write_call(self,call):
        #pysam coordinates are 0-based, VCF 1-based - therefore +1 is added to call.pos
        end=call.end + 1
        pos=call.pos + 1

        #Determine genotypes columns
        ac=0 #Allele count
        supvec=[]
        sample_genotypes=[]
        for internal_id, sample_id in self.config.sample_ids_vcf:
            if internal_id in call.genotypes and call.genotypes[internal_id]!=None:
                gt_curr=call.genotypes[internal_id]
                sample_genotypes.append(format_genotype(gt_curr))
                if gt_curr[0]!="." and gt_curr[4]>0: #Not non-genotype and has supporting reads
                    ac+=sum(call.genotypes[internal_id][:2])
                    supp="1"
                else:
                    supp="0"
            else:
                sample_genotypes.append(format_genotype(self.default_genotype))
                supp="0"
            supvec.append(supp)

        if len(self.config.sample_ids_vcf) > 1:
            call.set_info("AC",ac)
            call.set_info("SUPP_VEC","".join(supvec))
            if ac==0:
                call.filter="GT"

        #Output core SV attributes
        infos={"SVTYPE": call.svtype,
               "SVLEN": call.svlen,
               "END": end,
               "SUPPORT": call.support,
               "RNAMES": call.rnames if self.config.output_rnames else None,
               "COVERAGE": f"{call.coverage_upstream},{call.coverage_start},{call.coverage_center},{call.coverage_end},{call.coverage_downstream}",
               "STRAND": ("+" if call.fwd>0 else "") + ("-" if call.rev>0 else ""),
               "NM": call.nm}

        if call.svtype=="BND":
            infos["SVLEN"]=None
            infos["END"]=None

        infos_ordered=["PRECISE" if call.precise else "IMPRECISE"]
        infos_ordered.extend(format_info(k,infos[k]) for k in self.info_order if infos[k]!=None)
        info_str=";".join(infos_ordered)

        #Output call specific additional information
        for k in sorted(call.info):
            if call.info[k]==None:
                continue
            info_str+=";" + format_info(k,call.info[k])

        #if call.id==None:
        #    call.id=f"Sniffles2.{call.svtype}.{self.call_count+1:06}"

        #Resolve DEL sequence
        if not self.config.symbolic and call.svtype=="DEL" and self.reference_handle != None and abs(call.svlen) <= self.config.max_del_seq_len:
            try:
                call.ref=self.reference_handle.fetch(call.contig,call.pos,call.pos-call.svlen)
                call.alt="N"
            except KeyError:
                call.ref="N"
                call.alt=f"<{call.svtype}>"
            except ValueError:
                call.ref="N"
                call.alt=f"<{call.svtype}>"

        if self.config.symbolic:
            call.ref="N"
            call.alt=f"<{call.svtype}>"

        self.write_raw("\t".join(str(v) for v in [call.contig,pos,self.config.id_prefix+call.id,call.ref,call.alt,call.qual,call.filter,info_str,self.genotype_format]+sample_genotypes))
        self.call_count+=1

    def read_svs_iter(self):
        self.header_str=""
        line_index=0
        for line in self.handle:
            try:
                if isinstance(line,bytes):
                    line=line.decode("utf-8")
                line_index+=1
                line_strip=line.strip()
                if line_strip=="" or line_strip[0]=="#":
                    if line_strip[0]=="#":
                        self.header_str+=line_strip+"\n"
                    continue
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO=line.split("\t")[:8]
                info_dict={}
                for info_item in INFO.split(";"):
                    if "=" in info_item:
                        key,value=info_item.split("=")
                    else:
                        key,value=info_item,True
                    info_dict[key]=value
                call=sv.SVCall(contig=CHROM,
                               pos=int(POS)-1,
                               id=line_index,
                               ref=REF,
                               alt=ALT,
                               qual=QUAL,
                               filter=FILTER,
                               info=info_dict,
                               svtype=None,
                               svlen=None,
                               end=None,
                               rnames=None,
                               qc=True,
                               postprocess=None,
                               genotypes=None,
                               precise=None,
                               support=0,
                               fwd=0,
                               rev=0,
                               nm=-1)
                if len(call.alt)>len(call.ref):
                    call.svtype="INS"
                    call.svlen=len(call.alt)
                    call.end=call.pos
                else:
                    call.svtype="DEL"
                    call.svlen=-len(call.ref)
                    call.end=call.pos+call.svlen

                if "SVTYPE" in info_dict:
                    call.svtype=info_dict["SVTYPE"]
                    if call.svtype=="TRA":
                        call.svtype="BND"

                if "SVLEN" in info_dict:
                    call.svlen=int(info_dict["SVLEN"])
                if "END" in info_dict:
                    call.end=int(info_dict["END"])

                if call.svtype=="BND":
                    bnd_parts=call.alt.replace("]","[").split("[")
                    if len(bnd_parts)>2:
                        mate_contig,mate_ref_start=bnd_parts[1].split(":")
                        call.bnd_info=sv.SVCallBNDInfo(mate_contig=mate_contig, mate_ref_start=int(mate_ref_start), is_first=(call.alt[0]=="N"), is_reverse=("]" in call.alt))
                    else:
                        raise ValueError("BND ALT not formatted according to VCF 4.2 specifications")

                call.raw_vcf_line=line_strip
                call.raw_vcf_line_index=line_index
                yield call
            except Exception as e:
                util.fatal_error(f"Error parsing input VCF: Line {line_index}: {e}")

    def rewrite_genotype(self,svcall):
        parts_no_gt=svcall.raw_vcf_line.split("\t")[:8]
        gt_format=self.config.genotype_format
        if svcall.genotype_match_sv != None:
            if len(svcall.genotype_match_sv.genotypes)>0:
                gt=svcall.genotype_match_sv.genotypes[0]
            else:
                gt=svcall.genotypes[0]
        else:
            gt=svcall.genotypes[0] #0/0 or ./. depending on whether input SV had coverage in sample
        #parts=parts_no_gt + [gt_format,vcf.format_genotype(gt)]
        #gt_vcf=svcall.raw_vcf_line.split("\t")[9].split(":")[0]
        #parts= parts_no_gt + [gt_vcf] + [gt_format,vcf.format_genotype(gt)]
        parts=parts_no_gt + [gt_format,format_genotype(gt)]
        #parts[7]="NA"
        #parts[3]=f"REF_{len(parts[3])}"
        #parts[4]=f"ALT_{len(parts[4])}"
        self.write_raw("\t".join(parts))

    def rewrite_header_genotype(self,orig_header):
        header_lines=orig_header.split("\n")
        header_lines.insert(1,'##genotypeFileDate="'+self.config.start_date+'"')
        header_lines.insert(1,'##genotypeCommand="'+self.config.command+'"')
        header_lines.insert(1,f"##genotypeSource={self.config.version}_{self.config.build}")
        self.write_raw("\n".join(header_lines),endl="")

    def close(self):
        self.handle.close()
