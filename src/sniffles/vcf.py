#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     15.08.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
import logging
from collections import Counter

import pysam
import os

from sniffles import sv
from sniffles import util
from sniffles.config import SnifflesConfig


log = logging.getLogger(__name__)


def format_info(k, v):
    if isinstance(v, float):
        return f"{k}={v:.3f}"
    elif isinstance(v, list):
        return f"{k}={','.join(v)}"
    elif v is None:
        v = '.'

    return f"{k}={v}"


def unpack_phase(phase, svid="") -> tuple:
    try:
        hp_i, ps = phase
    except TypeError:
        if phase is None:
            hp_i, ps = ".", "."
        else:
            log.debug(f"Single not 'None'-valued phase: {phase}|{svid}")
            hp_i, ps = phase, phase
    ps = ps if ps is not None else "."
    return hp_i, ps


def format_genotype(gt, is_phased):
    """
    hp_i is the index of the haplotype in config.phase_identifiers:
    HP:1 => index 0 => phased genotype in the form of 1|0
    HP:2 => index 1 => phased genotype in the form of 0|1
    is_phased was added to modify the output as phasing adds the PS tag to FORMAT
    """
    gt_multi_sample_fields = 6
    if len(gt) == gt_multi_sample_fields:
        a, b, qual, dr, dv, phase = gt
        hp_i, ps = unpack_phase(phase)
        if hp_i is not None and (a, b) == (0, 1) and is_phased:
            gt_sep = "|"
            if hp_i == 0:
                a, b = b, a
        else:
            gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{ps}" if is_phased else f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}"
    else:
        a, b, qual, dr, dv, phase, svid = gt
        hp_i, ps = unpack_phase(phase, svid)
        if hp_i is not None and (a, b) == (0, 1) and is_phased:
            gt_sep = "|"
            if hp_i == 0:
                a, b = b, a
        else:
            gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{ps}:{svid}" if is_phased \
            else f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{svid}"


class VCF:
    def __init__(self, config: SnifflesConfig, handle):
        self.config = config
        self.handle = handle
        self.call_count = 0
        self.info_order = ["SVTYPE", "SVLEN", "END", "SUPPORT", "RNAMES", "COVERAGE", "STRAND"]
        if config.qc_nm_measure:
            self.info_order.append("NM")

        if config.dev_emit_sv_lengths:
            self.info_order.append("SVLENGTHS")

        self.default_genotype = config.genotype_none

        # Add phasing if needed
        self.genotype_format = config.genotype_format
        if config.phase:
            self.genotype_format += ":PS"
            # it has it added already as None
        if config.mode == "combine":
            self.genotype_format += ":ID"
            self.default_genotype += tuple(["NULL"])

        self.reference_handle = None
        self.header_str = ""

    def open_reference(self, generate_index: bool = True):
        if self.config.reference is None:
            return

        if not os.path.exists(self.config.reference + ".fai") and not os.path.exists(self.config.reference + ".gzi") and generate_index:
            log.warning(f"Fasta index for {self.config.reference} not found. Generating with pysam.faidx (this may take a while)")
            pysam.faidx(self.config.reference)

        try:
            self.reference_handle = pysam.FastaFile(self.config.reference)
        except:
            log.error(f'Unable to open reference file {self.config.reference}')

    def write_header(self, contigs_lengths):
        self.write_header_line("fileformat=VCFv4.2")
        self.write_header_line(f"source={self.config.version}_{self.config.build}")
        self.write_header_line('command="' + self.config.command + '"')
        self.write_header_line('fileDate="' + self.config.start_date + '"')
        for contig, contig_len in contigs_lengths:
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
        self.write_header_line('FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase-block, zero if none or not phased">')
        self.write_header_line('FORMAT=<ID=ID,Number=1,Type=String,Description="Individual sample SV ID for multi-sample output">')

        self.write_header_line('FILTER=<ID=PASS,Description="All filters passed">')
        self.write_header_line('FILTER=<ID=GT,Description="Genotype filter">')
        self.write_header_line('FILTER=<ID=SUPPORT_MIN,Description="Minimum read support filter">')
        self.write_header_line('FILTER=<ID=STDEV_POS,Description="SV Breakpoint standard deviation filter">')
        self.write_header_line('FILTER=<ID=STDEV_LEN,Description="SV length standard deviation filter">')
        self.write_header_line('FILTER=<ID=COV_MIN,Description="Minimum coverage filter">')
        self.write_header_line('FILTER=<ID=COV_MIN_GT,Description="Minimum coverage filter (missing genotype)">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_DEL,Description="Coverage change filter for DEL">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_DUP,Description="Coverage change filter for DUP">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_INS,Description="Coverage change filter for INS">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage fractional change filter: upstream-start">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage fractional change filter: start-center">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage fractional change filter: center-end">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage fractional change filter: end-downstream">')
        self.write_header_line('FILTER=<ID=COV_VAR,Description="Coverage variance exceeded">')
        self.write_header_line('FILTER=<ID=MOSAIC_VAF,Description="Mosaic variant allele fraction filter">')
        self.write_header_line('FILTER=<ID=NOT_MOSAIC_VAF,Description="Variant allele fraction filter for non-mosaic">')
        self.write_header_line('FILTER=<ID=ALN_NM,Description="Length adjusted mismatch filter">')
        self.write_header_line('FILTER=<ID=STRAND_BND,Description="Strand support filter for BNDs">')
        self.write_header_line('FILTER=<ID=STRAND,Description="Strand support filter for germline SVs">')
        self.write_header_line('FILTER=<ID=STRAND_MOSAIC,Description="Strand support filter for mosaic SVs">')
        self.write_header_line('FILTER=<ID=SVLEN_MIN,Description="SV length filter">')
        self.write_header_line('FILTER=<ID=SVLEN_MIN_MOSAIC,Description="SV length filter for mosaic SVs">')
        self.write_header_line('INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Structural variation with precise breakpoints">')
        self.write_header_line('INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Structural variation with imprecise breakpoints">')
        self.write_header_line('INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Structural variation classified as putative mosaic">')
        self.write_header_line('INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variation">')
        if self.config.dev_emit_sv_lengths:
            self.write_header_line('INFO=<ID=SVLENGTHS,Number=.,Type=Integer,Description="Lengths of structural variation (all)">')
        self.write_header_line('INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variation">')
        self.write_header_line('INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">')
        self.write_header_line('INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the structural variation">')
        self.write_header_line('INFO=<ID=SUPPORT_INLINE,Number=1,Type=Integer,Description="Number of reads supporting an INS/DEL SV (non-split events only)">')
        self.write_header_line('INFO=<ID=SUPPORT_SA,Number=1,Type=Integer,Description="Number of reads supporting a DEL SV through supplementary alignments (split events)">')
        self.write_header_line('INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Number of soft-clipped reads putatively supporting the long insertion SV">')
        self.write_header_line('INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">')
        self.write_header_line('INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Standard deviation of structural variation start position">')
        self.write_header_line('INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Standard deviation of structural variation length">')
        self.write_header_line('INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Coverages near upstream, start, center, end, downstream of structural variation">')
        self.write_header_line('INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of supporting reads for structural variant">')
        self.write_header_line('INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count, summed up over all samples">')
        self.write_header_line('INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="List of read support for all samples">')
        self.write_header_line('INFO=<ID=CONSENSUS_SUPPORT,Number=1,Type=Integer,Description="Number of reads that support the generated insertion (INS) consensus sequence">')
        self.write_header_line('INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads (if enabled with --output-rnames)">')
        self.write_header_line('INFO=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction">')
        self.write_header_line('INFO=<ID=COVERAGE_VAR,Number=1,Type=Float,Description="Variance of coverage across large events">')
        self.write_header_line('INFO=<ID=NM,Number=.,Type=Float,Description="Mean number of query alignment length adjusted mismatches of supporting reads">')
        self.write_header_line('INFO=<ID=PHASE,Number=.,Type=String,Description="Phasing information derived from supporting reads, represented as list of: HAPLOTYPE,PHASESET,HAPLOTYPE_SUPPORT,PHASESET_SUPPORT,HAPLOTYPE_FILTER,PHASESET_FILTER">')

        if self.config.combine_population:
            self.write_header_line('INFO=<ID=POPULATION_AF,Number=1,Type=Float,Description="Population Allele Frequency">')
            self.write_header_line('INFO=<ID=POPULATION_SIZE,Number=1,Type=Integer,Description="Size of genotyped population for this variant">')

        samples_header = "\t".join(sample_id for _, sample_id in self.config.sample_ids_vcf)
        self.write_raw(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_header}")

    def write_raw(self, text, endl="\n"):
        if self.config.vcf_output_bgz:
            self.handle.write(text.encode())
            self.handle.write(endl.encode())
        else:
            self.handle.write(text)
            self.handle.write(endl)

    def write_header_line(self, text):
        self.write_raw("##" + text)

    def write_call(self, call: sv.SVCall) -> int:
        # pysam coordinates are 0-based, VCF 1-based
        # but VCF also requires the index of the base before the SV to be reported,
        # so we are fine without offsetting
        end = call.end
        pos = call.pos if call.pos > 0 else 1

        # Determine genotypes columns
        ac = 0  # Allele count
        supvec = []
        sample_genotypes = []
        for internal_id, _ in self.config.sample_ids_vcf:
            if internal_id in call.genotypes and call.genotypes[internal_id] is not None:
                gt_curr = call.genotypes[internal_id]
                sample_genotypes.append(format_genotype(gt_curr, self.config.phase))
                if gt_curr[0] != "." and gt_curr[4] > 0:  # Not non-genotype and has supporting reads
                    ac += sum(call.genotypes[internal_id][:2])
                    supp = "1"
                else:
                    supp = "0"
            else:
                sample_genotypes.append(format_genotype(self.default_genotype, self.config.phase))
                supp = "0"
            supvec.append(supp)

        if len(self.config.sample_ids_vcf) > 1:
            call.set_info("AC", ac)
            call.set_info("SUPP_VEC", svec := "".join(supvec))

            if int(svec) == 0:
                log.debug(f'Dropped {call} due to all zero support vector.')
                return 0

            if ac == 0:
                call.filter = "GT"

        # Check if svlen == len(alt) in INS, NOT len(alt)-1 as the ref base is added after this step
        if "INS" == call.svtype:
            if call.svlen != len(call.alt) and not self.config.symbolic and "<INS>" != call.alt:
                log.debug(f"Updating SVLEN for INS to match sequence length: {call.svlen} v {len(call.alt)}")
                call.svlen = len(call.alt)
            if call.svlen < self.config.minsvlen:
                return 0

        # Output core SV attributes
        infos = {
            "SVTYPE": call.svtype,
            "SVLEN": call.svlen,
            "SVLENGTHS": ",".join(map(str, call.svlens)) if call.svlens else None,
            "END": end,
            "SUPPORT": call.support,
            "RNAMES": call.rnames if self.config.output_rnames else None,
            "COVERAGE": f"{call.coverage_upstream},{call.coverage_start},{call.coverage_center},{call.coverage_end},"
                        f"{call.coverage_downstream}",
            "STRAND": ("+" if call.fwd > 0 else "") + ("-" if call.rev > 0 else ""),
            "NM": call.nm
        }

        if call.svtype == "BND":
            infos["SVLEN"] = None
            infos["SVLENGTHS"] = None
            infos["END"] = None

        if call.svtype == "DEL":
            # END is POS + length of REF allele - 1
            infos["END"] = end - 1

        infos_ordered = ["PRECISE" if call.precise else "IMPRECISE"]
        af = call.get_info("VAF")
        af = af if af is not None else 0
        sv_is_mosaic = af <= self.config.mosaic_af_max
        if sv_is_mosaic and self.config.mosaic:
            infos_ordered.append("MOSAIC")
        infos_ordered.extend(format_info(k, infos[k]) for k in self.info_order if infos[k] is not None)
        info_str = ";".join(infos_ordered)

        # Output call specific additional information
        for k in sorted(call.info):
            if call.info[k] is None:
                continue
            info_str += ";" + format_info(k, call.info[k])

        # if call.id==None:
        #    call.id=f"Sniffles2.{call.svtype}.{self.call_count+1:06}"

        # Resolve DEL sequence
        if (not self.config.symbolic and call.svtype == "DEL" and self.reference_handle is not None
                and abs(call.svlen) <= self.config.max_del_seq_len):
            try:
                # VCF requires inclusion of the last reference base before the SV
                call.ref = self.reference_handle.fetch(call.contig, call.pos - 1, call.pos - call.svlen)
                call.alt = call.ref[0]
            except KeyError:
                call.ref = "N"
                call.alt = f"<{call.svtype}>"
            except ValueError:
                call.ref = "N"
                call.alt = f"<{call.svtype}>"
            else:
                if 'N' in call.ref and (pct_n := Counter(call.ref)['N'] / len(call.ref)) > self.config.max_unknown_pct:
                    # don't emit calls with too many N bases
                    log.debug(f'Not emitting {call.id} (length {call.svlen}) due to {pct_n*100:.2f}% N bases in reference.')
                    return 0

        if self.config.symbolic:
            call.ref = "N"
            call.alt = f"<{call.svtype}>"
        else:
            if self.reference_handle is not None and call.ref == 'N':
                # Fetch the base before the SV
                try:
                    call.ref = self.reference_handle.fetch(call.contig, start := max(0, call.pos - 1), start + 1)
                except (KeyError, ValueError):
                    ...
                else:
                    if call.svtype == "INS" and call.alt != '<INS>':
                        call.alt = call.ref + call.alt
                    elif call.svtype == 'BND' and call.alt != '<BND>':
                        call.alt = (call.ref + call.alt[1:]) if call.alt.startswith('N') else call.alt[:-1] + call.ref

        call.qual = max(0, min(60, call.qual)) if call.qual is not None else None

        self.write_raw("\t".join(str(v) for v in [call.contig, pos, self.config.id_prefix + call.id, call.ref,
                                                  call.alt, call.qual if call.qual is not None else '.', call.filter,
                                                  info_str, self.genotype_format] + sample_genotypes))
        self.call_count += 1
        return 1

    def read_svs_iter(self):
        self.header_str = ""
        line_index = 0
        for line in self.handle:
            try:
                if isinstance(line, bytes):
                    line = line.decode("utf-8")
                line_index += 1
                line_strip = line.strip()
                if line_strip == "" or line_strip[0] == "#":
                    if line_strip[0] == "#":
                        self.header_str += line_strip + "\n"
                    continue
                CHROM, POS, _, REF, ALT, QUAL, FILTER, INFO = line.split("\t")[:8]
                info_dict = {}
                for info_item in INFO.split(";"):
                    if "=" in info_item:
                        key, value = info_item.split("=")
                    else:
                        key, value = info_item, True
                    info_dict[key] = value
                call = sv.SVCall(contig=CHROM,
                                 pos=int(POS) - 1,
                                 id=line_index,
                                 ref=REF,
                                 alt=ALT,
                                 qual=int(QUAL) if QUAL != '.' else None,
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
                if len(call.alt) > len(call.ref):
                    call.svtype = "INS"
                    call.svlen = len(call.alt)
                    call.end = call.pos
                else:
                    call.svtype = "DEL"
                    call.svlen = -len(call.ref)
                    call.end = call.pos + call.svlen

                if "SVTYPE" in info_dict:
                    call.svtype = info_dict["SVTYPE"]
                    if call.svtype == "TRA":
                        call.svtype = "BND"

                if "SVLEN" in info_dict:
                    call.svlen = int(info_dict["SVLEN"])
                if "SVLENGTHS`" in info_dict:
                    call.svlens = info_dict["SVLENGTHS"]

                if "END" in info_dict:
                    call.end = int(info_dict["END"])

                if call.svtype == "BND":
                    bnd_parts = call.alt.replace("]", "[").split("[")
                    if len(bnd_parts) > 2:
                        mate_contig, mate_ref_start = bnd_parts[1].split(":")
                        call.bnd_info = sv.SVCallBNDInfo(mate_contig=mate_contig, mate_ref_start=int(mate_ref_start),
                                                         is_first=(call.alt[0] == "N"), is_reverse=("]" in call.alt))
                    else:
                        raise ValueError("BND ALT not formatted according to VCF 4.2 specifications")

                call.raw_vcf_line = line_strip
                call.raw_vcf_line_index = line_index
                yield call
            except Exception as e:
                util.fatal_error(f"Error parsing input VCF: Line {line_index}: {e}")

    def rewrite_genotype(self, svcall):
        parts_no_gt = svcall.raw_vcf_line.split("\t")[:8]
        gt_format = self.config.genotype_format
        if svcall.genotype_match_sv != None:
            if len(svcall.genotype_match_sv.genotypes) > 0:
                gt = svcall.genotype_match_sv.genotypes[0]
            else:
                gt = svcall.genotypes[0]
        else:
            gt = svcall.genotypes[0]  # 0/0 or ./. depending on whether input SV had coverage in sample
        # parts=parts_no_gt + [gt_format,vcf.format_genotype(gt)]
        # gt_vcf=svcall.raw_vcf_line.split("\t")[9].split(":")[0]
        # parts= parts_no_gt + [gt_vcf] + [gt_format,vcf.format_genotype(gt)]
        parts = parts_no_gt + [gt_format, format_genotype(gt, self.config.phase)]
        # parts[7]="NA"
        # parts[3]=f"REF_{len(parts[3])}"
        # parts[4]=f"ALT_{len(parts[4])}"
        self.write_raw("\t".join(parts))

    def rewrite_header_genotype(self,orig_header):
        header_lines=orig_header.split("\n")
        header_lines.insert(1,'##genotypeFileDate="'+self.config.start_date+'"')
        header_lines.insert(1,'##genotypeCommand="'+self.config.command+'"')
        header_lines.insert(1,f"##genotypeSource={self.config.version}_{self.config.build}")

        # Fix missing genotype headers errors when reading input vcf in genotype mode.
        # This error will happen when input vcf does not have GT,GQ,DR,DV headers.
        # Some tools such as turvari will throw errors when processing the output vcf.
        has_gt_headers = {
            "GT": False,
            "GQ": False,
            "DR": False,
            "DV": False,
        }
        for header_line in header_lines:
            for gt in has_gt_headers.keys():
                if "##FORMAT=<ID="+gt+"," in header_line:
                    has_gt_headers[gt] = True
        
        if not has_gt_headers["GT"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        if not has_gt_headers["GQ"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        if not has_gt_headers["DR"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        if not has_gt_headers["DV"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')

        self.write_raw("\n".join(header_lines), endl="")

    def close(self):
        self.handle.close()
