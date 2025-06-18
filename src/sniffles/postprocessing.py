#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     27.08.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
from enum import IntEnum
import math
import logging

from sniffles import util
from sniffles import consensus
from sniffles.config import SnifflesConfig
from sniffles.sv import SVCall

log = logging.getLogger('sniffles.postprocessing')

class CoverageMode(IntEnum):
    """A cheap(er) class to represent the coverage mode used when setting SVCall attributes."""
    START = 1
    CENTER = 2
    END = 3
    UPSTREAM = 4
    DOWNSTREAM = 5

    @property
    def to_attr(self):
        return {
            CoverageMode.START: "coverage_start",
            CoverageMode.CENTER: "coverage_center",
            CoverageMode.END: "coverage_end",
            CoverageMode.UPSTREAM: "coverage_upstream",
            CoverageMode.DOWNSTREAM: "coverage_downstream",
        }[self]


def annotate_sv(svcall: SVCall, config):
    if config.phase:
        phase = phase_sv(svcall, config)
    else:
        phase = (None, None)

    genotype_sv(svcall, config, phase)

    if svcall.svtype == "INS" and not config.symbolic:
        merged_leads = [lead for lead in svcall.postprocess.cluster.leads if lead.seq is not None]

        if len(merged_leads):
            best_lead = merged_leads[0]
            best_index = 0
            best_diff = abs(len(best_lead.seq) - svcall.svlen) + abs(best_lead.ref_start - svcall.pos) * 1.5
            for i, ld in enumerate(merged_leads):
                if i == 0:
                    continue
                curr_diff = abs(len(ld.seq) - svcall.svlen) + abs(ld.ref_start - svcall.pos) * 1.5
                if curr_diff < best_diff:
                    best_lead = ld
                    best_index = i
                    best_diff = curr_diff

            merged_leads.pop(best_index)
            # merged_leads_new=list()

            # for lead in merged_leads:
            #    curr_lend_iff=abs(len(ld.seq) - len(best_lead.seq)) + abs(ld.ref_start - best_lead.ref_start)
            #    #if curr_len_diff < 25:
            #    merged_leads_new.append(lead)
            # merged_leads=merged_leads_new

            if len(merged_leads) >= config.consensus_min_reads and not config.no_consensus:
                kmer_len = config.consensus_kmer_len
                skip = config.consensus_kmer_skip_base + int(len(best_lead.seq)*config.consensus_kmer_skip_seqlen_mult)
                skip_repetitive = skip

                svcall.alt = consensus.novel_from_reads(best_lead, merged_leads, klen=kmer_len, skip=skip,
                                                        skip_repetitive=skip_repetitive)
            else:
                svcall.alt = best_lead.seq


def add_request(svcall_i: int, field: CoverageMode, pos, requests_for_coverage: dict[int, set[tuple[int, int]]], config):
    """
    Add a request for one of the five coverage fields to the given SVCall
    """
    bin_index = int(pos / config.coverage_binsize) * config.coverage_binsize
    if bin_index not in requests_for_coverage:
        requests_for_coverage[bin_index] = set()
    requests_for_coverage[bin_index].add( (svcall_i, field) )


def coverage(calls, lead_provider, config):
    requests_for_coverage_attrs, requests_for_coverage_sample_starts, requests_for_coverage_sample_ends = coverage_build_requests(calls, config)
    return coverage_fulfill(calls, requests_for_coverage_attrs, requests_for_coverage_sample_starts, requests_for_coverage_sample_ends, lead_provider, config)


def coverage_build_requests(calls, config: SnifflesConfig):
    """
    Updated: changed coverage calculation to be inside the SV for INV/DEL/DUP, outside (like it is now) for INS/BND?
    BIN: ---
    INS/BND    |-------------V---------------|
              U---       S--- E---            D---
                            C---
    DUP/DEL    |------------===============------------|
    INV       U---          S---       E---            D---
                                  C---
    """
    requests_for_coverage_attrs = {}
    requests_for_coverage_sample_starts = {}
    requests_for_coverage_sample_ends = {}

    for sv_i, svcall in enumerate(calls):
        start = svcall.pos
        if svcall.svtype == "INS":
            end = start + 1
        else:
            end = svcall.pos + abs(svcall.svlen)

        if svcall.svtype in ("DEL", 'INV', 'DUP') and abs(svcall.svlen) >= config.long_del_length:
            # Sampling more intervals for large deletions
            clean_start = start
            if start < 0:
                # Excusing a negative start, we'll find the first positive position we'd encounter striding our large coverage interval and start there instead
                clean_start = start + ((-start // config.large_coverage_sample_interval) + 1) * config.large_coverage_sample_interval
                log.warning(f"[sv] Encountered SV={svcall.id} with negative pos={start}. Setting to first viable positive pos={clean_start} for coverage postprocessing.")
            first_bin = (clean_start // config.coverage_binsize) * config.coverage_binsize
            last_bin_pos = clean_start + ((end - clean_start - 1) // config.large_coverage_sample_interval) * config.large_coverage_sample_interval
            last_bin = (last_bin_pos // config.coverage_binsize) * config.coverage_binsize

            if first_bin not in requests_for_coverage_sample_starts:
                requests_for_coverage_sample_starts[first_bin] = set()
            if last_bin not in requests_for_coverage_sample_ends:
                requests_for_coverage_sample_ends[last_bin] = set()
            requests_for_coverage_sample_starts[first_bin].add(sv_i)
            requests_for_coverage_sample_ends[last_bin].add(sv_i)

        if svcall.svtype in ["INS", "BND"]:
            add_request(sv_i, CoverageMode.START, start - config.coverage_binsize, requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.CENTER, int((start + end - config.coverage_binsize) / 2),
                        requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.END, end + config.coverage_binsize, requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.UPSTREAM, start - config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.DOWNSTREAM, end + config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage_attrs, config)
        else:
            add_request(sv_i, CoverageMode.START, start, requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.CENTER, int((start + end - config.coverage_binsize) / 2),
                        requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.END, end - config.coverage_binsize, requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.UPSTREAM, start - config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage_attrs, config)
            add_request(sv_i, CoverageMode.DOWNSTREAM, end + config.coverage_binsize * config.coverage_updown_bins,
                        requests_for_coverage_attrs, config)

    return requests_for_coverage_attrs, requests_for_coverage_sample_starts, requests_for_coverage_sample_ends


def coverage_fulfill(calls, requests_for_coverage_attrs, requests_for_coverage_sample_starts, requests_for_coverage_sample_ends, lead_provider, config: SnifflesConfig):
    if len(requests_for_coverage_attrs) == 0 and len(requests_for_coverage_sample_starts) == 0:
        return -1, -1

    start_bin = lead_provider.covrtab_min_bin
    end_bin = int(lead_provider.end / config.coverage_binsize) * config.coverage_binsize
    coverage_fwd = 0
    coverage_rev = 0
    coverage_fwd_total = 0
    coverage_rev_total = 0
    n = 0

    # We can efficiently store and query for large coverage sampling intervals (LCSI) by
    #  maintaining a set for each possible bin modulus which stores the indices of SVCalls
    #  that are open to be sampled at all bins with the same modulus.
    # The diagram below illustrates a series of bins and how this approach would work.
    # Consider an LCSI every 4 bins (actual positions ignored for simplicity).
    # Here SV_i starts at a position between bin 1 and 2 and ends at a position
    #  between bin 8 and 9. As the coverage postprocessing forces sampling to start and
    #  end on a coverage bin, sampling of SV_i will start at bin 1 and end at bin 8.
    #  Bin 1 has a modulus of 1 with respect to the LCSI of 4, as does 5, 9 and so on.
    # For all the possible moduli (in this case, [0, 4]) we can maintain a set of
    # SVCall indices that are open to be sampled at the next bin with the same modulus.
    #
    #   id   0    1    2    3    4    5    6    7    8    9   10   11    12 ..
    # id%4   0    1    2    3    0    1    2    3    4    1    2    3    0 ...
    # Bins   o----o----o----o----o----o----o----o----o----o----o----o----o ...
    # LCSI%0 X-------------------X-------------------X-------------------X ...
    # SV_i        ===S-------------------------------X==E
    #             ^ First Bin has modulus 1          ^ Last Bin
    #               we'll sample all modulus 1 bins before last bin
    # LCSI%1       X-------------------X-------------------X-------------->...
    # SV_i LCSI    X                   X
    bin_moduli = {}

    for bin_pos in range(start_bin, end_bin + config.coverage_binsize, config.coverage_binsize):
        bin_modulus = bin_pos % config.large_coverage_sample_interval
        n += 1

        coverage_fwd += lead_provider.covrtab_fwd.get(bin_pos, 0)  # TODO -- should this be reset
        coverage_rev += lead_provider.covrtab_rev.get(bin_pos, 0)  # TODO -- should this be reset
        coverage_total_curr = coverage_fwd + coverage_rev

        # handle attr style coverage requests for this bin position
        for sv_i, mode in requests_for_coverage_attrs.get(bin_pos, []):
            setattr(calls[sv_i], mode.to_attr, coverage_total_curr)

        # open
        # if this bin position opens SVs, add each SVCall index to the set of SVs to be
        # updated for any bin with the same modulo
        if bin_pos in requests_for_coverage_sample_starts:
            if bin_modulus not in bin_moduli:
                bin_moduli[bin_modulus] = set()
            bin_moduli[bin_modulus].update(requests_for_coverage_sample_starts[bin_pos])

        # sample welfordly
        for sv_i in bin_moduli.get(bin_modulus, []):
            calls[sv_i].forward_difference_sampler.push(coverage_total_curr)

        # close
        # this is the last bin, do not process any more overlaps on the next iteration
        if bin_pos in requests_for_coverage_sample_ends:
            bin_moduli[bin_modulus].difference_update(requests_for_coverage_sample_ends[bin_pos])

        coverage_fwd_total += coverage_fwd
        coverage_rev_total += coverage_rev

    average_coverage_fwd = coverage_fwd_total / float(n) if n > 0 else 0
    average_coverage_rev = coverage_rev_total / float(n) if n > 0 else 0
    return average_coverage_fwd, average_coverage_rev


def qc_sv_support(svcall, coverage_global, config):
    if config.mosaic:
        return True
    if config.minsupport == "auto":
        if not qc_support_auto(svcall, coverage_global, config):
            svcall.filter = "SUPPORT_MIN"
            return False
    else:
        if not qc_support_const(svcall, config):
            svcall.filter = "SUPPORT_MIN"
            return False
    return True


def rescale_support(svcall, config: SnifflesConfig):
    if svcall.svtype != "INS" or svcall.svlen < config.long_ins_length:
        return svcall.support
    else:
        base = svcall.support
        scale_factor = config.long_ins_rescale_mult * (float(svcall.svlen) / config.long_ins_length)
        return round(base * (config.long_ins_rescale_base + scale_factor))


def qc_support_auto(svcall, coverage_global, config):
    support = rescale_support(svcall, config)
    # if svcall.svtype=="INS":
    #    coverage_list=[svcall.coverage_center]
    # else:
    coverage_list = [each_coverage for each_coverage in [svcall.coverage_upstream, svcall.coverage_downstream] if
                     each_coverage != 0]
    if len(coverage_list) == 0:
        coverage_list = [each_coverage for each_coverage in [svcall.coverage_start, svcall.coverage_center,
                                                             svcall.coverage_end] if each_coverage != 0]
    if len(coverage_list) == 0:
        coverage_regional = coverage_global
    else:
        coverage_regional = round(sum(coverage_list) / len(coverage_list))
        if coverage_regional == 0:
            coverage_regional = coverage_global
    coverage_global_weight = (1.0 - config.minsupport_auto_regional_coverage_weight)
    coverage = (coverage_regional * config.minsupport_auto_regional_coverage_weight +
                coverage_global * coverage_global_weight)
    min_support = round(config.minsupport_auto_base + config.minsupport_auto_mult * coverage)
    return support >= min_support
    # return True


def qc_support_const(svcall, config):
    # svcall.set_info("MINSUPPORT",config.minsupport)
    return svcall.support >= config.minsupport


def qc_sv(svcall: SVCall, config: SnifflesConfig):
    af = svcall.get_info("VAF")
    af = af if af is not None else 0
    sv_is_mosaic = af <= config.mosaic_af_max

    if config.qc_stdev:
        stdev_pos = svcall.get_info("STDEV_POS")
        if stdev_pos > config.qc_stdev_abs_max:
            svcall.filter = "STDEV_POS"
            return False
        if svcall.svtype != "BND" and stdev_pos / abs(svcall.svlen) > 2.0:
            svcall.filter = "STDEV_POS"
            return False

        stdev_len = svcall.get_info("STDEV_LEN")
        if stdev_len is not None:
            if svcall.svtype != "BND" and stdev_len / abs(svcall.svlen) > 1.0:
                svcall.filter = "STDEV_LEN"
                return False
            if stdev_len > config.qc_stdev_abs_max:
                svcall.filter = "STDEV_LEN"
                return False

    if config.mosaic and sv_is_mosaic:
        if svcall.support < config.mosaic_min_reads:
            svcall.filter = "SUPPORT_MIN"
            return False

    if abs(svcall.svlen) < config.minsvlen:
        svcall.filter = "SVLEN_MIN"
        return False

    if svcall.svtype == "BND":
        if config.qc_bnd_filter_strand and len(set(lead.strand for lead in svcall.postprocess.cluster.leads)) < 2:
            svcall.filter = "STRAND_BND"
            return False
    elif not (config.mosaic and sv_is_mosaic) and config.qc_strand:
        is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
        if not is_long_ins and len(set(leads.strand for leads in svcall.postprocess.cluster.leads)) < 2:
            svcall.filter = "STRAND"
            return False
    elif (config.mosaic and sv_is_mosaic) and config.mosaic_qc_strand:
        is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
        if (not is_long_ins and len(set(leads.strand for leads in svcall.postprocess.cluster.leads)) < 2
                and svcall.support >= config.mosaic_use_strand_thresholds):
            svcall.filter = "STRAND_MOSAIC"
            return False

    if config.mosaic and sv_is_mosaic:
        if (svcall.svtype == "INV" or svcall.svtype == "DUP") and svcall.svlen < config.mosaic_qc_invdup_min_length:
            svcall.filter = "SVLEN_MIN_MOSAIC"
            return False

    if svcall.coverage_center < config.qc_coverage and svcall.svtype not in ["DEL", "INS"]:
        if (svcall.svtype == "INV" and svcall.svlen) > config.long_inv_length and not (config.mosaic and sv_is_mosaic):
            pass
        else:
            svcall.filter = "COV_MIN"
            return False

    upstream_downstream_max_coverage_diff = 0.7  # 70%
    upstream_downstream_diff = 0.5  # add to config
    if (svcall.svtype == "DEL" and config.long_del_length != -1 and abs(svcall.svlen) >= config.long_del_length and
            not config.mosaic):
        scaled_long_del_coverage = config.long_del_coverage/2.0   # 0.66/2 = 0.33
        if svcall.coverage_center > (svcall.coverage_upstream + svcall.coverage_downstream) * scaled_long_del_coverage:
            # check if slopped coverage, that often happens over large spans
            if svcall.coverage_upstream > svcall.coverage_center > svcall.coverage_downstream:
                if svcall.coverage_downstream/svcall.coverage_upstream < upstream_downstream_max_coverage_diff:
                    svcall.filter = "COV_CHANGE_DEL"
                    return False
            elif svcall.coverage_upstream < svcall.coverage_center < svcall.coverage_downstream:
                if svcall.coverage_upstream/svcall.coverage_downstream < upstream_downstream_max_coverage_diff:
                    svcall.filter = "COV_CHANGE_DEL"
                    return False
            else:
                pass
        if svcall.coverage_upstream > svcall.coverage_downstream:
            if (upstream_downstream_diff > svcall.coverage_downstream/svcall.coverage_upstream or
                    svcall.coverage_center > svcall.coverage_downstream):
                svcall.filter = "COV_CHANGE_DEL"
                return False
        elif svcall.coverage_upstream < svcall.coverage_downstream:
            if (upstream_downstream_diff > svcall.coverage_upstream/svcall.coverage_downstream or
                    svcall.coverage_upstream < svcall.coverage_center):
                svcall.filter = "COV_CHANGE_DEL"
                return False
        else:
            pass
    elif (svcall.svtype == "DUP" and config.long_dup_length != -1 and abs(svcall.svlen) >= config.long_dup_length and
          not config.mosaic):
        scaled_long_dup_coverage = config.long_dup_coverage / 2.0
        if svcall.coverage_center < (svcall.coverage_upstream + svcall.coverage_downstream) * scaled_long_dup_coverage:
            # check if slopped coverage, that often happens over large spans
            if svcall.coverage_upstream > svcall.coverage_center > svcall.coverage_downstream:
                if svcall.coverage_downstream/svcall.coverage_upstream < upstream_downstream_max_coverage_diff:
                    svcall.filter = "COV_CHANGE_DUP"
                    return False
            elif svcall.coverage_upstream < svcall.coverage_center < svcall.coverage_downstream:
                if svcall.coverage_upstream/svcall.coverage_downstream < upstream_downstream_max_coverage_diff:
                    svcall.filter = "COV_CHANGE_DUP"
                    return False
            else:
                pass
            if svcall.coverage_upstream > svcall.coverage_downstream:
                if (upstream_downstream_diff > svcall.coverage_downstream / svcall.coverage_upstream or
                        svcall.coverage_center < svcall.coverage_downstream):
                    svcall.filter = "COV_CHANGE_DUP"
                    return False
            elif svcall.coverage_upstream < svcall.coverage_downstream:
                if (upstream_downstream_diff > svcall.coverage_upstream / svcall.coverage_downstream or
                        svcall.coverage_upstream > svcall.coverage_center):
                    svcall.filter = "COV_CHANGE_DUP"
                    return False
            else:
                pass
    elif svcall.svtype == "INS" and (svcall.coverage_upstream < config.qc_coverage or
                                     svcall.coverage_downstream < config.qc_coverage):
        svcall.filter = "COV_CHANGE_INS"
        return False

    qc, val = svcall.qc_coverage_samples()
    svcall.set_info('COVERAGE_VAR', val)
    if not qc:
        svcall.filter = f'COV_VAR'
        return False

    qc_coverage_max_change_frac = config.qc_coverage_max_change_frac
    if config.mosaic and sv_is_mosaic:
        qc_coverage_max_change_frac = config.mosaic_qc_coverage_max_change_frac
    if qc_coverage_max_change_frac != -1.0:
        if svcall.coverage_upstream != 0:
            u = float(svcall.coverage_upstream)
        else:
            u = 1.0

        if svcall.coverage_start != 0:
            s = float(svcall.coverage_start)
        else:
            s = 1.0

        if svcall.coverage_center != 0:
            c = float(svcall.coverage_center)
        else:
            c = 1.0

        if svcall.coverage_end != 0:
            e = float(svcall.coverage_end)
        else:
            e = 1.0

        if svcall.coverage_downstream != 0:
            d = float(svcall.coverage_downstream)
        else:
            d = 1.0

        if abs(u - s) / max(u, s) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_US"
            return False

        if abs(s - c) / max(s, c) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_SC"
            return False

        if abs(c - e) / max(c, e) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_CE"
            return False

        if abs(e - d) / max(e, d) > qc_coverage_max_change_frac:
            svcall.filter = "COV_CHANGE_FRAC_ED"
            return False

    return True


def qc_sv_post_annotate(svcall: SVCall, config: SnifflesConfig):
    af = svcall.get_info("VAF")
    af = af if af is not None else 0
    sv_is_mosaic = af <= config.mosaic_af_max

    if (svcall.coverage_center < config.qc_coverage and
            (len(svcall.genotypes) == 0 or (svcall.genotypes[0][0] != "." and
                                            svcall.genotypes[0][0] + svcall.genotypes[0][1] < 2))):
        svcall.filter = "COV_MIN_GT"
        return False

    qc_nm = config.qc_nm
    qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if config.mosaic and sv_is_mosaic:
        qc_nm = config.mosaic_qc_nm
        qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if qc_nm and svcall.nm > qc_nm_threshold and (len(svcall.genotypes) == 0 or svcall.genotypes[0][1] == 0):
        svcall.filter = "ALN_NM"
        return False

    if config.mosaic:
        if sv_is_mosaic and (af < config.mosaic_af_min or af > config.mosaic_af_max):
            svcall.filter = "MOSAIC_VAF"
            return False
        elif not sv_is_mosaic and not config.mosaic_include_germline:
            svcall.filter = "NOT_MOSAIC_VAF"
            return False

    return True


def binomial_coef(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


def genotype_sv(svcall: SVCall, config, phase: tuple | None = None):
    from sniffles.genotyping import GENOTYPER_BY_TYPE, Genotyper

    GENOTYPER_BY_TYPE.get(svcall.svtype, Genotyper)(svcall, config, phase).calculate()


def phase_sv(svcall, config):
    reads_phases = {lead.read_id[0]: (lead.read_id[1], lead.read_id[2]) for lead in svcall.postprocess.cluster.leads}
    hp_list = util.most_common(hp for hp, ps in reads_phases.values())
    ps_list = util.most_common(ps for hp, ps in reads_phases.values())

    hp_support, hp = hp_list[0]
    ps_support, ps = ps_list[0]
    if hp is None:
        hp = "NULL"
    if ps is None:
        ps = "NULL"

    other_hp_support = sum(other_supp for other_supp, other_hp in hp_list if other_hp != hp and other_hp != "NULL")
    other_ps_support = sum(other_supp for other_supp, other_ps in ps_list if other_ps != ps and other_ps != "NULL")

    hp_filter = "FAIL"
    if (float(other_hp_support) / (hp_support + other_hp_support) < config.phase_conflict_threshold and hp != "NULL" and
            hp_support > 0):
        hp_filter = "PASS"

    ps_filter = "FAIL"
    if (float(other_ps_support) / (ps_support + other_ps_support) < config.phase_conflict_threshold and ps != "NULL" and
            ps_support > 0):
        ps_filter = "PASS"

    svcall.set_info("PHASE", f"{hp},{ps},{hp_support},{ps_support},{hp_filter},{ps_filter}")
    return (config.phase_identifiers.index(hp) if hp in config.phase_identifiers else None if hp_filter == "PASS" else
            None, ps if "PASS" == ps_filter else None)
