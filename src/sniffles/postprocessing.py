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
import math
import logging

from sniffles import util
from sniffles import consensus
from sniffles.config import SnifflesConfig
from sniffles.leadprov import LeadProvider
from sniffles.sv import SVCall

log = logging.getLogger('sniffles.postprocessing')


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


def coverage(calls: list[SVCall], lead_provider: 'LeadProvider') -> float:
    """
    Annotate the given SVCalls with coverage information based on the lead provider's coverage data.
    Returns the average coverage across the contig.

    INS/BND    |-------------V---------------|
              U---       S--- E---            D---
                            C---
    DUP/DEL    |------------===============------------|
    INV       U---          S---       E---            D---

    """
    cv = lead_provider.coverage
    config = SnifflesConfig.GLOBAL

    for svcall in calls:
        start = svcall.pos
        if svcall.svtype == "INS":
            end = start + 1
        else:
            end = svcall.pos + abs(svcall.svlen)

        if svcall.svtype in ("INS", "BND"):
            try:
                svcall.coverage_start = int(cv[start - config.coverage_binsize])
            except IndexError:
                logging.debug(f'Coverage index out of range: start {start - config.coverage_binsize} outside 0..{len(cv)} for {svcall}')
            try:
                svcall.coverage_center = int(cv[start])
            except IndexError:
                logging.debug(f'Coverage index out of range: center {start} outside 0..{len(cv)} for {svcall}')
            try:
                svcall.coverage_end = int(cv[end + config.coverage_binsize])
            except IndexError:
                logging.debug(f'Coverage index out of range: end {end + config.coverage_binsize} outside 0..{len(cv)} for {svcall}')
        else:
            try:
                svcall.coverage_start = int(cv[start])
            except IndexError:
                logging.debug(f'Coverage index out of range: start {start} outside 0..{len(cv)} for {svcall}')
            try:
                svcall.coverage_center = int(cv[int((start + end) / 2)])
            except IndexError:
                logging.debug(f'Coverage index out of range: center {int((start + end) / 2)} outside 0..{len(cv)} for {svcall}')
            try:
                svcall.coverage_end = int(cv[end - config.coverage_binsize])
            except IndexError:
                logging.debug(f'Coverage index out of range: end {end - config.coverage_binsize} outside 0..{len(cv)} for {svcall}')

        try:
            svcall.coverage_upstream = int(cv[start - config.coverage_binsize * config.coverage_updown_bins])
        except IndexError:
            logging.debug(f'Coverage index out of range: upstream {start - config.coverage_binsize * config.coverage_updown_bins} outside 0..{len(cv)} for {svcall}')
        try:
            svcall.coverage_downstream = int(cv[end + config.coverage_binsize * config.coverage_updown_bins])
        except IndexError:
            logging.debug(f'Coverage index out of range: downstream {end + config.coverage_binsize * config.coverage_updown_bins} outside 0..{len(cv)} for {svcall}')

    return lead_provider.coverage.mean()


def qc_sv_support(svcall, coverage_global, config) -> bool:
    dev_sv_filter: list = []
    if config.dev_filter:
        if type(svcall.filter) is type(""):
            if "-" in svcall.filter:
                dev_sv_filter = svcall.filter.split("-")
            else:
                dev_sv_filter = [svcall.filter]

    if config.minsupport == "auto":
        if not qc_support_auto(svcall, coverage_global, config):
            if config.dev_filter:
                dev_sv_filter.append("SUPPORT_MIN")
            else:
                svcall.filter = "SUPPORT_MIN"
                return False
    else:
        if not qc_support_const(svcall, config):
            if config.dev_filter:
                dev_sv_filter.append("SUPPORT_MIN")
            else:
                svcall.filter = "SUPPORT_MIN"
                return False

    if config.dev_filter:
        svcall.filter = "-".join(dev_sv_filter)
    return True


def rescale_support(svcall, config: SnifflesConfig) -> int:
    """
    Rescale support for long insertions.
    """
    if svcall.svtype != "INS" or svcall.svlen < config.long_ins_length:
        return svcall.support
    else:
        base = svcall.support
        scale_factor = config.long_ins_rescale_mult * (float(svcall.svlen) / config.long_ins_length)
        return round(base * (config.long_ins_rescale_base + scale_factor))


def qc_support_auto(svcall, coverage_global, config):
    support = rescale_support(svcall, config)

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
    coverage_ = (coverage_regional * config.minsupport_auto_regional_coverage_weight +
                 coverage_global * coverage_global_weight)
    min_support = round(config.minsupport_auto_base + config.minsupport_auto_mult * coverage_)
    return support >= min_support


def qc_support_const(svcall, config):
    # svcall.set_info("MINSUPPORT",config.minsupport)
    return svcall.support >= config.minsupport


def qc_sv(svcall: SVCall, config: SnifflesConfig):
    dev_sv_filter: list = []
    if config.dev_filter:
        if type(svcall.filter) is type(""):
            if "-" in svcall.filter:
                dev_sv_filter = svcall.filter.split("-")
            else:
                dev_sv_filter = [svcall.filter]

    if config.qc_stdev:
        stdev_pos = svcall.get_info("STDEV_POS")
        if stdev_pos > config.qc_stdev_abs_max:
            if config.dev_filter:
                dev_sv_filter.append("STDEV_POS")
            else:
                svcall.filter = "STDEV_POS"
                return False
        if svcall.svtype not in ("BND", "SINGLE_LEFT", "SINGLE_RIGHT") and stdev_pos / abs(svcall.svlen) > 2.0:
            if config.dev_filter:
                dev_sv_filter.append("STDEV_POS")
            else:
                svcall.filter = "STDEV_POS"
        if svcall.svtype not in ("BND", "SINGLE_LEFT", "SINGLE_RIGHT") and stdev_pos / abs(svcall.svlen) > 2.0:
            svcall.filter = f'{svcall.filter}-STDEV_POS' if config.dev_filter else "STDEV_POS"
            if not config.dev_filter:
                return False

        stdev_len = svcall.get_info("STDEV_LEN")
        if stdev_len is not None and stdev_len != 0:
            if svcall.svtype != "BND" and stdev_len / abs(svcall.svlen) > 1.0:
                if config.dev_filter:
                    dev_sv_filter.append("STDEV_LEN")
                else:
                    svcall.filter = "STDEV_LEN"
                    return False
            if stdev_len > config.qc_stdev_abs_max:
                if config.dev_filter:
                    dev_sv_filter.append("STDEV_LEN")
                else:
                    svcall.filter = "STDEV_LEN"
                    return False

    if svcall.is_single_break:
        svcall.filter = 'SINGLE_BREAK'
        return False

    support_overwrite_svlen = 10
    if abs(svcall.svlen) < config.minsvlen:
        if svcall.support < support_overwrite_svlen:
            if config.dev_filter:
                dev_sv_filter.append("SVLEN_MIN")
            else:
                svcall.filter = "SVLEN_MIN"
                return False

    if svcall.svtype == "BND":
        if config.qc_bnd_filter_strand and len(set(lead.strand for lead in svcall.postprocess.cluster.leads)) < 2:
            if config.dev_filter:
                dev_sv_filter.append("STRAND_BND")
            else:
                svcall.filter = "STRAND_BND"
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
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DEL")
                    else:
                        svcall.filter = "COV_CHANGE_DEL"
                        return False
            elif svcall.coverage_upstream < svcall.coverage_center < svcall.coverage_downstream:
                if svcall.coverage_upstream/svcall.coverage_downstream < upstream_downstream_max_coverage_diff:
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DEL")
                    else:
                        svcall.filter = "COV_CHANGE_DEL"
                        return False
            else:
                pass
        if svcall.coverage_upstream > svcall.coverage_downstream:
            if (upstream_downstream_diff > svcall.coverage_downstream/svcall.coverage_upstream or
                    svcall.coverage_center > svcall.coverage_downstream):
                if config.dev_filter:
                    dev_sv_filter.append("COV_CHANGE_DEL")
                else:
                    svcall.filter = "COV_CHANGE_DEL"
                    return False
        elif svcall.coverage_upstream < svcall.coverage_downstream:
            if (upstream_downstream_diff > svcall.coverage_upstream/svcall.coverage_downstream or
                    svcall.coverage_upstream < svcall.coverage_center):
                if config.dev_filter:
                    dev_sv_filter.append("COV_CHANGE_DEL")
                else:
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
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DUP")
                    else:
                        svcall.filter = "COV_CHANGE_DUP"
                        return False
            elif svcall.coverage_upstream < svcall.coverage_center < svcall.coverage_downstream:
                if svcall.coverage_upstream/svcall.coverage_downstream < upstream_downstream_max_coverage_diff:
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DUP")
                    else:
                        svcall.filter = "COV_CHANGE_DUP"
                        return False
            else:
                pass
            if svcall.coverage_upstream > svcall.coverage_downstream:
                if (upstream_downstream_diff > svcall.coverage_downstream / svcall.coverage_upstream or
                        svcall.coverage_center < svcall.coverage_downstream):
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DUP")
                    else:
                        svcall.filter = "COV_CHANGE_DUP"
                        return False
            elif svcall.coverage_upstream < svcall.coverage_downstream:
                if (upstream_downstream_diff > svcall.coverage_upstream / svcall.coverage_downstream or
                        svcall.coverage_upstream > svcall.coverage_center):
                    if config.dev_filter:
                        dev_sv_filter.append("COV_CHANGE_DUP")
                    else:
                        svcall.filter = "COV_CHANGE_DUP"
                        return False
            else:
                pass
    elif svcall.svtype == "INS" and (svcall.coverage_upstream < config.qc_coverage or
                                     svcall.coverage_downstream < config.qc_coverage):
        if config.dev_filter:
            dev_sv_filter.append("COV_CHANGE_INS")
        else:
            svcall.filter = "COV_CHANGE_INS"
            return False

    qc, val = svcall.qc_coverage_samples()
    svcall.set_info('COVERAGE_VAR', val)
    if not qc:
        svcall.filter = f'COV_VAR'
        if not config.dev_filter:
            return False

    qc_coverage_max_change_frac = config.qc_coverage_max_change_frac
    # if config.mosaic and sv_is_mosaic:
    #     qc_coverage_max_change_frac = config.mosaic_qc_coverage_max_change_frac
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
            if config.dev_filter:
                dev_sv_filter.append("COV_CHANGE_FRAC_US")
            else:
                svcall.filter = "COV_CHANGE_FRAC_US"
                return False

        if abs(s - c) / max(s, c) > qc_coverage_max_change_frac:
            if config.dev_filter:
                dev_sv_filter.append("COV_CHANGE_FRAC_SC")
            else:
                svcall.filter = "COV_CHANGE_FRAC_SC"
                return False

        if abs(c - e) / max(c, e) > qc_coverage_max_change_frac:
            if config.dev_filter:
                dev_sv_filter.append("COV_CHANGE_FRAC_CE")
            else:
                svcall.filter = "COV_CHANGE_FRAC_CE"
                return False

        if abs(e - d) / max(e, d) > qc_coverage_max_change_frac:
            if config.dev_filter:
                dev_sv_filter.append("COV_CHANGE_FRAC_ED")
            else:
                svcall.filter = "COV_CHANGE_FRAC_ED"
                return False

    if config.dev_filter:
        svcall.filter = "-".join(dev_sv_filter)
    return True


def qc_sv_post_annotate(svcall: SVCall, config: SnifflesConfig, coverage_average_total: float):
    dev_sv_filter: list = []
    if config.dev_filter:
        if type(svcall.filter) is type(""):
            if "-" in svcall.filter:
                dev_sv_filter = svcall.filter.split("-")
            else:
                dev_sv_filter = [svcall.filter]

    af = svcall.get_info("VAF")
    af = af if af is not None else 0
    sv_is_mosaic = af <= config.mosaic_af_max

    if (svcall.coverage_center < config.qc_coverage and
            (len(svcall.genotypes) == 0 or (svcall.genotypes[0][0] != "." and
                                            svcall.genotypes[0][0] + svcall.genotypes[0][1] < 2))):
        if config.dev_filter:
            dev_sv_filter.append("COV_MIN_GT")
        else:
            svcall.filter = "COV_MIN_GT"
            return False

    if config.mosaic and not sv_is_mosaic:
        if not qc_sv_support(svcall, coverage_average_total, config):
            if not config.dev_filter:
                return False

    qc_nm = config.qc_nm
    qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if config.mosaic and sv_is_mosaic:
        qc_nm = config.mosaic_qc_nm
        qc_nm_threshold = config.qc_nm_threshold * config.qc_nm_mult
    if qc_nm and svcall.nm > qc_nm_threshold and (len(svcall.genotypes) == 0 or svcall.genotypes[0][1] == 0):
        if config.dev_filter:
            dev_sv_filter.append("ALN_NM")
        else:
            svcall.filter = "ALN_NM"
            return False

    if not config.mosaic and sv_is_mosaic:
        if config.dev_filter:
            dev_sv_filter.append("MOSAIC_VAF")
        else:
            svcall.filter = "MOSAIC_VAF"
            return False

    if config.mosaic and sv_is_mosaic:
        # SUPPORT
        stdev_pos = svcall.info["STDEV_POS"] if "STDEV_POS" in svcall.info else None
        stdev_len = svcall.info["STDEV_LEN"] if "STDEV_LEN" in svcall.info else None
        svlen = svcall.info["SVLEN"] if "SVLEN" in svcall.info else 1
        min_mosaic_support = config.mosaic_min_reads
        max_stdev_to_svlen_ratio = 0.1
        max_stdev_pos_difference = 5
        if stdev_pos is not None and stdev_len is not None:
            filter_low_supp = ((not svcall.precise or stdev_len/abs(svcall.svlen) > max_stdev_to_svlen_ratio or
                                stdev_pos > max_stdev_pos_difference) and abs(svlen) <= config.max_svlen_mosaic)
            min_mosaic_support = config.mosaic_min_reads if (svcall.svtype in ["BND", "INV"] or filter_low_supp) \
                else config.mosaic_min_reads-1
        # if svcall.support < config.mosaic_min_reads:
        if svcall.support < min_mosaic_support:
            if config.dev_filter:
                dev_sv_filter.append("SUPPORT_MIN")
            else:
                svcall.filter = "SUPPORT_MIN"
                return False
        # SVLEN
        if "BND" != svcall.svtype:
            if abs(svcall.svlen) > config.max_svlen_mosaic:
                if config.dev_filter:
                    dev_sv_filter.append("SVLEN_MAX_MOSAIC")
                else:
                    svcall.filter = "SVLEN_MAX_MOSAIC"
                    return False

    if svcall.svtype != "BND":
        if not (config.mosaic and sv_is_mosaic) and config.qc_strand:
            is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
            if not is_long_ins and len(set(leads.strand for leads in svcall.postprocess.cluster.leads)) < 2:
                if config.dev_filter:
                    dev_sv_filter.append("STRAND")
                else:
                    svcall.filter = "STRAND"
                    return False
        elif (config.mosaic and sv_is_mosaic) and config.mosaic_qc_strand:
            is_long_ins = (svcall.svtype == "INS" and svcall.svlen >= config.long_ins_length)
            if (not is_long_ins and len(set(leads.strand for leads in svcall.postprocess.cluster.leads)) < 2
                    and svcall.support >= config.mosaic_use_strand_thresholds):
                if config.dev_filter:
                    dev_sv_filter.append("STRAND_MOSAIC")
                else:
                    svcall.filter = "STRAND_MOSAIC"
                    return False

    if config.mosaic and sv_is_mosaic:
        if (svcall.svtype == "INV" or svcall.svtype == "DUP") and svcall.svlen < config.mosaic_qc_invdup_min_length:
            if config.dev_filter:
                dev_sv_filter.append("SVLEN_MIN_MOSAIC")
            else:
                svcall.filter = "SVLEN_MIN_MOSAIC"
                return False

    if svcall.coverage_center < config.qc_coverage and svcall.svtype not in ["DEL", "INS"]:
        if (svcall.svtype == "INV" and svcall.svlen) > config.long_inv_length and not (config.mosaic and sv_is_mosaic):
            pass
        else:
            if config.dev_filter:
                dev_sv_filter.append("COV_MIN")
            else:
                svcall.filter = "COV_MIN"
                return False

    if config.mosaic:
        if sv_is_mosaic and (af < config.mosaic_af_min or af > config.mosaic_af_max):
            if config.dev_filter:
                dev_sv_filter.append("MOSAIC_VAF")
            else:
                svcall.filter = "MOSAIC_VAF"
                return False
        elif not sv_is_mosaic and not config.mosaic_include_germline:
            if config.dev_filter:
                dev_sv_filter.append("NOT_MOSAIC_VAF")
            else:
                svcall.filter = "NOT_MOSAIC_VAF"
                return False

    if config.dev_filter and "PASS" != svcall.filter:
        svcall.filter = "-".join(dev_sv_filter)
    return True


def binomial_coef(n, k):
    return math.factorial(n) / (math.factorial(k) * math.factorial(n - k))


def genotype_sv(svcall: SVCall, config, phase: tuple | None = None):
    from sniffles.genotyping import GENOTYPER_BY_TYPE, Genotyper

    GENOTYPER_BY_TYPE.get(svcall.svtype, Genotyper)(svcall, config, phase).calculate()

    # post HP assessment for GT, hom alt should skip hp_filter, but gt is after phase
    try:
        a, b, gq, dr, dv, phase = svcall.genotypes[0]
        if a == b and a == 1 and (phase := svcall.get_info("PHASE")):
            hp, ps, hp_supp, ps_supp, hp_filt, ps_filt = phase.split(",")
            if "NULL" != hp and "NULL" != ps:
                hp_filt, ps_filt = "PASS", "PASS"
                phase = (hp, ps)
                svcall.genotypes[0] = (a, b, gq, dr, dv, phase)
                svcall.set_info("PHASE", f"{hp},{ps},{hp_supp},{ps_supp},{hp_filt},{ps_filt}")
    except KeyError:
        pass


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
    hp_return = config.phase_identifiers.index(hp) if hp in config.phase_identifiers and hp_filter == "PASS" else None
    ps_return = ps if "PASS" == ps_filter else None
    return hp_return, ps_return
