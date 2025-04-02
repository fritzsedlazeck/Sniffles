#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     20.09.2024
# Author:      Hermann Romanek
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
"""
Genotyping
"""
import math
from dataclasses import dataclass
from typing import Any

from sniffles.postprocessing import rescale_support
from sniffles.sv import SVCall


class UnknownGenotypeError(Exception):
    """
    Unable to determine genotype
    """


def binomial_probability(k, n, p):
    try:
        # Binomial coef cancels out for likelihood ratios
        return (p ** k) * ((1.0 - p) ** (n - k))
    except OverflowError:
        return 1.0


def likelihood_ratio(q1, q2):
    if q1 / q2 > 0:
        try:
            return math.log(q1 / q2, 10)
        except ValueError:
            return 0
    else:
        return 0


class UnknownGenotype:
    ...


@dataclass
class Genotype:
    a: int
    b: int
    qual: int  # GQ, 0-60
    dr: int
    dv: int
    phase: Any

    UNKNOWN = UnknownGenotype()



class Genotyper:
    """
    Generic genotyping
    """
    _support: int
    _coverage: float

    def __init__(self, svcall: SVCall, config, phase: tuple | None):
        self.svcall = svcall
        self.config = config
        self.phase = phase if phase is not None else self._get_phase()

    def _get_phase(self) -> tuple | None:
        """
        Get phasing information from the genotyped SV, or None if it can't be determined.
        """
        try:
            return self.svcall.genotypes[0][5]
        except (KeyError, IndexError):
            return None

    def _calculate_support(self) -> int:
        return self.svcall.support

    def _calculate_coverage(self, support: int) -> int:
        return self._get_coverage_from_list()

    def _calculate_af(self, support: int, coverage: int) -> float:
        """
        Calculate allele fraction
        """
        return support / float(coverage)

    def _get_coverage_from_list(self, coverage_list: list = None) -> int:
        """
        Coverage here is NOT the same as in IGV, but the number of reads spanning the SV
        """
        svcall = self.svcall
        if coverage_list is None:
            coverage_list = [svcall.coverage_start, svcall.coverage_center, svcall.coverage_end]

        coverage_list = [each_coverage for each_coverage in coverage_list if each_coverage != 0]

        if len(coverage_list) > 0:
            if None in coverage_list:
                new_coverage_list = [cov_value for cov_value in coverage_list if cov_value is not None]
                if len(new_coverage_list) > 0:
                    return round(sum(new_coverage_list) / len(new_coverage_list))
                else:
                    raise UnknownGenotypeError()
            else:
                return round(sum(coverage_list) / len(coverage_list))
        else:
            raise UnknownGenotypeError()

    def _filter_by_z_score(self, z_score: float) -> bool:
        """
        Should this call be filtered due to z score?
        """
        return z_score < self.config.genotype_min_z_score and not self.config.mosaic

    def calculate(self):
        config = self.config
        normalization_target = 250
        hom_ref_p = config.genotype_error
        het_p = (1.0 / config.genotype_ploidy)  # - config.genotype_error
        hom_var_p = 1.0 - config.genotype_error
        svcall = self.svcall

        support = self._calculate_support()
        try:
            coverage = self._calculate_coverage(support)
        except UnknownGenotypeError:
            return

        if support > coverage:
            coverage = support

        af = self._calculate_af(support, coverage)

        genotype_p = [((0, 0), hom_ref_p),
                      ((0, 1), het_p),
                      ((1, 1), hom_var_p)]

        max_lead = max(support, coverage)
        if max_lead > normalization_target:
            norm = normalization_target / float(max_lead)
            normalized_support = round(support * norm)
            normalized_coverage = round(coverage * norm)
        else:
            normalized_support = support
            normalized_coverage = coverage

        genotype_likelihoods = []
        for gt, p in genotype_p:
            q = binomial_probability(normalized_support, normalized_coverage, p)
            genotype_likelihoods.append((gt, q))
        genotype_likelihoods.sort(key=lambda k: k[1], reverse=True)

        sum_likelihoods = sum(q for gt, q in genotype_likelihoods)
        normalized_likelihoods = [(gt, (q / sum_likelihoods)) for gt, q in genotype_likelihoods]

        gt1, q1 = normalized_likelihoods[0]
        gt2, q2 = normalized_likelihoods[1]
        qz = [q for gt, q in normalized_likelihoods if gt == (0, 0)][0]
        genotype_z_score = min(60, int((-10) * likelihood_ratio(qz, q1)))
        genotype_quality = min(60, int((-10) * likelihood_ratio(q2, q1)))

        if svcall.filter == "PASS" and self._filter_by_z_score(genotype_z_score):
            svcall.filter = "GT"
            svcall.qc = False

        a, b = gt1
        svcall.genotypes[0] = (a, b, genotype_quality, coverage - support, support, self.phase)
        svcall.set_info("VAF", af)


class InsertionGenotyper(Genotyper):
    def _calculate_support(self):
        """
        Rescale support for long insertions. This skews our support way higher than the number of reads we have!
        """
        return rescale_support(self.svcall, self.config)

    def _calculate_coverage(self, coverage_list: list = None) -> float:
        return self._get_coverage_from_list([self.svcall.coverage_center])

    def _filter_by_z_score(self, z_score: float) -> bool:
        """
        When detecting large insertions, we allow for below threshold z scores
        """
        flt = super()._filter_by_z_score(z_score)
        if flt and self.svcall.svlen >= self.config.long_ins_length and self.config.detect_large_ins:
            return False
        return flt


class DuplicationGenotyper(Genotyper):
    def _calculate_coverage(self, support: int) -> float:
        # Experimental other coverage calculation?
        # if False and svcall.coverage_start != 0 and svcall.coverage_end != 0:
        #     if svcall.coverage_start > svcall.coverage_end:
        #         coverage_list = [svcall.coverage_end]
        #     else:
        #         coverage_list = [svcall.coverage_start]
        svcall = self.svcall
        return self._get_coverage_from_list([svcall.coverage_start, svcall.coverage_end]) + round(support * 0.75)


class InversionGenotyper(Genotyper):

    def _calculate_coverage(self, support: int) -> int:
        svcall = self.svcall
        # check event length, for whole chromosome event do something different
        return self._get_coverage_from_list([svcall.coverage_upstream, svcall.coverage_downstream]) + round(support * 0.5)


class DeletionGenotyper(Genotyper):

    def _calculate_coverage(self, support: int) -> int:
        svcall = self.svcall
        if support_sa := svcall.get_info('SUPPORT_SA'):
            return self._get_coverage_from_list([svcall.coverage_start + support_sa, svcall.coverage_center + support_sa, svcall.coverage_end + support_sa])
        else:
            return super()._calculate_coverage(support)


GENOTYPER_BY_TYPE = {
    'INS': InsertionGenotyper,
    'DEL': DeletionGenotyper,
    'DUP': DuplicationGenotyper,
    'INV': InversionGenotyper,
}
