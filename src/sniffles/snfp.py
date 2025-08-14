#!/usr/bin/env python3
"""
Population SNF
"""
import logging
import math
from dataclasses import dataclass, asdict
from typing import Optional, Callable

from pysam import VariantFile
try:
    from edlib import align
except ImportError:
    align: Optional[Callable] = None

from sniffles.config import SnifflesConfig
# from sniffles.genotyping import Genotype
from sniffles.snf import SNFileBase
from sniffles.sv import SVCall


log = logging.getLogger(__name__)


@dataclass
class PopulationVariant:
    contig: str
    pos: int
    id: str
    alt: str

    svtype: str
    svlen: int
    end: int

    af: float  # population allele fraction
    genotyped_sample_count: int  # number of genotyped samples for this variant
    variant_sample_count: int  # number of samples in the population carrying this variant

    @staticmethod
    def _calculate_frequency(genotypes: dict[int, tuple]) -> tuple[float, int, int]:
        """
        Calculate population AF from a dict of genotypes. Keys are sample ids, values are genotype tuples.
        Returns a tuple with population AF, total number of genotyped samples and number of samples carrying this SV.
        """
        total_alleles = 0
        variant_alleles = 0
        genotyped_samples = 0
        variant_samples = 0
        gp = SnifflesConfig.GLOBAL.genotype_ploidy

        for gt in genotypes.values():
            if gt[0] == '.':
                continue
            else:
                genotyped_samples += 1
                variant_number = gt[0] + gt[1]
                total_alleles += gp
                variant_alleles += variant_number
                if variant_number > 0:
                    variant_samples += 1

        return variant_alleles / total_alleles, genotyped_samples, variant_samples

    @classmethod
    def from_svcall(cls, svcall: SVCall) -> Optional['PopulationVariant']:
        """
        Generate a population variant from a merged SV call, calculating frequencies.
        """
        af, genotyped_samples, variant_samples = cls._calculate_frequency(svcall.genotypes)
        population_size = len(SnifflesConfig.GLOBAL.snf_input_info)

        if (pct_genotyped := (genotyped_samples / population_size)) < SnifflesConfig.GLOBAL.dev_population_min_gt:
            log.warning(f'Not emitting call {svcall} due to only {genotyped_samples}/{population_size} ({pct_genotyped*100:.2f}%) genotyped samples.')
            return None

        obj = cls(
            contig=svcall.contig,
            pos=svcall.pos,
            id=svcall.id,
            alt=svcall.alt,
            svtype=svcall.svtype,
            svlen=svcall.svlen,
            end=svcall.end,
            af=af,
            genotyped_sample_count=genotyped_samples,
            variant_sample_count=variant_samples,
        )
        return obj

    def match(self, svcall: SVCall) -> int | None:
        """
        Determine if svcall is the same variant. Returns the distance (smaller being better) or None if it
        doesn't match at all.
        """
        config = SnifflesConfig.GLOBAL
        dist = abs(self.pos - svcall.pos) + abs(abs(self.svlen) - abs(svcall.svlen))
        minlen = float(min(abs(self.svlen), abs(svcall.svlen)))
        if dist > config.combine_match * math.sqrt(minlen) or dist > config.combine_match_max:
            return None

        if self.svtype == 'INS' and (limit := config.combine_pctseq):
            distance = align(self.alt, svcall.alt)['editDistance']
            if (self.svlen - distance) / self.svlen <= limit:
                return None

        return dist


@dataclass
class PopulationInfo:
    version: int
    name: str
    description: str
    size: int


class PopulationSNF(SNFileBase):
    """
    A population SNF holds information on variants and their frquency in a population
    """
    _blocks = None

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._blocks = {}

    def _calculate_block_index(self, pos: int) -> int:
        return int(pos / self.config.snf_block_size) * self.config.snf_block_size

    def get_population_AF(self, svcall: SVCall) -> tuple[float, int] | None:
        """
        Returns the population AF and population size for the given SVCall. Returns None if the variant is not found in the population SNF.
        """
        if svcall.contig not in self._blocks:
            self._blocks[svcall.contig] = self.get_all_blocks(svcall.contig)

        block = str(self._calculate_block_index(svcall.pos))
        best_dist = None
        best_variant = None
        try:
            for pv in self._blocks[svcall.contig][block][svcall.svtype]:
                pv: PopulationVariant
                dist = pv.match(svcall)
                if dist is not None:
                    if best_dist is None or dist < best_dist:
                        best_dist = dist
                        best_variant = pv
        except KeyError:
            ...
        else:
            if best_variant is not None:
                return round(best_variant.af, 5), best_variant.genotyped_sample_count

        return None

    def _create_header(self, config: SnifflesConfig, main_index: dict, snf_candidate_count: int) -> dict:
        """
        Extend the standard snf header with population information
        """
        d = super()._create_header(config, main_index, snf_candidate_count)
        d['population'] = asdict(PopulationInfo(
            version=1,
            name='Population',
            description='A sample population',
            size=len(config.snf_input_info)
        ))
        return d

    def read_header(self):
        """
        Deserialize population information
        """
        super().read_header()
        try:
            self.header['population'] = PopulationInfo(**self.header['population'])
        except:  # noqa
            log.warning(f'Unable to deserialize population information from SNF header.', exc_info=True)

    def _calculate_contig_coverages(self, *args, **kwargs) -> dict:
        """
        Coverages are not used in population SNFs.
        """
        return {}

    def store(self, svcand: SVCall) -> None:
        """
        Convert SVCall to PopulationVariant and store it
        """
        if (variant := PopulationVariant.from_svcall(svcand)) is not None:
            super().store(variant)
        return variant is not None


class PopulationVCF:
    def __init__(self, vcf_path: str):
        self.vcf = VariantFile(vcf_path)

    def get_population_annotation(self, sv: SVCall):
        self.vcf.fetch(sv.contig, sv.pos - 500, sv.pos + sv.svlen + 500)
