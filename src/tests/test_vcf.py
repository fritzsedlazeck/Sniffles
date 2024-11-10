#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     02.01.2024
# Author:      Hermann Romanek
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#
from typing import Tuple, Callable
from unittest import TestCase
from unittest.mock import Mock

from sniffles.sv import SVCall
from sniffles.vcf import VCF


class TestVCFFormat(TestCase):
    """
    Unittests for writing VCF output
    """

    @staticmethod
    def get_config():
        """
        Generate a mock config
        """
        config = Mock()
        config.sample_ids_vcf = []
        config.output_rnames = True
        config.mosaic_af_max = 0.3
        config.mosaic = False
        config.id_prefix = 'Sniffles.'
        config.symbolic = False
        config.max_del_seq_len = 50000
        config.genotype_format = "GT:GQ:DR:DV"
        return config

    def verify_common_fields(self, *args, **kwargs) -> Tuple[int, str, str, str]:
        """
        Verifies common fields, returns pos, ref and alt fields
        """
        flds = args[0].split('\t')
        self.assertGreater(len(flds), 8)
        contig, pos, id, ref, alt, qual, filter, info = flds[:8]

        self.assertEqual('chr1', contig)
        self.assertEqual('Sniffles.unittest-1', id)

        return int(pos), ref, alt, info

    def parse_info(self, info: str) -> tuple[set[str], dict[str, str]]:
        """
        Simple parsing of an info string
        """
        flags = []
        flds = {}
        for fld in info.split(';'):
            if '=' in fld:
                key, value = fld.split('=')
                flds[key] = value
            else:
                flags.append(fld)
        return set(flags), flds

    def get_vcf(self, reference: str):
        vcf = VCF(self.get_config(), None)
        vcf.reference_handle = Mock()
        vcf.reference_handle.fetch = Mock(
            side_effect=lambda refname, start, end: reference[start:end]
        )
        return vcf

    def get_svcall(self, **kwargs):
        sv_kwargs = {
            'contig': 'chr1',
            'id': 'unittest-1',
            'qual': 10,
            'filter': 'PASS',
            'info': {},
            'genotypes': {},
            'precise': True,
            'support': 100,
            'rnames': ['ut'],
            'postprocess': None,
            'qc': True,
            'nm': -1,
            'fwd': 1,
            'rev': 1,
        }
        sv_kwargs.update(kwargs)
        return SVCall(**sv_kwargs)

    def test_spec_ins(self):
        """
        Test from VCF spec 4.2, chapter 5.2.2

        #CHROM POS ID REF ALT   QUAL FILTER INFO
        20     3   .  C   CTAG  .    PASS   DP=100

        Ref    a t C - - - g a
        1      a t C T A G g a
        """
        def verify(*args, **kwargs):
            pos, ref, alt, _ = self.verify_common_fields(*args, **kwargs)
            self.assertEqual(3, pos)
            self.assertEqual('C', ref)
            self.assertEqual('CTAG', alt)

        vcf = self.get_vcf('atCga')
        vcf.write_raw = Mock(side_effect=verify)

        vcf.write_call(self.get_svcall(
            svtype='INS',
            ref='N',
            alt='TAG',
            pos=3,
            svlen=3,
            end=3,
        ))

        vcf.reference_handle.fetch.assert_called_with('chr1', 2, 3)
        vcf.write_raw.assert_called()

    def test_spec_del(self):
        """
        Test from VCF spec 4.2, chapter 5.2.3

        #CHROM POS ID REF ALT QUAL FILTER INFO
        20     2   .  TCG T   .    PASS   DP=100

        Ref  a T C G a
        Alt  a T - - a
        """
        def verify(*args, **kwargs):
            pos, ref, alt, info = self.verify_common_fields(*args, **kwargs)
            self.assertEqual(2, pos)
            self.assertEqual('TCG', ref)
            self.assertEqual('T', alt)
            _, info_fields = self.parse_info(info)
            self.assertEqual(pos + abs(int(info_fields['SVLEN'])) - 1, int(info_fields['END']))

        vcf = self.get_vcf('aTCGa')
        vcf.write_raw = Mock(side_effect=verify)

        vcf.write_call(self.get_svcall(
            svtype='DEL',
            ref='N',
            alt='<DEL>',
            pos=2,
            svlen=-2,
            end=4,
        ))

        vcf.write_raw.assert_called()

    def test_del_issue31(self):
        """
        See https://github.com/fritzsedlazeck/Sniffles_dev/issues/31
        """
        def mock_fetch(refname, start, end):
            """
            Reference is extracted from:
            reference_file = pysam.FastaFile('GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz')
            reference_file.fetch('chr1', 964600, 964800)
            """
            reference = 'CAGTGGGGATGTGCTGCGGGGAGGGGGGCGCGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCGCGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCGCGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCGCGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCGCGGGTCCGCAGTGGGGAT'
            return reference[start-964600:end-964600]

        vcf = VCF(self.get_config(), None)
        vcf.reference_handle = Mock()
        vcf.reference_handle.fetch = Mock(
            side_effect=mock_fetch
        )

        def verify(*args, **kwargs):
            pos, ref, alt, info = self.verify_common_fields(*args, **kwargs)
            self.assertEqual(964631, pos)
            self.assertEqual('CGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCGCGGGTCCGCAGTGGGGATGTGCTGCCGGGAGGGGGGCG', ref)
            self.assertEqual('C', alt)

        vcf.write_raw = Mock(side_effect=verify)

        vcf.write_call(self.get_svcall(
            svtype='DEL',
            ref='N',
            alt='<DEL>',
            pos=964631,
            svlen=-75,
            end=964631-75,
        ))

        vcf.write_raw.assert_called()

    def test_unresolved_ins(self):
        """
        Tests output of an SV where we were unable to resolve the sequence.
          See https://github.com/fritzsedlazeck/Sniffles/issues/501
        """
        def verify(*args, **kwargs):
            pos, ref, alt, info = self.verify_common_fields(*args, **kwargs)
            self.assertEqual(2, pos)
            self.assertEqual('T', ref)
            self.assertEqual('<INS>', alt)

        vcf = self.get_vcf('T'*50)
        vcf.write_raw = Mock(side_effect=verify)

        vcf.write_call(self.get_svcall(
            svtype='INS',
            ref='N',
            alt='<INS>',
            pos=2,
            svlen=20,
            end=22,
        ))

        vcf.write_raw.assert_called()
