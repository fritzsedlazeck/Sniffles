"""
These tests check Sniffles BND calls/leads against the HG008 and HG002 dataset.

Requires the BAM file to be present at the path specified in BNDLeadTestCase.bam_file. The BAM file should be indexed and contain the necessary tags (e.g. NM).
"""
from typing import Tuple
from unittest import TestCase

import pysam

from sniffles.config import SnifflesConfig
from sniffles.leadprov import Lead


class BNDLeadTestCase(TestCase):
    """
    Base class for testing leads of BNDs.
    """
    bam_file = 'data/hg008.bam'

    @classmethod
    def _get_common_args(cls) -> Tuple:
        """
        Get common required arguments for sniffles
        """
        return '--input', cls.bam_file, '--vcf', 'out.vcf'

    @classmethod
    def setUpClass(cls):
        cls.config = SnifflesConfig(*cls._get_common_args())
        cls.bam = pysam.AlignmentFile(cls.bam_file)

    def _get_read(self, contig: str, pos: int, query_name: str) -> pysam.AlignedSegment:
        region = f'{contig}:{pos}-{pos+1}' if not self.bam_file.endswith('.sam') else None
        for read in self.bam.fetch(region=region):
            if read.query_name == query_name:
                return read
        raise ValueError(f'Read {query_name} not found in region {region}')

    def _asserts(self, lead: Lead):
        ...

    def do_test(self, read: pysam.AlignedSegment):
        with self.subTest('New version'):
            lead = Lead.for_bnd(0, read)
            self._asserts(lead)
        

class TestBNDLeadsOrange(BNDLeadTestCase):
    """
    Truth:
    chr1    23272628        SV_1    G       G]chr5:52747359]
    """
    def _asserts(self, lead):
        print(lead)
        print(lead.bnd_info)
        self.assertEqual(lead.contig, 'chr1')
        self.assertEqual(lead.ref_start, 23_272_628)
        self.assertEqual(lead.bnd_info.mate_contig, 'chr5')
        self.assertEqual(lead.bnd_info.mate_ref_start, 52_747_359)
        self.assertTrue(lead.bnd_info.is_first)
        self.assertTrue(lead.bnd_info.is_reverse)

    def test_LeadsPrimaryForward(self):
        read = self._get_read(
            'chr1',
            23_272_628,
            'fcdb7746-5405-4548-9d72-3a0c81903e1c'
        )
        self.assertFalse(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsPrimaryReverse(self):
        """
        Lead for a primary read on the reverse strand
        """
        read = self._get_read(
            'chr1',
            23_272_628,
            '4c68b01d-b732-49f3-9e4a-6f1594ac5f0a'
        )
        self.assertFalse(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryForward(self):
        """
        Lead for a supplementary read on the forward strand
        """
        read = self._get_read(
            'chr1',
            23_272_628,
            '5089c480-4eae-4c61-87f8-7278dea0daaa'
        )
        self.assertTrue(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryReverse(self):
        """
        Lead for a supplementary read on the reverse strand
        """
        read = self._get_read(
            'chr1',
            23_272_628,
            '5647a0ed-80f2-4c6f-bbe4-937d95ac327b'
        )
        self.assertTrue(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)


class TestBNDLeadsGreen(BNDLeadTestCase):
    """
    Test leads of BNDs in the HG008 dataset, case green

    Truth:
    chr18   21493610        SV_136  T       [chr20:25499120[T
    """
    def _asserts(self, lead):
        print(lead)
        print(lead.bnd_info)
        self.assertEqual(lead.contig, 'chr18')
        self.assertEqual(lead.ref_start, 21_493_610)
        self.assertEqual(lead.bnd_info.mate_contig, 'chr20')
        self.assertEqual(lead.bnd_info.mate_ref_start, 25_499_120)
        self.assertFalse(lead.bnd_info.is_first)
        self.assertTrue(lead.bnd_info.is_reverse)

    def test_LeadsPrimaryForward(self):
        """
        Lead for a primary read on the forward strand
        """
        read = self._get_read(
            'chr18',
            21_493_610,
            '7c40fcdd-2d5a-4302-aead-a5ed5bd452a3'
        )
        self.assertFalse(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsPrimaryReverse(self):
        """
        Lead for a primary read on the reverse strand
        """
        read = self._get_read(
            'chr18',
            21_493_610,
            '7297cbb7-714c-4586-998a-017051004b25'
        )
        self.assertFalse(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryForward(self):
        """
        Lead for a supplementary read on the forward strand
        """
        read = self._get_read(
            'chr18',
            21_493_610,
            '42353033-1bbd-4a0c-84dc-cbd6068295f3'
        )
        self.assertTrue(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryReverse(self):
        """
         Lead for a supplementary read on the reverse strand
        """
        read = self._get_read(
            'chr18',
            21_493_610,
            '90398957-a526-49ad-be1b-2665c1b8189e'
        )
        self.assertTrue(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)


class TestBNDLeadsRedLeft(BNDLeadTestCase):
    """
    Test Leads of BNDs in the HG008 dataset, case red
    - strand of both sides of the BND has to be the same

    Truth:
      Left side:  (should be 28481424?)
      - chr18   28481423        SV_138  C       C[chrX:95812869[
      Right side:
      - chrX    95812869        SV_204  G       ]chr18:28481423]G
    """
    def _asserts(self, lead):
        print(lead)
        print(lead.bnd_info)
        self.assertEqual(lead.contig, 'chr18')
        self.assertEqual(lead.ref_start, 28_481_424)
        self.assertEqual(lead.bnd_info.mate_contig, 'chrX')
        self.assertEqual(lead.bnd_info.mate_ref_start, 95812869)
        self.assertTrue(lead.bnd_info.is_first)
        self.assertFalse(lead.bnd_info.is_reverse)

    def test_LeadsPrimaryForward(self):
        read = self._get_read(
            'chr18',
            28_481_423,
            '49485b61-facf-4f8b-81ab-4ff0f1241ec8'
        )
        self.assertFalse(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsPrimaryReverse(self):
        """
        Lead for a primary read on the reverse strand
        """
        read = self._get_read(
            'chr18',
            28_481_423,
            '48d9d042-886f-41e5-916c-77a52bd75f29'
        )
        self.assertFalse(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryForward(self):
        """
        Lead for a supplementary read on the forward strand
        """
        read = self._get_read(
            'chr18',
            28_481_423,
            '04920d3b-9413-4c38-9394-9a888bb7f6cb'
        )
        self.assertTrue(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryReverse(self):
        """
         Lead for a supplementary read on the reverse strand
        """
        read = self._get_read(
            'chr18',
            28_481_423,
            '4812c8e2-daa8-440c-be1f-7bb15f87b99a'
        )
        self.assertTrue(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)


class TestBNDLeadsRedRight(BNDLeadTestCase):
    """
    Truth:
    chrX    95812869        SV_204  G       ]chr18:28481423]G
    """
    def _asserts(self, lead):
        print(lead)
        print(lead.bnd_info)
        self.assertEqual(lead.contig, 'chrX')
        self.assertEqual(lead.ref_start, 95812869)
        self.assertEqual(lead.bnd_info.mate_contig, 'chr18')
        self.assertEqual(lead.bnd_info.mate_ref_start, 28_481_424)
        self.assertFalse(lead.bnd_info.is_first)
        self.assertTrue(lead.bnd_info.is_reverse)

    def test_LeadsPrimaryForward(self):
        read = self._get_read(
            'chrX',
            95812869,
            '04920d3b-9413-4c38-9394-9a888bb7f6cb'
        )
        self.assertFalse(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsPrimaryReverse(self):
        """
        Lead for a primary read on the reverse strand
        """
        read = self._get_read(
            'chrX',
            95812869,
            '4812c8e2-daa8-440c-be1f-7bb15f87b99a'
        )
        self.assertFalse(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryForward(self):
        """
        Lead for a supplementary read on the forward strand
        """
        read = self._get_read(
            'chrX',
            95812869,
            '49485b61-facf-4f8b-81ab-4ff0f1241ec8'
        )
        self.assertTrue(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)

    def test_LeadsSupplementaryReverse(self):
        """
        Lead for a supplementary read on the reverse strand
        """
        read = self._get_read(
            'chrX',
            95812869,
            '48d9d042-886f-41e5-916c-77a52bd75f29'
        )
        self.assertTrue(read.is_supplementary)
        self.assertTrue(read.is_reverse)

        self.do_test(read)


class TestBNDLeadsRedRightHG002(BNDLeadTestCase):
    """
    Test a BND lead in the HG002 dataset, case red right
        chr1:72.346.157
    """
    bam_file = 'data/hg002.bam'

    def _asserts(self, lead):
        print(lead)
        print(lead.bnd_info)
        self.assertEqual(lead.contig, 'chr1')
        self.assertEqual(lead.ref_start, 72_346_157)
        self.assertEqual(lead.bnd_info.mate_contig, 'chr1')
        self.assertEqual(lead.bnd_info.mate_ref_start, 72_300_641)
        self.assertFalse(lead.bnd_info.is_first)
        self.assertTrue(lead.bnd_info.is_reverse)

    def test_LeadsPrimaryForward(self):
        """
        Lead for a primary read on the forward strand
        """
        read = self._get_read(
            'chr1',
            72_346_157,
            '1a370ebb-0928-48e1-b8d3-ae8473e35654'
        )
        self.assertFalse(read.is_supplementary)
        self.assertFalse(read.is_reverse)

        self.do_test(read)
