from dataclasses import dataclass
from typing import Tuple
from unittest import TestCase
from unittest.mock import MagicMock

from sniffles.cluster import Cluster
from sniffles.leadprov import Lead
from sniffles.sv import SVCall, SVCallBNDInfo, resolve_bnd


@dataclass
class TestInfo:
    id: str
    contig: str
    pos: int
    mate_contig: str
    mate_pos: int
    is_first: bool
    is_reverse: bool
    expected_alt: str = None
    expected_orientation: str = None


class TestBND(TestCase):
    """
    Tests for BND representation in VCF and CSV output
    """
    def _get_test(self, info: TestInfo) -> Tuple[SVCall, Cluster]:
        svcall = SVCall(
            contig=info.contig,
            pos=info.pos,
            id=info.id,
            ref='N',
            alt='<BND>',
            qual=60,
            filter='PASS',
            info={},
            svtype='BND',
            svlen=0,
            end=info.pos,
            genotypes={},
            precise=True,
            support=10,
            rnames=None,
            qc=True,
            nm=-1,
            postprocess=MagicMock()
        )
        cluster = Cluster(
            id=info.id[-1:],
            svtype='BND',
            contig=info.contig,
            start=info.pos,
            end=info.pos,
            seed=info.pos,
            leads=[
                ld := Lead(
                    read_id=1,
                    read_qname='read1',
                    contig=info.contig,
                    ref_start=info.pos,
                    ref_end=info.pos,
                    qry_start=1000,
                    qry_end=1000,
                    strand='+',
                    mapq=60,
                    nm=100,
                )
            ],
            repeat=False,
            leads_long=None
        )
        ld.bnd_info = SVCallBNDInfo(
            mate_contig=info.mate_contig,
            mate_ref_start=info.mate_pos,
            is_first=info.is_first,
            is_reverse=info.is_reverse,
        )
        return svcall, cluster

    def test_resolve_bnd(self) -> None:
        """
        Test BND representation. See VCF 4.2 spec, chapter 5.4
        """
        tests = [
            TestInfo(
                'bnd_W', 'chr2', 321681, 'chr17', 198982,
                is_first=True, is_reverse=True,
                expected_alt='N]chr17:198982]',
                expected_orientation='++',
            ),
            TestInfo(
                'bnd_V', 'chr2', 321682, 'chr13', 123456,
                is_first=False, is_reverse=True,
                expected_alt=']chr13:123456]N',
                expected_orientation='-+',
            ),
            TestInfo(
                'bnd_U', 'chr13', 123456, 'chr2', 321682,
                is_first=True, is_reverse=False,
                expected_alt='N[chr2:321682[',
                expected_orientation='+-',
            ),
            TestInfo(
                'bnd_X', 'chr13', 123457, 'chr17', 198983,
                is_first=False, is_reverse=False,
                expected_alt='[chr17:198983[N',
                expected_orientation='--',
            ),
            TestInfo(
                'bnd_Y', 'chr17', 198982, 'chr2', 321681,
                is_first=True, is_reverse=True,
                expected_alt='N]chr2:321681]',
                expected_orientation='++',
            ),
            TestInfo(
                'bnd_Z', 'chr17', 198983, 'chr13', 123457,
                is_first=False, is_reverse=False,
                expected_alt='[chr13:123457[N',
                expected_orientation='--',
            ),
        ]

        for ti in tests:
            with self.subTest(ti.id):
                svcall, cluster = self._get_test(ti)
                resolve_bnd(svcall, cluster)

                self.assertEqual(svcall.alt, ti.expected_alt)
                self.assertIn(f'CHR2', svcall.info)
                self.assertEqual(svcall.info['CHR2'], ti.mate_contig)
                csv_fields = svcall._to_csv_line()
                self.assertEqual(
                    csv_fields[:7],
                    ('BND', ti.expected_orientation[0], ti.contig, str(ti.pos), ti.expected_orientation[1], ti.mate_contig, str(ti.mate_pos))
                )
