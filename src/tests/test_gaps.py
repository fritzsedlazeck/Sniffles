import unittest

from tests.common import get_svcall


class TestGaps(unittest.TestCase):
    def test_Gaps(self):
        from sniffles.util import UnresolvedReferenceRegions

        unresolved = UnresolvedReferenceRegions('chr1', [(0, 10_000), (20_000, 30_000)])

        with self.subTest('Outside gap region/normal call'):
            c = get_svcall(
                pos=15_000,
                end=15_100,
            )
            self.assertEqual(0.0, unresolved.get_N_percentage(c))

        with self.subTest('Inside gap region'):
            c = get_svcall(
                pos=5_000,
                end=5_100,
            )
            self.assertEqual(100.0, unresolved.get_N_percentage(c))

        with self.subTest('Partially inside gap region'):
            c = get_svcall(
                pos=9_900,
                end=10_100,
            )
            self.assertAlmostEqual(50.0, unresolved.get_N_percentage(c), delta=0.001)

        with self.subTest('Enveloping gap region'):
            c = get_svcall(
                pos=19_990,
                end=30_010,
            )
            res = unresolved.get_N_percentage(c)
            self.assertAlmostEqual(99.8, res, delta=0.001)

        with self.subTest('Spanning multiple gaps'):
            c = get_svcall(
                pos=0,
                end=50_000,
            )
            self.assertAlmostEqual(40.0, unresolved.get_N_percentage(c), delta=0.001)

    def test_NoGapRegions(self):
        from sniffles.util import UnresolvedReferenceRegions

        unresolved = UnresolvedReferenceRegions('chr1', [])

        c = get_svcall(
            pos=0,
            end=50_000,
        )
        self.assertEqual(0.0, unresolved.get_N_percentage(c))
