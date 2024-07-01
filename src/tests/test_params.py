from unittest import TestCase

from sniffles.region import Region
from tests.common import get_config


class TestParams(TestCase):

    def test_contig_processing(self):
        from sniffles.util import should_process_contig

        cfg = get_config()

        with self.subTest(f'Normal contig'):
            self.assertTrue(should_process_contig("chr1", 248956422, cfg))

        with self.subTest(f'Exclude short contig'):
            self.assertFalse(should_process_contig("fragment", 123456, cfg))

        with self.subTest(f'Include short contig given by -c'):
            self.assertTrue(should_process_contig("fragment", 123456, get_config('-c', 'fragment')))

        with self.subTest(f'Include short contig given by regions'):
            cfg.regions_by_contig["fragment"] = Region("fragment", 0, 123456)
            self.assertTrue(should_process_contig("fragment", 123456, cfg))

