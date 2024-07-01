import logging
import unittest
from typing import Tuple
from unittest.mock import patch, mock_open

from sniffles.config import SnifflesConfig
from sniffles.region import Region


class TestRegions(unittest.TestCase):
    """
    Tests --regions parameter
    """
    @staticmethod
    def _get_common_args() -> Tuple:
        """
        Get common required arguments for sniffles
        """
        return '--input', 'input.bam', '--vcf', 'out.vcf'

    def test_GoodFile(self):
        data = """
# comment line is ok
chr1\t100\t200\n
chr1\t500\t600\n
chr3\t500\t600\n
        """
        with patch("builtins.open", mock_open(read_data=data)) as mock_file:
            config = SnifflesConfig(*self._get_common_args(), '--regions', 'regions.bed')

        self.assertDictEqual(
            config.regions_by_contig,
            {'chr1': [Region('chr1', 100, 200), Region('chr1', 500, 600)], 'chr3': [Region('chr3', 500, 600)]}
        )
        mock_file.assert_called_with("regions.bed", "r")

    def test_ContigConflict(self):
        """
        Expect exception if both --contig and --regions are provided
        """
        with self.assertRaises(SystemExit):
            SnifflesConfig(*self._get_common_args(), '--regions', 'regions.bed', '-c', 'chr6')

    def test_FileNotFound(self):
        with self.assertRaises(FileNotFoundError):
            SnifflesConfig(*self._get_common_args(), '--regions', 'regions.bed')

    def test_FileInvalidLines(self):
        data = """
... <- invalid line
chr1\t100\t200\n  valid line

"""
        with self.assertLogs(level=logging.WARNING), patch("builtins.open", mock_open(read_data=data)) as mock_file:
            config = SnifflesConfig(*self._get_common_args(), '--regions', 'regions.bed')

        self.assertDictEqual(
            config.regions_by_contig,
            {'chr1': [Region('chr1', 100, 200)]}
        )
