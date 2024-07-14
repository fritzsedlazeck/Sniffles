#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created:     18.10.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#

import os
import sys
import datetime
import argparse
from collections import defaultdict

from typing import Union, Optional

from sniffles import util
from sniffles.region import Region

VERSION = "Sniffles2"
BUILD = "2.4"
SNF_VERSION = "S2_rc4"


class ArgFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action_invocation(self, action):
        if action.metavar is not None:
            return super()._format_action_invocation(action)
        else:
            parts = []
            if action.option_strings:
                parts.extend(action.option_strings)
            elif action.nargs == 0:
                parts.append(action.dest)
            else:
                metavar, = self._metavar_formatter(action, action.dest)(1)
                parts.append('%s' % metavar)
            return ', '.join(parts)


class ExampleAction(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print("sniffles example commands:")
        print("""
Call SVs for a single sample
=============================
     -> sniffles --input sorted_indexed_alignments.bam --vcf output.vcf

   ... OR, with CRAM input and bgzipped+tabix indexed VCF output:
     -> sniffles --input sample.cram --vcf output.vcf.gz

   ... OR, producing only a SNF file with SV candidates:
     -> sniffles --input sample1.bam --snf sample1.snf

   ... OR, simultaneously produce a single-sample VCF and SNF file:
     -> sniffles --input sample1.bam --vcf sample1.vcf.gz --snf sample1.snf

   ... OR, with tandem repeat annotations, reference (for DEL sequences) and mosaic mode for detecting rare SVs:
     -> sniffles --input sample1.bam --vcf sample1.vcf.gz --tandem-repeats tandem_repeats.bed --reference genome.fa --mosaic

Multi-sample calling
====================
   Step 1. Create .snf for each sample:
     -> sniffles --input sample1.bam --snf sample1.snf
   Step 2. Combined calling:
     -> sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf

   ... OR, using a .tsv file containing a list of .snf files and sample ids (one sample per line):
   Step 2. Combined calling:
     -> sniffles --input snf_files_list.tsv --vcf multisample.vcf

Determine genotypes for a set of known SVs (force calling)
==========================================================
     -> sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf
""")
        parser.exit()

class SnifflesConfig(argparse.Namespace):
    header = f"Sniffles2: A fast structural variant (SV) caller for long-read sequencing data\n Version {BUILD}\n Contact: sniffles@romanek.at"
    usage = "sniffles --input SORTED_INPUT.bam [--vcf OUTPUT.vcf] [--snf MERGEABLE_OUTPUT.snf] [--threads 4] [--mosaic]\n\n" +  \
            header +  \
            "\n\n Use --help for full parameter information\n" + \
            " Use --example for detailed usage information\n\n"

    parser: argparse.ArgumentParser

    quiet: bool

    dev_call_region: Union[str, dict]

    @property
    def sort(self):
        """
        Output is always sorted
        """
        return not self.no_sort

    input: str
    vcf: str
    snf: str
    reference: str
    tandem_repeats: str
    phase: bool
    threads: int
    contig: Optional[str]
    run_id: str

    @property
    def vcf_output_bgz(self) -> Optional[bool]:
        """
        Should the output vcf file be compressed?
        """
        if self.vcf:
            path, ext = os.path.splitext(self.vcf)
            return ext == ".gz" or ext == ".bgz"


    def add_main_args(self, parser):
        main_args = parser.add_argument_group("Common parameters")
        main_args.add_argument("-i", "--input", metavar="IN", type=str, required=True, nargs="+",
                               help=("For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format) "
                                     "file containing aligned reads. - OR - For multi-sample calling: Multiple .snf files "
                                     "(generated before by running Sniffles2 for individual samples with --snf)"))
        main_args.add_argument("-v", "--vcf", metavar="OUT.vcf", type=str, required=False,
                               help=("VCF output filename to write the called and refined SVs to. If the given filename ends "
                                     "with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it."))
        main_args.add_argument("--snf", metavar="OUT.snf", type=str, required=False,
                               help="Sniffles2 file (.snf) output filename to store candidates for later multi-sample calling")
        main_args.add_argument("--reference", metavar="REF.fa", type=str, default=None,
                               help=("(Optional) Reference sequence the reads were aligned against. To enable output of "
                                     "deletion SV sequences, this parameter must be set."))
        main_args.add_argument("--tandem-repeats", metavar="IN.bed", type=str, default=None,
                               help="(Optional) Input .bed file containing tandem repeat annotations for the reference genome.")
        main_args.add_argument("--regions", metavar="REG.bed", type=str,
                               help="(Optional) Only process the specified regions.", default=None)
        main_args.add_argument("-c", "--contig", default=None, type=str, action="append",
                                help="(Optional) Only process the specified contigs. May be given more than once.")
        main_args.add_argument("--phase", action="store_true",
                                help="Determine phase for SV calls (requires the input alignments to be phased)")
        main_args.add_argument("-t", "--threads", type=int, default=4,
                               help="Number of parallel threads to use (%(default)s)")


    minsupport: Union[str, int]
    minsvlen: int
    minsvlen_screen_ratio: float

    def add_filter_args(self, parser):
        filter_args = parser.add_argument_group("SV Filtering parameters")
        filter_args.add_argument("--minsupport", type=str, default="auto",
                                 help="Min number of supporting reads for a SV to be reported (%(default)s)")
        filter_args.add_argument("--minsupport-auto-mult", type=float, default=None,
                                 help="Coverage based auto-minsupport multiplier for germline mode (0.1/0.025)")
        filter_args.add_argument("--minsvlen", type=int, default=50,
                                 help="Min SV length in bp (%(default)s)")
        filter_args.add_argument("--minsvlen-screen-ratio", type=float, default=0.9,
                                 help="Min length for SV candidates as fraction of --minsvlen (%(default)s)")
        filter_args.add_argument("--mapq", type=int, default=argparse.SUPPRESS,
                                 help="Alignments with mapping quality lower than this value will be ignored")
        filter_args.add_argument("--no-qc", "--qc-output-all", action="store_true",
                                 help="Output all SV candidates, disregarding quality control steps")
        # This was previously a --qc-stdev true argument but now a flag
        filter_args.add_argument("--qc-stdev", action="store_false",
                                 help="Apply filtering based on SV start position and length standard deviation")
        filter_args.add_argument("--qc-stdev-abs-max", type=int, default=500,
                                 help="Max standard deviation for SV length and size in bp (%(default)s)")
        # This was previously an argument but now a flag
        filter_args.add_argument("--qc-strand", action="store_true",
                                 help="Apply filtering based on strand support of SV calls")
        filter_args.add_argument("--qc-coverage", type=int, default=1,
                                 help="Min surrounding region coverage of SV calls (%(default)s)")
        # Trimmed the description, wiki should hold more detail
        filter_args.add_argument("--long-ins-length", type=int, default=2500,
                                 help="Insertion SVs longer than this are subjected to more sensitive filtering (%(default)s)")
        filter_args.add_argument("--long-del-length",type=int, default=50000,
                                 help=("Deletion SVs longer than this are subjected to central coverage drop-based filtering. "
                                       "Not applicable for --mosaic (%(default)s)"))
        filter_args.add_argument("--long-inv-length", type=int, default=10000,
                                 help=("Inversion SVs longer than this value are not subjected to central coverage "
                                       "drop-based filtering (%(default)s)"))
        filter_args.add_argument("--long-del-coverage", type=float, default=0.66,
                                 help=("Long deletions with central coverage higher than this value will be filtered. "
                                      "Not applicable for --mosaic (%(default)s)"))
        filter_args.add_argument("--long-dup-length", type=int, default=50000,
                                 help=("Duplication SVs longer than this value are subjected to central coverage increase-based "
                                       "filtering. Not applicable for --mosaic (%(default)s)"))
        # This was previously an argument but now a flag
        filter_args.add_argument("--qc-bnd-filter-strand", action="store_false",
                                 help="Filter breakends that do not have support for both strands")
        filter_args.add_argument("--bnd-min-split-length", type=int, default=1000,
                                 help="Min length of read splits to be considered for breakends (%(default)s)")
        filter_args.add_argument("--long-dup-coverage", type=float, default=1.33,
                                 help=("Long duplications with central coverage lower than this value will be filtered. "
                                       "Not applicable for --mosaic (%(default)s)"))
        filter_args.add_argument("--max-splits-kb", type=float, default=0.1,
                                 help="Additional number of splits per kilobase read sequence allowed before reads are ignored (%(default)s)")
        filter_args.add_argument("--max-splits-base", metavar="N", type=int, default=3,
                                 help="Base number of splits allowed before reads are ignored (%(default)s)")
        filter_args.add_argument("--min-alignment-length", type=int, default=argparse.SUPPRESS,
                                 help="Reads with alignments shorter than this length in bp will be ignored")
        filter_args.add_argument("--phase-conflict-threshold", type=float, default=0.1,
                                 help="Max fraction of conflicting reads permitted for SV phase information to be labelled as PASS. Only for --phase (%(default)s)")
        # This was previously an argument but now a flag
        filter_args.add_argument("--detect-large-ins", action="store_false",
                                 help="Infer insertions that are longer than most reads and therefore are spanned by few alignments only.")
        # filter_args.add_argument("--large-ins-threshold", metavar="N", type=int, help="Min clipping at read ends to be considered a potential large insertion (only with --detect-large-ins)", default=5000)

    cluster_binsize: int
    cluster_binsize_combine_mult: int

    def add_cluster_args(self, parser):
        cluster_args = parser.add_argument_group("SV Clustering parameters")
        cluster_args.add_argument("--cluster-binsize", type=int, default=100,
                                  help="Initial screening bin size in bp (%(default)s)")
        cluster_args.add_argument("--cluster-r", type=float, default=2.5,
                                  help="Multiplier for SV start position standard deviation criterion in cluster merging (%(default)s)")
        cluster_args.add_argument("--cluster-repeat-h", type=float, default=1.5,
                                 help="Multiplier for mean SV length criterion for tandem repeat cluster merging (%(default)s)")
        cluster_args.add_argument("--cluster-repeat-h-max", type=float, default=1000,
                                  help="Max. merging distance based on SV length criterion for tandem repeat cluster merging (%(default)s)")
        cluster_args.add_argument("--cluster-merge-pos", type=int, default=150,
                                 help="Max. merging distance for insertions and deletions on the same read and cluster in non-repeat regions (%(default)s)")
        cluster_args.add_argument("--cluster-merge-len", type=float, default=0.33,
                                 help="Max. size difference for merging SVs as fraction of SV length (%(default)s)")
        cluster_args.add_argument("--cluster-merge-bnd", type=int, default=1000,
                                  help="Max. merging distance for breakend SV candidates (%(default)s)")

    genotype_ploidy: int
    genotype_vcf: str

    def add_genotype_args(self, parser):
        genotype_args = parser.add_argument_group("SV Genotyping parameters")
        genotype_args.add_argument("--genotype-ploidy", type=int, default=2,
                                   help="Sample ploidy (%(default)s)")
        genotype_args.add_argument("--genotype-error", type=float, default=0.05,
                                   help="Estimated false positive rate for leads (%(default)s)")
        # I assume this is also used in the output VCF, so I removed notes about .snf
        genotype_args.add_argument("--sample-id", type=str, default="SAMPLE",
                                  help="Custom ID for this sample (%(default)s))")
        genotype_args.add_argument("--genotype-vcf", metavar="IN.vcf", type=str,
                                   help="Forced calling input.vcf")

    def add_multi_args(self, parser):
        multi_args = parser.add_argument_group("Multi-Sample Calling / Combine parameters")
        multi_args.add_argument("--combine-high-confidence", type=float, default=0.0,
                                help="Min fraction of passed QC samples an SV needs (%(default)s)")
        multi_args.add_argument("--combine-low-confidence", type=float, default=0.2,
                                help="Min fraction of present samples an SV needs (%(default)s)")
        multi_args.add_argument("--combine-low-confidence-abs", type=int, default=2,
                                help="Min number of present samples an SV needs (%(default)s)")
        multi_args.add_argument("--combine-null-min-coverage", type=int, default=5,
                                help="Min coverage for a genotype to be reported as 0/0 instead of ./. (%(default)s)")
        # The formula and detail should be in wiki
        multi_args.add_argument("--combine-match", type=int, default=250,
                                help=("Multiplier for maximum deviation of multiple SV's start/end position for them to "
                                     "be combined across samples. Given by max_dev=M*sqrt(min(SV_length_a,SV_length_b)), "
                                     "where M is this parameter (%(default)s)"))
        multi_args.add_argument("--combine-match-max", type=int, default=1000,
                                 help="Upper limit for the max deviation computed for --combine-match, in bp (%(default)s)")
        multi_args.add_argument("--combine-separate-intra", action="store_true", 
                                help="Disable combination of SVs within the same sample")
        multi_args.add_argument("--combine-output-filtered", action="store_true",
                                 help="Include low-confidence / mosaic SVs in multi-calling")
        multi_args.add_argument("--combine-pair-relabel", action="store_true",
                                help="Override low-quality genotypes when combining paired samples")
        multi_args.add_argument("--combine-pair-relabel-threshold", default=20, type=int,
                                help="Genotype quality minimum before relabeling (%(default)s)")
        multi_args.add_argument("--combine-close-handles", action="store_true",
                                help="Close .SNF file handles after each use to avoid opened files ulimit when merging many samples.")
        multi_args.add_argument("--combine-pctseq", default=0.7, type=float,
                                help="Min alignment distance as percent of SV length to be merged. 0=off (%(default)s)")
        multi_args.add_argument("--combine-max-inmemory-results", default=20, type=int, help=argparse.SUPPRESS)
        # multi_args.add_argument("--combine-exhaustive", help="(DEV) Disable performance optimization in multi-calling", default=False, action="store_true")
        # multi_args.add_argument("--combine-relabel-rare", help="(DEV)", default=False, action="store_true")
        # multi_args.add_argument("--combine-with-missing", help="(DEV)", default=False, action="store_true")

    allow_overwrite: bool

    def add_postprocess_args(self, parser):
        postprocess_args = parser.add_argument_group("Output formatting parameters")
        postprocess_args.add_argument("--output-rnames", action="store_true",
                                      help="Output names supporting reads in INFO/RNAME")
        postprocess_args.add_argument("--no-consensus", action="store_true",
                                      help="Disable consensus sequence generation for insertion SV calls")
        postprocess_args.add_argument("--no-sort", action="store_true",
                                      help="Do not sort output VCF")
        postprocess_args.add_argument("--no-progress", action="store_true",
                                      help="Disable progress display")
        postprocess_args.add_argument("--quiet", action="store_true",
                                      help="Disable any non-error logging")
        postprocess_args.add_argument("--max-del-seq-len", type=int, default=50000,
                                      help="Max deletion sequence length in output before writing as symbolic <DEL> (%(default)s)")
        postprocess_args.add_argument("--symbolic", action="store_true",
                                      help="Output all SVs as symbolic")
        postprocess_args.add_argument("--allow-overwrite", action="store_true",
                                      help="Allow overwriting existing output files")

    mosaic_include_germline: bool
    mosaic_qc_nm: bool
    # TODO some better rules here
    mosaic_min_reads: int = 3
    mosaic_use_strand_thresholds: int = 10

    def add_mosaic_args(self, parser):
        mosaic_args = parser.add_argument_group("Mosaic/somatic calling mode parameters")
        mosaic_args.add_argument("--mosaic", action="store_true",
                                 help="Turn on mosaic calling")
        mosaic_args.add_argument("--mosaic-af-max", type=float, default=0.2,
                                 help="Max allele frequency for which SVs are considered mosaic (%(default)s)")
        mosaic_args.add_argument("--mosaic-af-min", type=float, default=0.05, 
                                 help="Min allele frequency for mosaic SVs to be output (%(default)s)")
        mosaic_args.add_argument("--mosaic-qc-invdup-min-length", type=int,  default=500,
                                 help="Min SV length for mosaic inversion and duplication SVs (%(default)s)")
        # This flag is always True
        mosaic_args.add_argument("--mosaic-qc-nm", default=True, action="store_true", help=argparse.SUPPRESS)
        mosaic_args.add_argument("--mosaic-qc-nm-mult", type=float, default=1.66, help=argparse.SUPPRESS)
        mosaic_args.add_argument("--mosaic-qc-coverage-max-change-frac", type=float, default=0.1,
                                 help="Max relative coverage change across breakpoints (%(default)s)")
        # This parameter is now a flag
        mosaic_args.add_argument("--mosaic-qc-strand", action="store_false", 
                                 help="Apply filtering based on strand support of calls")
        mosaic_args.add_argument("--mosaic-include-germline", action="store_true",
                                 help="Report germline SVs as well in mosaic mode")

    qc_nm: bool
    combine_consensus: bool
    low_memory: bool

    def add_developer_args(self, parser):
        developer_args = parser.add_argument_group("Developer parameters")
        developer_args.add_argument("--dev-cache", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-cache-dir", metavar="PATH", type=str, default=None, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-debug-svtyping", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-keep-lowqual-splits", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-call-region", metavar="REGION", type=str, default=None, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-dump-clusters", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-merge-inline", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-seq-cache-maxlen", metavar="N", type=int, default=50000, help=argparse.SUPPRESS)
        developer_args.add_argument("--consensus-max-reads", metavar="N", type=int, default=20, help=argparse.SUPPRESS)
        developer_args.add_argument("--consensus-max-reads-bin", metavar="N", type=int, default=10, help=argparse.SUPPRESS)
        developer_args.add_argument("--combine-consensus", help="Output the consensus genotype of all samples", default=False, action="store_true")
        developer_args.add_argument("--dev-dump-coverage", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-no-resplit", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-no-resplit-repeat", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-skip-snf-validation", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--low-memory", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--repeat", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--qc-nm", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--qc-nm-mult", metavar="F", type=float, default=1.66, help=argparse.SUPPRESS)
        developer_args.add_argument("--qc-coverage-max-change-frac", help="Max relative coverage change across SV breakpoints", metavar="F", type=float, default=-1)
        developer_args.add_argument("--coverage-updown-bins", metavar="N", type=int, default=5, help=argparse.SUPPRESS)
        developer_args.add_argument("--coverage-shift-bins", metavar="N", type=int, default=3, help=argparse.SUPPRESS)
        developer_args.add_argument("--coverage-shift-bins-min-aln-length", metavar="N", type=int, default=1000, help=argparse.SUPPRESS)
        developer_args.add_argument("--cluster-binsize-combine-mult", metavar="N", type=int, default=5, help=argparse.SUPPRESS)
        developer_args.add_argument("--cluster-resplit-binsize", metavar="N", type=int, default=20, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-trace-read", default=False, metavar="read_id", type=str, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-split-max-query-distance-mult", metavar="N", type=int, default=5, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-disable-interblock-threads", default=False, help=argparse.SUPPRESS, action="store_true")
        developer_args.add_argument("--dev-combine-medians", default=False, help=argparse.SUPPRESS, action="store_true")
        developer_args.add_argument("--dev-monitor-memory", metavar="N", type=int, default=0, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-monitor-filename", metavar="memory.csv", type=str, help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-debug-log", default=False, action="store_true", help=argparse.SUPPRESS)
        developer_args.add_argument("--dev-progress-log", default=False, action="store_true", help=argparse.SUPPRESS)

        # developer_args.add_argument("--qc-strand", help="(DEV)", default=False, action="store_true")

    def __init__(self, *args, **kwargs):
        super().__init__(**kwargs)

        parser = argparse.ArgumentParser(description="", usage=self.usage,
            formatter_class=lambda prog: ArgFormatter(prog, width=100))
        parser.add_argument('--example', action=ExampleAction, help='Show example usage and exit')
        parser.add_argument("--version", action="version", version=f"Sniffles2, Version {BUILD}")

        self.add_main_args(parser)
        self.add_filter_args(parser)
        self.add_cluster_args(parser)
        self.add_genotype_args(parser)
        self.add_multi_args(parser)
        self.add_postprocess_args(parser)
        self.add_mosaic_args(parser)
        self.add_developer_args(parser)

        parser.parse_args(args=args or None, namespace=self)

        if self.quiet:
            sys.stdout = open(os.devnull, "w")

        self.start_date = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
        self.run_id = f'{os.environ.get("SLURM_JOB_ID") or os.getpid()}'

        self.task_count_multiplier = 0

        self.version = VERSION
        self.build = BUILD
        self.snf_format_version = SNF_VERSION
        self.command = " ".join(sys.argv)

        if self.dev_call_region is not None:
            region_contig, region_startend = self.dev_call_region.replace(",", "").split(":")
            start, end = region_startend.split("-")
            self.dev_call_region = dict(contig=region_contig, start=int(start), end=int(end))

        if self.contig and self.regions:
            util.fatal_error('Please provide either --contig or --regions, not both.')

        if self.regions is not None:
            regions = defaultdict(list)
            with open(self.regions, 'r') as f:
                for line in f.readlines():
                    if line.startswith('#') or line.strip() == '':
                        continue
                    else:
                        r = Region.from_bed_line(line)
                        if r is not None:
                            regions[r.contig].append(r)
            self.regions_by_contig = regions
        else:
            self.regions_by_contig = {}

        # "--minsvlen" parameter is for final output filtering
        # for intermediate steps, a lower threshold is used to account for sequencing, mapping imprecision
        self.minsvlen_screen = int(self.minsvlen_screen_ratio * self.minsvlen)
        # config.minsupport_screen=max(1,int(0.333*config.minsupport*(config.cluster_binsize/100.0)))

        if self.minsupport != "auto":
            self.minsupport = int(self.minsupport)

        if not hasattr(self, 'mapq'):
            self.mapq = 0 if self.no_qc else 20
        if not hasattr(self, 'min_alignment_length'):
            self.min_alignment_length = 0 if self.no_qc else 1000

        # --minsupport auto defaults
        self.minsupport_auto_base = 1.5
        self.minsupport_auto_regional_coverage_weight = 0.75

        if self.minsupport_auto_mult is None:
            self.minsupport_auto_mult = 0.1

        self.coverage_binsize = self.cluster_binsize
        self.coverage_binsize_combine = self.cluster_binsize * self.cluster_binsize_combine_mult

        # INS Consensus parameters
        # config.consensus_max_reads=20
        # config.consensus_max_reads_bin=10
        self.consensus_min_reads = 4
        self.consensus_kmer_len = 6
        self.consensus_kmer_skip_base = 3
        self.consensus_kmer_skip_seqlen_mult = 1.0 / 500.0
        self.consensus_low_threshold = 0.0  # 0.15

        # Large INS
        self.long_ins_rescale_base = 1.66
        self.long_ins_rescale_mult = 0.33

        # BND
        self.bnd_cluster_length = 1000

        # Genotyping
        self.genotype_format = "GT:GQ:DR:DV"
        self.genotype_none = (".", ".", 0, 0, 0, (None, None))
        self.genotype_null = (0, 0, 0, 0, 0, (None, None))
        self.genotype_min_z_score = 5
        if self.genotype_ploidy != 2:
            util.fatal_error("Currently only --genotype-ploidy 2 is supported")

        # SNF
        self.snf_block_size = 10 ** 5

        # Combine
        self.combine_exhaustive = False
        self.combine_relabel_rare = False
        self.combine_overlap_abs = 2500
        self.combine_min_size = 100

        # Misc
        self.precise = 25  # Max. sum of pos and length stdev for SVs to be labelled PRECISE
        self.tandem_repeat_region_pad = 500
        self.id_prefix = "Sniffles2."
        self.phase_identifiers = ["1", "2"]

        self.dev_profile = False

        self.workdir = os.getcwd()

        # Mosaic
        if self.mosaic_include_germline:
            self.mosaic = True

        self.qc_nm_measure = self.qc_nm
        if self.mosaic:
            # config.qc_coverage_max_change_frac=config.mosaic_qc_coverage_max_change_frac
            self.qc_nm_measure = self.qc_nm_measure or self.mosaic_qc_nm
            # config.qc_nm_mult=config.mosaic_qc_nm_mult
            # config.qc_strand=config.mosaic_qc_strand
