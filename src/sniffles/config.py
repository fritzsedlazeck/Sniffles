#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 18.10.2021
# Author:      Moritz Smolka
# Maintainer:  Hermann Romanek
# Contact:     sniffles@romanek.at
#

import os
import sys
import datetime
import argparse

from typing import Union, Optional

from sniffles import util
from sniffles.region import Region

VERSION = "Sniffles2"
BUILD = "2.4.1"
SNF_VERSION = "S2_rc4"


class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


def tobool(v):
    if v is True or v is False:
        return v
    elif v.strip().lower() == "true" or v.strip() == "1":
        return True
    elif v.strip().lower() == "false" or v.strip() == "0":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value (True | False) required for argument")


class SnifflesConfig(argparse.Namespace):
    header = f"Sniffles2: A fast structural variant (SV) caller for long-read sequencing data\n Version {BUILD}\n Contact: sniffles@romanek.at"
    example = """ Usage example A - Call SVs for a single sample:
       sniffles --input sorted_indexed_alignments.bam --vcf output.vcf

       ... OR, with CRAM input and bgzipped+tabix indexed VCF output:
         sniffles --input sample.cram --vcf output.vcf.gz

       ... OR, producing only a SNF file with SV candidates for later multi-sample calling:
         sniffles --input sample1.bam --snf sample1.snf

       ... OR, simultaneously producing a single-sample VCF and SNF file for later multi-sample calling:
         sniffles --input sample1.bam --vcf sample1.vcf.gz --snf sample1.snf

       ... OR, with additional options to specify tandem repeat annotations (for improved call accuracy), reference (for DEL sequences) and mosaic mode for detecting rare SVs:
         sniffles --input sample1.bam --vcf sample1.vcf.gz --tandem-repeats tandem_repeats.bed --reference genome.fa --mosaic

    Usage example B - Multi-sample calling:
       Step 1. Create .snf for each sample: sniffles --input sample1.bam --snf sample1.snf
       Step 2. Combined calling: sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf

       ... OR, using a .tsv file containing a list of .snf files, and custom sample ids in an optional second column (one sample per line):
       Step 2. Combined calling: sniffles --input snf_files_list.tsv --vcf multisample.vcf

    Usage example C - Determine genotypes for a set of known SVs (force calling):
       sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf
       """
    usage = "sniffles --input SORTED_INPUT.bam [--vcf OUTPUT.vcf] [--snf MERGEABLE_OUTPUT.snf] [--threads 4] [--mosaic]\n\n" + header + "\n\n" + example + "\n\n Use --help for full parameter/usage information\n \n"

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
        main_args.add_argument("-i", "--input", metavar="IN", type=str, help="For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format) file containing aligned reads. - OR - For multi-sample calling: Multiple .snf files (generated before by running Sniffles2 for individual samples with --snf)", required=True, nargs="+")
        main_args.add_argument("-v", "--vcf", metavar="OUT.vcf", type=str, help="VCF output filename to write the called and refined SVs to. If the given filename ends with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.", required=False)
        main_args.add_argument("--snf", metavar="OUT.snf", type=str, help="Sniffles2 file (.snf) output filename to store candidates for later multi-sample calling", required=False)
        main_args.add_argument("--reference", metavar="reference.fasta", type=str, help="(Optional) Reference sequence the reads were aligned against. To enable output of deletion SV sequences, this parameter must be set.", default=None)
        main_args.add_argument("--tandem-repeats", metavar="IN.bed", type=str, help="(Optional) Input .bed file containing tandem repeat annotations for the reference genome.", default=None)
        main_args.add_argument("--phase", help="Determine phase for SV calls (requires the input alignments to be phased)", default=False, action="store_true")
        main_args.add_argument("-t", "--threads", metavar="N", type=int, help="Number of parallel threads to use (speed-up for multi-core CPUs)", default=4)
        main_args.add_argument("-c", "--contig", default=None, type=str, help="(Optional) Only process the specified contigs. May be given more than once.", action="append")
        main_args.add_argument("--regions", metavar="REGIONS.bed", type=str, help="(Optional) Only process the specified regions.", default=None)

    minsupport: Union[str, int]
    minsvlen: int
    minsvlen_screen_ratio: float

    def add_filter_args(self, parser):
        filter_args = parser.add_argument_group("SV Filtering parameters")
        filter_args.add_argument("--minsupport", metavar="auto", type=str, help="Minimum number of supporting reads for a SV to be reported (default: automatically choose based on coverage)", default="auto")
        filter_args.add_argument("--minsupport-auto-mult", metavar="0.1/0.025", type=float, help="Coverage based minimum support multiplier for germline mode (only for auto minsupport) ", default=None)
        filter_args.add_argument("--minsvlen", metavar="N", type=int, help="Minimum SV length (in bp)", default=50)
        filter_args.add_argument("--minsvlen-screen-ratio", metavar="N", type=float, help="Minimum length for SV candidates (as fraction of --minsvlen)", default=0.9)
        filter_args.add_argument("--mapq", metavar="N", type=int, help="Alignments with mapping quality lower than this value will be ignored", default=argparse.SUPPRESS)
        filter_args.add_argument("--no-qc", "--qc-output-all", help="Output all SV candidates, disregarding quality control steps.", default=False, action="store_true")
        filter_args.add_argument("--qc-stdev", help="Apply filtering based on SV start position and length standard deviation", metavar="True", type=tobool, default=True)
        filter_args.add_argument("--qc-stdev-abs-max", help="Maximum standard deviation for SV length and size (in bp)", metavar="N", type=int, default=500)
        filter_args.add_argument("--qc-strand", help="Apply filtering based on strand support of SV calls", metavar="False", type=tobool, default=False)
        filter_args.add_argument("--qc-coverage", help="Minimum surrounding region coverage of SV calls", metavar="N", type=int, default=1)
        filter_args.add_argument("--long-ins-length", help="Insertion SVs longer than this value are considered as hard to detect based on the aligner and read length and subjected to more sensitive filtering.", metavar="2500", type=int, default=2500)
        filter_args.add_argument("--long-del-length", help="Deletion SVs longer than this value are subjected to central coverage drop-based filtering (Not applicable for --mosaic)", metavar="50000", type=int, default=50000)
        filter_args.add_argument("--long-inv-length", help="Inversion SVs longer than this value are not subjected to central coverage drop-based filtering", metavar="10000", type=int, default=10000)
        filter_args.add_argument("--long-del-coverage", help="Long deletions with central coverage (in relation to upstream/downstream coverage) higher than this value will be filtered (Not applicable for --mosaic)", metavar="0.66", type=float, default=0.66)
        filter_args.add_argument("--long-dup-length", help="Duplication SVs longer than this value are subjected to central coverage increase-based filtering (Not applicable for --mosaic)", metavar="50000", type=int, default=50000)
        filter_args.add_argument("--qc-bnd-filter-strand", help="Filter breakends that do not have support for both strands", type=tobool, default=True)
        filter_args.add_argument("--bnd-min-split-length", help="Minimum length of read splits to be considered for breakends", type=int, default=1000)
        filter_args.add_argument("--long-dup-coverage", help="Long duplications with central coverage (in relation to upstream/downstream coverage) lower than this value will be filtered (Not applicable for --mosaic)", metavar="1.33", type=float, default=1.33)
        filter_args.add_argument("--max-splits-kb", metavar="N", type=float, help="Additional number of splits per kilobase read sequence allowed before reads are ignored", default=0.1)
        filter_args.add_argument("--max-splits-base", metavar="N", type=int, help="Base number of splits allowed before reads are ignored (in addition to --max-splits-kb)", default=3)
        filter_args.add_argument("--min-alignment-length", metavar="N", type=int, help="Reads with alignments shorter than this length (in bp) will be ignored", default=argparse.SUPPRESS)
        filter_args.add_argument("--phase-conflict-threshold", metavar="F", type=float, help="Maximum fraction of conflicting reads permitted for SV phase information to be labelled as PASS (only for --phase)", default=0.1)
        filter_args.add_argument("--detect-large-ins", help="Infer insertions that are longer than most reads and therefore are spanned by few alignments only.", metavar="True", type=tobool, default=True)
        # filter_args.add_argument("--large-ins-threshold", metavar="N", type=int, help="Minimum clipping at read ends to be considered a potential large insertion (only with --detect-large-ins)", default=5000)

    cluster_binsize: int
    cluster_binsize_combine_mult: int

    def add_cluster_args(self, parser):
        cluster_args = parser.add_argument_group("SV Clustering parameters")
        cluster_args.add_argument("--cluster-binsize", metavar="N", type=int, help="Initial screening bin size in bp", default=100)
        cluster_args.add_argument("--cluster-r", metavar="R", type=float, help="Multiplier for SV start position standard deviation criterion in cluster merging", default=2.5)
        cluster_args.add_argument("--cluster-repeat-h", metavar="H", type=float, help="Multiplier for mean SV length criterion for tandem repeat cluster merging", default=1.5)
        cluster_args.add_argument("--cluster-repeat-h-max", metavar="N", type=float, help="Max. merging distance based on SV length criterion for tandem repeat cluster merging", default=1000)
        cluster_args.add_argument("--cluster-merge-pos", metavar="N", type=int, help="Max. merging distance for insertions and deletions on the same read and cluster in non-repeat regions", default=150)
        cluster_args.add_argument("--cluster-merge-len", metavar="F", type=float, help="Max. size difference for merging SVs as fraction of SV length", default=0.33)
        cluster_args.add_argument("--cluster-merge-bnd", metavar="N", type=int, help="Max. merging distance for breakend SV candidates.", default=1000)

    genotype_ploidy: int
    genotype_vcf: str

    def add_genotype_args(self, parser):
        genotype_args = parser.add_argument_group("SV Genotyping parameters")
        genotype_args.add_argument("--genotype-ploidy", metavar="N", type=int, help="Sample ploidy (currently fixed at value 2)", default=2)
        genotype_args.add_argument("--genotype-error", metavar="N", type=float, help="Estimated false positive rate for leads (relating to total coverage)", default=0.05)
        genotype_args.add_argument("--sample-id", type=str, help="Custom ID for this sample, used for later multi-sample calling (stored in .snf)", default=None)
        genotype_args.add_argument("--genotype-vcf", metavar="IN.vcf", type=str, help="Determine the genotypes for all SVs in the given input .vcf file (forced calling). Re-genotyped .vcf will be written to the output file specified with --vcf.", default=None)

    def add_multi_args(self, parser):
        multi_args = parser.add_argument_group("Multi-Sample Calling / Combine parameters")
        multi_args.add_argument("--combine-high-confidence", metavar="F", type=float, help="Minimum fraction of samples in which a SV needs to have individually passed QC for it to be reported in combined output (a value of zero will report all SVs that pass QC in at least one of the input samples)", default=0.0)
        multi_args.add_argument("--combine-low-confidence", metavar="F", type=float, help="Minimum fraction of samples in which a SV needs to be present (failed QC) for it to be reported in combined output", default=0.2)
        multi_args.add_argument("--combine-low-confidence-abs", metavar="N", type=int, help="Minimum absolute number of samples in which a SV needs to be present (failed QC) for it to be reported in combined output", default=2)
        multi_args.add_argument("--combine-null-min-coverage", metavar="N", type=int, help="Minimum coverage for a sample genotype to be reported as 0/0 (sample genotypes with coverage below this threshold at the SV location will be output as ./.)", default=5)
        multi_args.add_argument("--combine-match", metavar="N", type=int, help="Multiplier for maximum deviation of multiple SV's start/end position for them to be combined across samples. Given by max_dev=M*sqrt(min(SV_length_a,SV_length_b)), where M is this parameter.", default=250)
        multi_args.add_argument("--combine-match-max", metavar="N", type=int, help="Upper limit for the maximum deviation computed for --combine-match, in bp.", default=1000)
        multi_args.add_argument("--combine-separate-intra", help="Disable combination of SVs within the same sample", default=False, action="store_true")
        multi_args.add_argument("--combine-output-filtered", help="Include low-confidence / mosaic SVs in multi-calling", default=False, action="store_true")
        multi_args.add_argument("--combine-pair-relabel", help="Override low-quality genotypes when combining 2 samples (may be used for e.g. tumor-normal comparisons)", default=False, action="store_true")
        multi_args.add_argument("--combine-pair-relabel-threshold", help="Genotype quality below which a genotype call will be relabeled", default=20, type=int)
        multi_args.add_argument("--combine-close-handles", help="Close .SNF file handles after each use. May lower performance, but may be required when maximum number of file handles supported by OS is reached when merging many samples.", default=False, action="store_true")
        multi_args.add_argument("--combine-pctseq", default=0.7, type=float, help="Minimum alignment distance as percent of SV length to be merged. Set to 0 to disable alignments for merging.")
        multi_args.add_argument("--combine-max-inmemory-results", default=20, type=int, help=argparse.SUPPRESS)
        # multi_args.add_argument("--combine-exhaustive", help="(DEV) Disable performance optimization in multi-calling", default=False, action="store_true")
        # multi_args.add_argument("--combine-relabel-rare", help="(DEV)", default=False, action="store_true")
        # multi_args.add_argument("--combine-with-missing", help="(DEV)", default=False, action="store_true")

    allow_overwrite: bool

    def add_postprocess_args(self, parser):
        postprocess_args = parser.add_argument_group("SV Postprocessing, QC and output parameters")
        postprocess_args.add_argument("--output-rnames", help="Output names of all supporting reads for each SV in the RNAMEs info field", default=False, action="store_true")
        postprocess_args.add_argument("--no-consensus", help="Disable consensus sequence generation for insertion SV calls (may improve performance)", default=False, action="store_true")
        postprocess_args.add_argument("--no-sort", help="Do not sort output VCF by genomic coordinates (may slightly improve performance)", default=False, action="store_true")
        postprocess_args.add_argument("--no-progress", help="Disable progress display", default=False, action="store_true")
        postprocess_args.add_argument("--quiet", help="Disable all logging, except errors", default=False, action="store_true")
        postprocess_args.add_argument("--max-del-seq-len", metavar="N", type=int, help="Maximum deletion sequence length to be output. Deletion SVs longer than this value will be written to the output as symbolic SVs.", default=50000)
        postprocess_args.add_argument("--symbolic", help="Output all SVs as symbolic, including insertions and deletions, instead of reporting nucleotide sequences.", default=False, action="store_true")
        postprocess_args.add_argument("--allow-overwrite", help="Allow overwriting output files if already existing", default=False, action="store_true")

    mosaic_include_germline: bool
    mosaic_qc_nm: bool
    # TODO some better rules here
    mosaic_min_reads: int = 3
    mosaic_use_strand_thresholds: int = 10

    def add_mosaic_args(self, parser):
        mosaic_args = parser.add_argument_group("Mosaic calling mode parameters")
        mosaic_args.add_argument("--mosaic", help="Set Sniffles run mode to detect rare, somatic and mosaic SVs", default=False, action="store_true")
        mosaic_args.add_argument("--mosaic-af-max", help="Maximum allele frequency for which SVs are considered mosaic", metavar="F", default=0.2, type=float)
        mosaic_args.add_argument("--mosaic-af-min", help="Minimum allele frequency for mosaic SVs to be output", metavar="F", default=0.05, type=float)
        mosaic_args.add_argument("--mosaic-qc-invdup-min-length", help="Minimum SV length for mosaic inversion and duplication SVs", metavar="N", default=500, type=int)
        mosaic_args.add_argument("--mosaic-qc-nm", default=True, action="store_true", help=argparse.SUPPRESS)
        mosaic_args.add_argument("--mosaic-qc-nm-mult", metavar="F", type=float, default=1.66, help=argparse.SUPPRESS)
        mosaic_args.add_argument("--mosaic-qc-coverage-max-change-frac", help="Maximum relative coverage change across SV breakpoints", metavar="F", type=float, default=0.1)
        mosaic_args.add_argument("--mosaic-qc-strand", help="Apply filtering based on strand support of SV calls", metavar="True", type=tobool, default=True)
        mosaic_args.add_argument("--mosaic-include-germline", help="Report germline SVs as well in mosaic mode", default=False, action="store_true")

    qc_nm: bool
    combine_consensus: bool
    low_memory: bool

    def add_developer_args(self, parser):
        developer_args = parser.add_argument_group("Developer parameters")

        developer_args.add_argument("--dev-emit-sv-lengths", default=False, action="store_true", help=argparse.SUPPRESS)
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
        developer_args.add_argument("--qc-coverage-max-change-frac", help="Maximum relative coverage change across SV breakpoints", metavar="F", type=float, default=-1)
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

        parser = argparse.ArgumentParser(description="", epilog=self.example, formatter_class=lambda prog: ArgFormatter(prog, max_help_position=100, width=150), usage=self.usage)
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
