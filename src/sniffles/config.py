#!/usr/bin/env python3
#
# Sniffles2
# A fast structural variant caller for long-read sequencing data
#
# Created: 18.10.2021
# Author:  Moritz Smolka
# Contact: moritz.g.smolka@gmail.com
#

import os
import sys
import datetime
import argparse

from sniffles import util

VERSION="Sniffles2"
BUILD="2.0.3"
SNF_VERSION="S2_rc3"

class ArgFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass

def tobool(v):
    if v==True or v==False:
        return v
    elif v.lower()=="true" or v=="1":
        return True
    elif v.lower()=="false" or v=="0":
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value (True | False) required for argument")

def from_cmdline():
    header=f"Sniffles2: A fast structural variant (SV) caller for long-read sequencing data\n Version {BUILD}\n Contact: moritz.g.smolka@gmail.com"
    example=""" Usage example A - Call SVs for a single sample:
    sniffles --input sorted_indexed_alignments.bam --vcf output.vcf

    ... OR, with CRAM input and bgzipped+tabix indexed VCF output:
      sniffles --input sample.cram --vcf output.vcf.gz

    ... OR, producing only a SNF file with SV candidates for later multi-sample calling:
      sniffles --input sample1.bam --snf sample1.snf

    ... OR, simultaneously producing a single-sample VCF and SNF file for later multi-sample calling:
      sniffles --input sample1.bam --vcf sample1.vcf.gz --snf sample1.snf

    ... OR, with additional options to specify tandem repeat annotations (for improved call accuracy), reference (for DEL sequences) and non-germline mode for detecting rare SVs:
      sniffles --input sample1.bam --vcf sample1.vcf.gz --tandem-repeats tandem_repeats.bed --reference genome.fa --non-germline

 Usage example B - Multi-sample calling:
    Step 1. Create .snf for each sample: sniffles --input sample1.bam --snf sample1.snf
    Step 2. Combined calling: sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf

    ... OR, using a .tsv file containing a list of .snf files, and custom sample ids in an optional second column (one sample per line):
    Step 2. Combined calling: sniffles --input snf_files_list.tsv --vcf multisample.vcf

 Usage example C - Determine genotypes for a set of known SVs (force calling):
    sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf
    """
    usage="sniffles --input SORTED_INPUT.bam [--vcf OUTPUT.vcf] [--snf MERGEABLE_OUTPUT.snf] [--threads 4] [--non-germline]\n\n" + header + "\n\n" + example + "\n\n Use --help for full parameter/usage information\n \n"
    parser = argparse.ArgumentParser(description="", epilog=example, formatter_class=lambda prog: ArgFormatter(prog,max_help_position=100,width=150), usage=usage)
    parser.add_argument("--version", action="version", version=f"Sniffles2, Version {BUILD}")

    main_args = parser.add_argument_group("Common parameters")
    main_args.add_argument("-i","--input", metavar="IN", type=str, help="For single-sample calling: A coordinate-sorted and indexed .bam/.cram (BAM/CRAM format) file containing aligned reads. - OR - For multi-sample calling: Multiple .snf files (generated before by running Sniffles2 for individual samples with --snf)", required=True, nargs="+")
    main_args.add_argument("-v","--vcf", metavar="OUT.vcf", type=str, help="VCF output filename to write the called and refined SVs to. If the given filename ends with .gz, the VCF file will be automatically bgzipped and a .tbi index built for it.", required=False)
    main_args.add_argument("--snf", metavar="OUT.snf", type=str, help="Sniffles2 file (.snf) output filename to store candidates for later multi-sample calling", required=False)
    main_args.add_argument("--reference", metavar="reference.fasta", type=str, help="(Optional) Reference sequence the reads were aligned against. To enable output of deletion SV sequences, this parameter must be set.", default=None)
    main_args.add_argument("--tandem-repeats", metavar="IN.bed", type=str, help="(Optional) Input .bed file containing tandem repeat annotations for the reference genome.", default=None)
    main_args.add_argument("--non-germline", help="Call non-germline SVs (rare, somatic or mosaic SVs)", default=False, action="store_true")
    main_args.add_argument("--phase", help="Determine phase for SV calls (requires the input alignments to be phased)", default=False, action="store_true")
    main_args.add_argument("-t","--threads", metavar="N", type=int, help="Number of parallel threads to use (speed-up for multi-core CPUs)", default=4)

    filter_args = parser.add_argument_group("SV Filtering parameters")
    filter_args.add_argument("--minsupport", metavar="auto", type=str, help="Minimum number of supporting reads for a SV to be reported (default: automatically choose based on coverage)", default="auto")
    filter_args.add_argument("--minsupport-auto-mult", metavar="0.1/0.025", type=float, help="Coverage based minimum support multiplier for germline/non-germline modes (only for auto minsupport) ", default=None)
    filter_args.add_argument("--minsvlen", metavar="N", type=int, help="Minimum SV length (in bp)", default=35)
    filter_args.add_argument("--minsvlen-screen-ratio", metavar="N", type=float, help="Minimum length for SV candidates (as fraction of --minsvlen)", default=0.95)
    filter_args.add_argument("--mapq", metavar="N", type=int, help="Alignments with mapping quality lower than this value will be ignored", default=25)
    filter_args.add_argument("--no-qc", help="Output all SV candidates, disregarding quality control steps.", default=False, action="store_true")
    filter_args.add_argument("--qc-stdev", help="Apply filtering based on SV start position and length standard deviation", metavar="True", type=tobool, default=True)
    filter_args.add_argument("--qc-stdev-abs-max", help="Maximum standard deviation for SV length and size (in bp)", metavar="N", type=int, default=500)
    filter_args.add_argument("--qc-strand", help="Apply filtering based on strand support of SV calls", metavar="False", type=tobool, default=False)
    filter_args.add_argument("--qc-coverage", help="Minimum surrounding region coverage of SV calls", metavar="N", type=int, default=1)
    filter_args.add_argument("--long-ins-length", help="Insertion SVs longer than this value are considered as hard to detect based on the aligner and read length and subjected to more sensitive filtering.", metavar="2500", type=int, default=2500)
    filter_args.add_argument("--long-del-length", help="Deletion SVs longer than this value are subjected to central coverage drop-based filtering (Not applicable for --non-germline)", metavar="50000", type=int, default=50000)
    filter_args.add_argument("--long-del-coverage", help="Long deletions with central coverage (in relation to upstream/downstream coverage) higher than this value will be filtered (Not applicable for --non-germline)", metavar="0.66", type=float, default=0.66)
    filter_args.add_argument("--long-dup-length", help="Duplication SVs longer than this value are subjected to central coverage increase-based filtering (Not applicable for --non-germline)", metavar="50000", type=int, default=50000)
    filter_args.add_argument("--long-dup-coverage", help="Long duplications with central coverage (in relation to upstream/downstream coverage) lower than this value will be filtered (Not applicable for --non-germline)", metavar="1.33", type=float, default=1.33)
    filter_args.add_argument("--max-splits-kb", metavar="N", type=int, help="Additional number of splits per kilobase read sequence allowed before reads are ignored", default=0.1)
    filter_args.add_argument("--max-splits-base", metavar="N", type=int, help="Base number of splits allowed before reads are ignored (in addition to --max-splits-kb)", default=3)
    filter_args.add_argument("--min-alignment-length", metavar="N", type=int, help="Reads with alignments shorter than this length (in bp) will be ignored", default=1000)
    filter_args.add_argument("--phase-conflict-threshold", metavar="F", type=float, help="Maximum fraction of conflicting reads permitted for SV phase information to be labelled as PASS (only for --phase)", default=0.1)
    filter_args.add_argument("--detect-large-ins", help="Infer insertions that are longer than most reads and therefore are spanned by few alignments only.", metavar="True", type=tobool, default=True)
    #filter_args.add_argument("--large-ins-threshold", metavar="N", type=int, help="Minimum clipping at read ends to be considered a potential large insertion (only with --detect-large-ins)", default=5000)

    cluster_args = parser.add_argument_group("SV Clustering parameters")
    cluster_args.add_argument("--cluster-binsize", metavar="N", type=int, help="Initial screening bin size in bp", default=100)
    cluster_args.add_argument("--cluster-r", metavar="R", type=float, help="Multiplier for SV start position standard deviation criterion in cluster merging", default=2.5)
    cluster_args.add_argument("--cluster-repeat-h", metavar="H", type=float, help="Multiplier for mean SV length criterion for tandem repeat cluster merging", default=1.5)
    cluster_args.add_argument("--cluster-repeat-h-max", metavar="N", type=float, help="Max. merging distance based on SV length criterion for tandem repeat cluster merging", default=1000)
    cluster_args.add_argument("--cluster-merge-pos", metavar="N", type=int, help="Max. merging distance for insertions and deletions on the same read and cluster in non-repeat regions", default=150)
    cluster_args.add_argument("--cluster-merge-len", metavar="F", type=float, help="Max. size difference for merging SVs as fraction of SV length", default=0.33)
    cluster_args.add_argument("--cluster-merge-bnd", metavar="N", type=int, help="Max. merging distance for breakend SV candidates.", default=1500)

    genotype_args = parser.add_argument_group("SV Genotyping parameters")
    genotype_args.add_argument("--genotype-ploidy", metavar="N", type=int, help="Sample ploidy (currently fixed at value 2)", default=2)
    genotype_args.add_argument("--genotype-error", metavar="N", type=float, help="Estimated false positve rate for leads (relating to total coverage)", default=0.05)
    genotype_args.add_argument("--sample-id", type=str, help="Custom ID for this sample, used for later multi-sample calling (stored in .snf)", default=None)
    genotype_args.add_argument("--genotype-vcf", metavar="IN.vcf", type=str, help="Determine the genotypes for all SVs in the given input .vcf file (forced calling). Re-genotyped .vcf will be written to the output file specified with --vcf.", default=None)

    multi_args = parser.add_argument_group("Multi-Sample Calling / Combine parameters")
    multi_args.add_argument("--combine-high-confidence", metavar="F", type=float, help="Minimum fraction of samples in which a SV needs to have individually passed QC for it to be reported in combined output (a value of zero will report all SVs that pass QC in at least one of the input samples)", default=0.0)
    multi_args.add_argument("--combine-low-confidence", metavar="F", type=float, help="Minimum fraction of samples in which a SV needs to be present (failed QC) for it to be reported in combined output", default=0.2)
    multi_args.add_argument("--combine-low-confidence-abs", metavar="N", type=int, help="Minimum absolute number of samples in which a SV needs to be present (failed QC) for it to be reported in combined output", default=3)
    multi_args.add_argument("--combine-null-min-coverage", metavar="N", type=int, help="Minimum coverage for a sample genotype to be reported as 0/0 (sample genotypes with coverage below this threshold at the SV location will be output as ./.)", default=5)
    multi_args.add_argument("--combine-match", metavar="N", type=int, help="Multiplier for maximum deviation of multiple SV's start/end position for them to be combined across samples. Given by max_dev=M*sqrt(min(SV_length_a,SV_length_b)), where M is this parameter.", default=250)
    multi_args.add_argument("--combine-match-max", metavar="N", type=int, help="Upper limit for the maximum deviation computed for --combine-match, in bp.", default=1000)
    multi_args.add_argument("--combine-consensus", help="Output the consensus genotype of all samples", default=False, action="store_true")
    multi_args.add_argument("--combine-separate-intra", help="Disable combination of SVs within the same sample", default=False, action="store_true")
    multi_args.add_argument("--combine-output-filtered", help="Include low-confidence / putative non-germline SVs in multi-calling", default=False, action="store_true")
    #multi_args.add_argument("--combine-exhaustive", help="(DEV) Disable performance optimization in multi-calling", default=False, action="store_true")
    #multi_args.add_argument("--combine-relabel-rare", help="(DEV)", default=False, action="store_true")
    #multi_args.add_argument("--combine-with-missing", help="(DEV)", default=False, action="store_true")

    postprocess_args = parser.add_argument_group("SV Postprocessing, QC and output parameters")
    postprocess_args.add_argument("--output-rnames", help="Output names of all supporting reads for each SV in the RNAMEs info field", default=False, action="store_true")
    postprocess_args.add_argument("--no-consensus", help="Disable consensus sequence generation for insertion SV calls (may improve performance)", default=False, action="store_true")
    postprocess_args.add_argument("--no-sort", help="Do not sort output VCF by genomic coordinates (may slightly improve performance)", default=False, action="store_true")
    postprocess_args.add_argument("--no-progress", help="Disable progress display", default=False, action="store_true")
    postprocess_args.add_argument("--quiet", help="Disable all logging, except errors", default=False, action="store_true")
    postprocess_args.add_argument("--max-del-seq-len", metavar="N", type=int, help="Maximum deletion sequence length to be output. Deletion SVs longer than this value will be written to the output as symbolic SVs.", default=50000)
    postprocess_args.add_argument("--symbolic", help="Output all SVs as symbolic, including insertions and deletions, instead of reporting nucleotide sequences.", default=False, action="store_true")

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
    developer_args.add_argument("--dev-dump-coverage", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--dev-no-resplit", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--dev-skip-snf-validation", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--low-memory", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--repeat", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--qc-nm", default=False, action="store_true", help=argparse.SUPPRESS)
    developer_args.add_argument("--qc-nm-max", metavar="F", type=float, default=0.2, help=argparse.SUPPRESS)
    #developer_args.add_argument("--qc-strand", help="(DEV)", default=False, action="store_true")

    config=parser.parse_args()

    if config.quiet:
        sys.stdout=open(os.devnull,"w")

    config.start_date=datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

    config.sort=not config.no_sort
    #if config.low_memory:
    #    config.task_count_multiplier=64
    #else:
    #    config.task_count_multiplier=1
    config.task_count_multiplier=0

    config.version=VERSION
    config.build=BUILD
    config.snf_format_version=SNF_VERSION
    config.command=" ".join(sys.argv)

    if config.dev_call_region != None:
        region_contig,region_startend=config.dev_call_region.replace(",","").split(":")
        start,end=region_startend.split("-")
        config.dev_call_region=dict(contig=region_contig,start=int(start),end=int(end))

    #"--minsvlen" parameter is for final output filtering
    #for intermediate steps, a lower threshold is used to account for sequencing, mapping imprecision
    config.minsvlen_screen=int(config.minsvlen_screen_ratio*config.minsvlen)
    #config.minsupport_screen=max(1,int(0.333*config.minsupport*(config.cluster_binsize/100.0)))

    if config.minsupport!="auto":
        config.minsupport=int(config.minsupport)

    #--minsupport auto defaults
    config.minsupport_auto_base=1.5
    config.minsupport_auto_regional_coverage_weight=0.75

    if config.minsupport_auto_mult==None:
        if config.non_germline:
            config.minsupport_auto_mult=0.025
        else:
            config.minsupport_auto_mult=0.1

    if config.non_germline:
        config.qc_nm=True

    config.coverage_binsize=config.cluster_binsize
    config.coverage_binsize_combine=config.cluster_binsize*5

    config.coverage_updown_bins=5
    config.coverage_shift_bins=3

    #INS Consensus parameters
    #config.consensus_max_reads=20
    #config.consensus_max_reads_bin=10
    config.consensus_min_reads=4
    config.consensus_kmer_len=6
    config.consensus_kmer_skip_base=3
    config.consensus_kmer_skip_seqlen_mult=1.0/500.0
    config.consensus_low_threshold=0.0 #0.15

    #Large INS
    config.long_ins_rescale_base=1.66
    config.long_ins_rescale_mult=0.33

    #BND
    config.bnd_cluster_length=1000
    config.bnd_cluster_resplit=0

    #Genotyping
    config.genotype_format="GT:GQ:DR:DV"
    config.genotype_none=(".",".",0,0,0)
    config.genotype_null=(0,0,0,0,0)
    config.genotype_min_z_score=5
    if config.genotype_ploidy!=2:
        util.fatal_error("Currently only --genotype-ploidy 2 is supported")

    #SNF
    config.snf_block_size=10**5
    config.snf_combine_keep_open=True #Keep file handles open during .snf combining (might be an issue if the number of .snf files to merge is very large)

    #Combine
    config.combine_exhaustive=False
    config.combine_relabel_rare=False
    config.combine_overlap_abs=2500
    config.combine_min_size=100

    #Misc
    config.precise=25 #Max. sum of pos and length stdev for SVs to be labelled PRECISE
    config.resplit_binsize=20
    config.tandem_repeat_region_pad=500
    config.id_prefix="Sniffles2."

    config.dev_profile=False

    config.workdir=os.getcwd()

    return config
