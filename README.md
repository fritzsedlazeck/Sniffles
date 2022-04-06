# Sniffles2
A fast structural variant caller for long-read sequencing, Sniffles2 accurately detect SVs on germline, somatic and population-level for PacBio and Oxford Nanopore read data.

## Quick Start: Germline SV calling using Sniffles2
To call SVs from long read alignments (PacBio / ONT), you can use:

`sniffles -i mapped_input.bam -v output.vcf`

For improved calling in repetitive regions, Sniffles2 accepts a tandem repeat annotations file using the option `--tandem-repeats annotations.bed`. Sniffles2 compatible tandem repeat annotations for human references can be downloaded from the [annotations/ folder](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations).

(see sniffles --help or below for full usage information).

## Installation
You can install Sniffles2 using pip or conda using:

`pip install sniffles`

or

`conda install sniffles=2.0 `

If you previously installed Sniffles1 using conda and want to upgrade to Sniffles2, you can use:

`conda update sniffles=2.0`

## Requirements
* Python >= 3.7
* pysam

#### Tested on:
* python==3.9.5
* pysam==0.16.0.1

## Citation
Please cite our paper at:

https://www.nature.com/articles/s41592-018-0001-7

A new preprint for the new methods and improvements introduced with Sniffles2 is here: https://www.biorxiv.org/content/10.1101/2022.04.04.487055v1 

## Use-Cases / Modes

### A. General (all Modes)
* To output deletion (DEL SV) sequences, the reference genome (.fasta) must be specified using e.g. `--reference reference.fasta`.
* Sniffles2 supports optionally specifying tandem repeat region annotations (.bed), which can improve calling in these regions `--tandem-repeats annotations.bed`. Sniffles2 compatible tandem repeat annotations for human references can be found in the [annotations/ folder](https://github.com/fritzsedlazeck/Sniffles/tree/master/annotations).
* Sniffles2 is fully parallelized and uses 4 threads by default. This value can be adapted using e.g. `--threads 4` as option. Memory requirements will increase with the number of threads used.
* To output read names in SNF and VCF files, the `--output-rnames` option is required.

### B. Multi-Sample SV Calling (Trios, Populations)
Multi-sample SV calling using Sniffles2 population mode works in two steps:

1. Call SV candidates and create an associated .snf file for each sample: `sniffles --input sample1.bam --snf sample1.snf`
2. Combined calling using multiple .snf files into a single .vcf: `sniffles --input sample1.snf sample2.snf ... sampleN.snf --vcf multisample.vcf`

Alternatively, for step 2. you can supply a .tsv file, containing a list of .snf files, and custom sample ids in an optional second column (one sample per line), .e.g.:
2. Combined calling using a .tsv as sample list: `sniffles --input snf_files_list.tsv --vcf multisample.vcf`

### C. Non-Germline SV Calling (Somatic)
To call non-germline SVs (i.e. somatic/mosaic) SVs, the *--non-germline* option should be added, i.e.:

`sniffles --input mapped_input.bam --vcf output.vcf --non-germline`

### D. Genotyping a known set of SVs (Force Calling)
Example command, to determine the genotype of each SV in *input_known_svs.vcf* for *sample.bam* and write the re-genotyped SVs to *output_genotypes.vcf*:

`sniffles --input sample.bam --genotype-vcf input_known_svs.vcf --vcf output_genotypes.vcf`

## Quick Tips

### Input / Output
* .bam or .cram files containing long read alignments (i.e. from minimap2 or ngmlr) are supported as input
* .vcf.gz (bgzipped+tabix indexed) output is supported
* Simultaneous output of both .vcf and .snf file (for multi-sample calling) is supported
