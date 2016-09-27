# Sniffles
Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs using evidence from split-read alignments, high-mismatch regions, and coverage analysis. Please note the current version of Sniffles requires output from BWA-MEM! If you experience problems or have suggestions please contact: fritz.sedlazeck@gmail.com

**************************************

INSTALL for Linux:

Download Sniffles:
```
git clone https://github.com/fritzsedlazeck/Sniffles
```

  cd Sniffles
  
  mkdir build
  
  cd build
  
  cmake ..
  
  make
 
  cd ../bin/
  
**************************************

USAGE:

Please use BWA-MEM to align the reads using the parameters: “–M  –x pacbio” . Reads do not need to be error corrected.

Sniffles requires a sorted bam file. 
 

```
./sniffles -m reads.bam -v calls.vcf
```

Options:

     ./sniffles  -m <string> [-s <int>] [--max_num_splits <int>] [-q <int>]
               [-l <int>] [-v <string>] [--bede <string>] [-c <int>] [-t
               <int>] [-d <int>] [-n <int>] [--use_MD_Cigar] [--]
               [--version] [-h]


Where: 

   -m <string>,  --mapped_reads <string>
     (required)  Bam File

   -s <int>,  --min_support <int>
     Minimum number of reads that support a SV. Default: 10

   --max_num_splits <int>
     Maximum number of splits per read to be still taken into account. 
     Default: 4

   -q <int>,  --minmapping_qual <int>
     Minimum Mapping Quality. Default: 20

   -l <int>,  --min_length <int>
     Minimum length of SV to be reported. Default:0

   -v <string>,  --vcf <string>
     VCF output file name

   --bede <string>
     Simplified format of bedpe Format.

   -c <int>,  --min_cigar_event <int>
     Minimum Cigar Event (e.g. Insertion, deletion) to take into account.
     Default:50 

   -t <int>,  --threads <int>
     Number of threads to use. Default: 3

   -d <int>,  --max_distance <int>
     Maximum distance to group SV together. Default: 1kb

   -n <int>,  --num_reads_report <int>
     Report up to N reads that support the SV. Default: 0

   --use_MD_Cigar
     Enables Sniffles to use the alignment information to screen for
     suspicious regions.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Sniffles version 0.0.1
   
   
   
   *****************************************************
   
   
# Parameter explanation:

   
   ```
   -d <int>,  --max_distance <int> :
   ```
   This parameter determines when two SVs calls (e.g. based on two reads) are considered the same. Sniffles assesses the breakpoints and determines if the leftmost breakpoint of both calls are within -d and if the right most breakpoints are within -d distance. If so the two calls are merged and the SVs is then supported by e.g. 2 reads. This parameter is by default 1kbp. However, if you have a high heterogeneity or high polyploidy you might expect similar events in close proximity. Then lowering the parameter might give you that. 
   ```
   -c <int>,  --min_cigar_event <int> :
   ```
   A cigar or MD string are ways the mapper report the alignment in a SAM format (see SAM format for more details). These events can also be noisy. The --min_cigar_event is to control from when on Sniffles detects e.g. deletions or insertions. It is the number of consecutive indel events to be observed before this region becomes a potential SVs. If you trust your aligner you can lower this parameter to also detect very small events. 
  
  ```
 --max_num_splits <int> :
 ```
 I often observe that at very noisy regions e.g. bwa mem tends to split reads. These reads are often falsely aligned. The parameter helps to reduce the falsely called region. Usually a read is split 0-1 times. You can have more complicated SVs were a  read is split more often. If you fully trust your mapper then you could increase this parameter, but be aware that the potential of false calls increases as well.


