# Sniffles
Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs using evidence from split-read alignments, high-mismatch regions, and coverage analysis. Please note the current version of Sniffles requires output from BWA-MEM with the optional SAM attributes enabled! If you experience problems or have suggestions please contact: fritz.sedlazeck@gmail.com

**************************************

INSTALL:

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
```
./sniffles -m reads.bam -v calls.vcf
```

Options:

  ./sniffles  [--maria_format <string>] [--min_reads_phase <int>] [-l
               <int>] [--regions <string>] [--corridor <int>] [-r <string>]
               [-t <int>] [-b <string>] [-d <int>] [--max_num_splits <int>]
               [-s <int>] [-n <int>] [-c <int>] [-v <string>] [-q <int>] -m
               <string> [--use_MD_Cigar] [--re-align] [--] [--version]
               [-h]


Where: 

   --maria_format <string>
     Simplified format of bede Format.

   --min_reads_phase <int>
     Minimum reads overlapping two SV to phase them together. Default: 1

   -l <int>,  --min_length <int>
     Minimum length of SV to be reported. Default:0

   --regions <string>
     List of regions CHR:start-stop; to check

   --corridor <int>
     Maximum size of corridor for realignment. Default: 2000

   -r <string>,  --reference <string>
     Reference fasta sequence for realign step

   -t <int>,  --threads <int>
     Number of threads to use. Default: 3

   -b <string>,  --bede <string>
     Bede output file name

   -d <int>,  --max_distance <int>
     Maximum distance to group SV together. Default: 1kb

   --max_num_splits <int>
     Maximum number of splits per read to be still taken into account.
     Default: 4

   -s <int>,  --min_support <int>
     Minimum number of reads that support a SV. Default: 10

   -n <int>,  --num_reads_report <int>
     Report up to N reads that support the SV. Default: 0

   -c <int>,  --min_cigar_event <int>
     Minimum Cigar Event (e.g. Insertion, deletion) to take into account.
     Default:50 

   -v <string>,  --vcf <string>
     VCF output file name

   -q <int>,  --minmapping_qual <int>
     Minimum Mapping Quality. Default: 20

   -m <string>,  --mapped_reads <string>
     (required)  Bam File

   --use_MD_Cigar
     Enables Sniffles to use the alignemtn information to screen for
     susbicious regions.

   --re-align
     Enables the realignment of reads at predicted SV sites. Leads to more
     accurate breakpoint predictions.

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Sniffles version 0.0.1
