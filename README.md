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
