# Sniffles
Sniffles is a structural variation caller using third generation sequencing (PacBio or Oxford Nanopore). It detects all types of SVs (10bp+) using evidence from split-read alignments, high-mismatch regions, and coverage analysis. Please note the current version of Sniffles requires sorted output from BWA-MEM (use -M and -x parameter), Minimap2 (sam file with Cigar & MD string) or NGMLR. If you experience problems or have suggestions please contact: fritz.sedlazeck@gmail.com


Please see our github wiki for more information (https://github.com/fritzsedlazeck/Sniffles/wiki)


# How to build Sniffles
<pre>wget https://github.com/fritzsedlazeck/Sniffles/archive/master.tar.gz -O Sniffles.tar.gz
tar xzvf Sniffles.tar.gz
cd Sniffles-master/
mkdir -p build/
cd build/
cmake ..
make

cd ../bin/sniffles*
./sniffles</pre>

Note Mac users often have to provide parameters to the cmake command:
<pre>cmake -D CMAKE_C_COMPILER=/opt/local/bin/gcc-mp-4.7 -D CMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.7 .. 
</pre>


**************************************
## NGMLR
Sniffles performs best with the mappings of NGMLR our novel long read mapping method. 
Please see:
https://github.com/philres/ngmlr

****************************************
## Citation:
Please see and cite our paper:
https://www.nature.com/articles/s41592-018-0001-7
  
**************************************
## Poster & Talks:

[Accurate and fast detection of complex and nested structural variations using long read technologies](http://schatzlab.cshl.edu/presentations/2016/2016.10.28.BIODATA.PacBioSV.pdf)
Biological Data Science, Cold Spring Harbor Laboratory, Cold Spring Harbor, NY, 26 - 29.10.2016

[NGMLR: Highly accurate read mapping of third generation sequencing reads for improved structural variation analysis](http://www.cibiv.at/~philipp_/files/gi2016_poster_phr.pdf) 
Genome Informatics 2016, Wellcome Genome Campus Conference Centre, Hinxton, Cambridge, UK, 19.09.-2.09.2016

**************************************
## Datasets used in the mansucript:
We provide the NGMLR aligned reads and the Sniffles calls for the data sets used:  

Arabidopsis trio: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Arabidopsis_trio) . 

Genome in the Bottle trio: 
+ Mappings: [ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/](ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/) . 

+ SV calls: [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/GiaB/)

NA12878: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/NA12878/) .  

SKBR3: 
+ [http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/](http://labshare.cshl.edu/shares/schatzlab/www-data/fsedlaze/Sniffles/Skbr3/) . 


