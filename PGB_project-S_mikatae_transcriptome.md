# PGB project - S. mikatae transcriptome analysis

## 1. Quality analysis
Before performing the quality analysis we are going to have a look and count the reads.

Output when using `head` on the first set of reads (s_mikatae_read1.fastq):

```
@HWI-D00733:25:C6KC4ANXX:8:1101:1499:2212 1:N:0:GTGAAA
TGGGACGTTCTTGTGTCTCTTTTCGTATCTGTTGTACTTAGGAATGTAAT
+
B@BBBGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGCGGGGGG
@HWI-D00733:25:C6KC4ANXX:8:1101:1988:2212 1:N:0:GTGAAA
TGTCAATCAAGTTATAAAGCCTTCTTGGCAGCATCAGCAGGAGAAACCTT
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGG
@HWI-D00733:25:C6KC4ANXX:8:1101:1762:2232 1:N:0:GTGAAA
GTAAGAAGCAGCTTGAGTAGAAGAACCAGAACCTAAGGAATGTTCTGTGT
```


```bash
%%bash

#Counting how many reads do we have in the first fastq.
READS1="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read1.fastq"

grep "^@" $READS1 | wc -l

```

     26165196



```bash
%%bash

#Counting how many reads do we have in the second fastq.
READS2="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read2.fastq"

grep "^@" $READS2 | wc -l

```

     26759790


Now we are going to use `fastqc` to analyse the quality of the reads.


```bash
%%bash


#First set of reads:
READS1="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read1.fastq"

"/Users/martahuertas/Desktop/Software/FastQC/"fastqc $READS1

#Second set of reads:
READS2="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read2.fastq"

"/Users/martahuertas/Desktop/Software/FastQC/"fastqc $READS2

```

    Analysis complete for s_mikatae_read2.fastq


    Started analysis of s_mikatae_read2.fastq
    Approx 5% complete for s_mikatae_read2.fastq
    Approx 10% complete for s_mikatae_read2.fastq
    Approx 15% complete for s_mikatae_read2.fastq
    Approx 20% complete for s_mikatae_read2.fastq
    Approx 25% complete for s_mikatae_read2.fastq
    Approx 30% complete for s_mikatae_read2.fastq
    Approx 35% complete for s_mikatae_read2.fastq
    Approx 40% complete for s_mikatae_read2.fastq
    Approx 45% complete for s_mikatae_read2.fastq
    Approx 50% complete for s_mikatae_read2.fastq
    Approx 55% complete for s_mikatae_read2.fastq
    Approx 60% complete for s_mikatae_read2.fastq
    Approx 65% complete for s_mikatae_read2.fastq
    Approx 70% complete for s_mikatae_read2.fastq
    Approx 75% complete for s_mikatae_read2.fastq
    Approx 80% complete for s_mikatae_read2.fastq
    Approx 85% complete for s_mikatae_read2.fastq
    Approx 90% complete for s_mikatae_read2.fastq
    Approx 95% complete for s_mikatae_read2.fastq


The results are similar in both cases and look like:

![s_mikatae_read1_quality_fastqc_results.png](attachment:s_mikatae_read1_quality_fastqc_results.png)

### 1.1. Remove adapters to improve quality

We are going to run `Trimmomatic` on our reads to get rid of the adapters that are giving poor quality to the first bases of the reads.


```bash
%%bash

TRIMMOMATIC="/Users/martahuertas/Desktop/Software/Trimmomatic-0.36"
READS1="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read1.fastq"
READS2="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read2.fastq"


java -jar $TRIMMOMATIC"/"trimmomatic-0.36.jar PE -phred33 $READS1 $READS2 s_mikatae_reads1_paired.fastq s_mikatae_reads1_unpaired.fastq s_mikatae_reads2_paired.fastq s_mikatae_reads2_unpaired.fastq ILLUMINACLIP:illumina_truseq_adapters.fa:2:33:20:2:true LEADING:36 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:35
```

    TrimmomaticPE: Started with arguments:
     -phred33 /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read1.fastq /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Raw/s_mikatae_read2.fastq s_mikatae_reads1_paired.fastq s_mikatae_reads1_unpaired.fastq s_mikatae_reads2_paired.fastq s_mikatae_reads2_unpaired.fastq ILLUMINACLIP:illumina_truseq_adapters.fa:2:33:20:2:true LEADING:36 TRAILING:32 SLIDINGWINDOW:4:30 MINLEN:35
    Multiple cores found: Using 4 threads
    java.io.FileNotFoundException: /Users/martahuertas/Desktop/PGB/Project/illumina_truseq_adapters.fa (No such file or directory)
    	at java.io.FileInputStream.open0(Native Method)
    	at java.io.FileInputStream.open(FileInputStream.java:195)
    	at java.io.FileInputStream.<init>(FileInputStream.java:138)
    	at org.usadellab.trimmomatic.fasta.FastaParser.parse(FastaParser.java:54)
    	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.loadSequences(IlluminaClippingTrimmer.java:110)
    	at org.usadellab.trimmomatic.trim.IlluminaClippingTrimmer.makeIlluminaClippingTrimmer(IlluminaClippingTrimmer.java:71)
    	at org.usadellab.trimmomatic.trim.TrimmerFactory.makeTrimmer(TrimmerFactory.java:32)
    	at org.usadellab.trimmomatic.Trimmomatic.createTrimmers(Trimmomatic.java:59)
    	at org.usadellab.trimmomatic.TrimmomaticPE.run(TrimmomaticPE.java:536)
    	at org.usadellab.trimmomatic.Trimmomatic.main(Trimmomatic.java:80)
    Input Read Pairs: 25791495 Both Surviving: 18347000 (71,14%) Forward Only Surviving: 2395887 (9,29%) Reverse Only Surviving: 1599527 (6,20%) Dropped: 3449081 (13,37%)
    TrimmomaticPE: Completed successfully


To know how many reads we have discarded:


```bash
%%bash
#First fastq file

READS1_clean="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads1_paired.fastq"

grep "^@" $READS1_clean | wc -l
```

     18347000



```bash
%%bash
#Second fastq file

READS2_clean="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads2_paired.fastq"

grep "^@" $READS2_clean | wc -l
```

     18347000


As we can see, we have only 71% of the previous reads remaining.

To verify that the read quality is now the desired, que can use `FastQC`again in the set of paired reads:


```bash
%%bash

#First set of reads:
READS1="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads1_paired.fastq"

"/Users/martahuertas/Desktop/Software/FastQC/"fastqc $READS1

#Second set of reads:
READS2="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads2_paired.fastq"

"/Users/martahuertas/Desktop/Software/FastQC/"fastqc $READS2
```

    Analysis complete for s_mikatae_reads2_paired.fastq


    Started analysis of s_mikatae_reads2_paired.fastq
    Approx 5% complete for s_mikatae_reads2_paired.fastq
    Approx 10% complete for s_mikatae_reads2_paired.fastq
    Approx 15% complete for s_mikatae_reads2_paired.fastq
    Approx 20% complete for s_mikatae_reads2_paired.fastq
    Approx 25% complete for s_mikatae_reads2_paired.fastq
    Approx 30% complete for s_mikatae_reads2_paired.fastq
    Approx 35% complete for s_mikatae_reads2_paired.fastq
    Approx 40% complete for s_mikatae_reads2_paired.fastq
    Approx 45% complete for s_mikatae_reads2_paired.fastq
    Approx 50% complete for s_mikatae_reads2_paired.fastq
    Approx 55% complete for s_mikatae_reads2_paired.fastq
    Approx 60% complete for s_mikatae_reads2_paired.fastq
    Approx 65% complete for s_mikatae_reads2_paired.fastq
    Approx 70% complete for s_mikatae_reads2_paired.fastq
    Approx 75% complete for s_mikatae_reads2_paired.fastq
    Approx 80% complete for s_mikatae_reads2_paired.fastq
    Approx 85% complete for s_mikatae_reads2_paired.fastq
    Approx 90% complete for s_mikatae_reads2_paired.fastq
    Approx 95% complete for s_mikatae_reads2_paired.fastq
    Approx 100% complete for s_mikatae_reads2_paired.fastq


Both results look similar and we can see that the poor quality bases have been eliminated. 

![s_mikatae_reads1_paired_fastqc.png](attachment:s_mikatae_reads1_paired_fastqc.png)

## 2. Align RNA-seq reads to genome

First we have to index the genome:


```bash
%%bash 

HISAT2="/Users/martahuertas/Desktop/Software/hisat2-2.2.1"
GENOME="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/s_mikatae.fa"
OUTPUT="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index"

$HISAT2"/"hisat2-build $GENOME $OUTPUT"/"s_mikatae_index
```

    Building DifferenceCoverSample
      Building sPrime
      Building sPrimeOrder
      V-Sorting samples
      V-Sorting samples time: 00:00:00
      Allocating rank array
      Ranking v-sort output
      Ranking v-sort output time: 00:00:00
      Invoking Larsson-Sadakane on ranks
      Invoking Larsson-Sadakane on ranks time: 00:00:00
      Sanity-checking and returning
    Building samples
    Reserving space for 12 sample suffixes
    Generating random suffixes
    QSorting 12 sample offsets, eliminating duplicates
    QSorting sample offsets, eliminating duplicates time: 00:00:00
    Multikey QSorting 12 samples
      (Using difference cover)
      Multikey QSorting samples time: 00:00:00
    Calculating bucket sizes
    Splitting and merging
      Splitting and merging time: 00:00:00
    Split 1, merged 6; iterating...
    Splitting and merging
      Splitting and merging time: 00:00:00
    Avg bucket size: 1.47006e+06 (target: 2205098)
    Getting block 1 of 8
      Reserving size (2205099) for bucket 1
      Calculating Z arrays for bucket 1
      Entering block accumulator loop for bucket 1:
      bucket 1: 10%
      bucket 1: 20%
      bucket 1: 30%
      bucket 1: 40%
      bucket 1: 50%
      bucket 1: 60%
      bucket 1: 70%
      bucket 1: 80%
      bucket 1: 90%
      bucket 1: 100%
      Sorting block of length 1494907 for bucket 1
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 1494908 for bucket 1
    Getting block 2 of 8
      Reserving size (2205099) for bucket 2
      Calculating Z arrays for bucket 2
      Entering block accumulator loop for bucket 2:
      bucket 2: 10%
      bucket 2: 20%
      bucket 2: 30%
      bucket 2: 40%
      bucket 2: 50%
      bucket 2: 60%
      bucket 2: 70%
      bucket 2: 80%
      bucket 2: 90%
      bucket 2: 100%
      Sorting block of length 1656980 for bucket 2
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 1656981 for bucket 2
    Getting block 3 of 8
      Reserving size (2205099) for bucket 3
      Calculating Z arrays for bucket 3
      Entering block accumulator loop for bucket 3:
      bucket 3: 10%
      bucket 3: 20%
      bucket 3: 30%
      bucket 3: 40%
      bucket 3: 50%
      bucket 3: 60%
      bucket 3: 70%
      bucket 3: 80%
      bucket 3: 90%
      bucket 3: 100%
      Sorting block of length 753383 for bucket 3
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 753384 for bucket 3
    Getting block 4 of 8
      Reserving size (2205099) for bucket 4
      Calculating Z arrays for bucket 4
      Entering block accumulator loop for bucket 4:
      bucket 4: 10%
      bucket 4: 20%
      bucket 4: 30%
      bucket 4: 40%
      bucket 4: 50%
      bucket 4: 60%
      bucket 4: 70%
      bucket 4: 80%
      bucket 4: 90%
      bucket 4: 100%
      Sorting block of length 1984563 for bucket 4
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 1984564 for bucket 4
    Getting block 5 of 8
      Reserving size (2205099) for bucket 5
      Calculating Z arrays for bucket 5
      Entering block accumulator loop for bucket 5:
      bucket 5: 10%
      bucket 5: 20%
      bucket 5: 30%
      bucket 5: 40%
      bucket 5: 50%
      bucket 5: 60%
      bucket 5: 70%
      bucket 5: 80%
      bucket 5: 90%
      bucket 5: 100%
      Sorting block of length 697559 for bucket 5
      (Using difference cover)
      Sorting block time: 00:00:01
    Returning block of 697560 for bucket 5
    Getting block 6 of 8
      Reserving size (2205099) for bucket 6
      Calculating Z arrays for bucket 6
      Entering block accumulator loop for bucket 6:
      bucket 6: 10%
      bucket 6: 20%
      bucket 6: 30%
      bucket 6: 40%
      bucket 6: 50%
      bucket 6: 60%
      bucket 6: 70%
      bucket 6: 80%
      bucket 6: 90%
      bucket 6: 100%
      Sorting block of length 2014418 for bucket 6
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 2014419 for bucket 6
    Getting block 7 of 8
      Reserving size (2205099) for bucket 7
      Calculating Z arrays for bucket 7
      Entering block accumulator loop for bucket 7:
      bucket 7: 10%
      bucket 7: 20%
      bucket 7: 30%
      bucket 7: 40%
      bucket 7: 50%
      bucket 7: 60%
      bucket 7: 70%
      bucket 7: 80%
      bucket 7: 90%
      bucket 7: 100%
      Sorting block of length 1820676 for bucket 7
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 1820677 for bucket 7
    Getting block 8 of 8
      Reserving size (2205099) for bucket 8
      Calculating Z arrays for bucket 8
      Entering block accumulator loop for bucket 8:
      bucket 8: 10%
      bucket 8: 20%
      bucket 8: 30%
      bucket 8: 40%
      bucket 8: 50%
      bucket 8: 60%
      bucket 8: 70%
      bucket 8: 80%
      bucket 8: 90%
      bucket 8: 100%
      Sorting block of length 1338034 for bucket 8
      (Using difference cover)
      Sorting block time: 00:00:00
    Returning block of 1338035 for bucket 8


    Settings:
      Output files: "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index.*.ht2"
      Line rate: 6 (line is 64 bytes)
      Lines per side: 1 (side is 64 bytes)
      Offset rate: 4 (one in 16)
      FTable chars: 10
      Strings: unpacked
      Local offset rate: 3 (one in 8)
      Local fTable chars: 6
      Local sequence length: 57344
      Local sequence overlap between two consecutive indexes: 1024
      Endianness: little
      Actual local endianness: little
      Sanity checking: disabled
      Assertions: disabled
      Random seed: 0
      Sizeofs: void*:8, int:4, long:8, size_t:8
    Input files DNA, FASTA:
      /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/s_mikatae.fa
    Reading reference sizes
      Time reading reference sizes: 00:00:00
    Calculating joined length
    Writing header
    Reserving space for joined string
    Joining reference sequences
      Time to join reference sequences: 00:00:00
      Time to read SNPs and splice sites: 00:00:00
    Using parameters --bmax 2205099 --dcv 1024
      Doing ahead-of-time memory usage test
      Passed!  Constructing with these parameters: --bmax 2205099 --dcv 1024
    Constructing suffix-array element generator
    Converting suffix-array elements to index image
    Allocating ftab, absorbFtab
    Entering GFM loop
    Exited GFM loop
    fchr[A]: 0
    fchr[C]: 3665295
    fchr[G]: 5885531
    fchr[T]: 8105456
    fchr[$]: 11760527
    Exiting GFM::buildToDisk()
    Returning from initFromVector
    Wrote 8122399 bytes to primary GFM file: /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index.1.ht2
    Wrote 2940136 bytes to secondary GFM file: /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index.2.ht2
    Re-opening _in1 and _in2 as input streams
    Returning from GFM constructor
    Returning from initFromVector
    Wrote 6060293 bytes to primary GFM file: /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index.5.ht2
    Wrote 2990286 bytes to secondary GFM file: /Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index.6.ht2
    Re-opening _in5 and _in5 as input streams
    Returning from HGFM constructor
    Headers:
        len: 11760527
        gbwtLen: 11760528
        nodes: 11760528
        sz: 2940132
        gbwtSz: 2940133
        lineRate: 6
        offRate: 4
        offMask: 0xfffffff0
        ftabChars: 10
        eftabLen: 0
        eftabSz: 0
        ftabLen: 1048577
        ftabSz: 4194308
        offsLen: 735033
        offsSz: 2940132
        lineSz: 64
        sideSz: 64
        sideGbwtSz: 48
        sideGbwtLen: 192
        numSides: 61253
        numLines: 61253
        gbwtTotLen: 3920192
        gbwtTotSz: 3920192
        reverse: 0
        linearFM: Yes
    Total time for call to driver() for forward index: 00:00:07


Then we can align our reads. Remember that, for this step, we are going to use the cleaned reads.


```bash
%%bash 

HISAT2="/Users/martahuertas/Desktop/Software/hisat2-2.2.1"
INDEX="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/Genome_index/s_mikatae_index"
PAIRED_READS1="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads1_paired.fastq"
PAIRED_READS2="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/Reads/Clean/s_mikatae_reads2_paired.fastq"

$HISAT2"/"hisat2 -x $INDEX -1 $PAIRED_READS1 -2 $PAIRED_READS2 -S "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/SAM/"s_mikatae_RNA.sam
```

    18347000 reads; of these:
      18347000 (100.00%) were paired; of these:
        376875 (2.05%) aligned concordantly 0 times
        17489018 (95.32%) aligned concordantly exactly 1 time
        481107 (2.62%) aligned concordantly >1 times
        ----
        376875 pairs aligned concordantly 0 times; of these:
          75950 (20.15%) aligned discordantly 1 time
        ----
        300925 pairs aligned 0 times concordantly or discordantly; of these:
          601850 mates make up the pairs; of these:
            471955 (78.42%) aligned 0 times
            95611 (15.89%) aligned exactly 1 time
            34284 (5.70%) aligned >1 times
    98.71% overall alignment rate


### 2.1. Convert to BAM

In order to use `stringtie` and other programs (`IGV` for instance), we are going to need the bam. More specifically, we are getting the sorted and indexed bam.


```bash
%%bash 
SAMTOOLS="/Users/martahuertas/Desktop/Software/samtools-1.14"
SAM="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/SAM/s_mikatae_RNA.sam"
BAM="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/SAM/s_mikatae_RNA.bam"
SORTED_BAM="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/SAM/s_mikatae_RNA_sorted.bam"

#First we obtain the bam file
$SAMTOOLS"/"samtools view -Su $SAM > $BAM

#And then we sort and index it. 
$SAMTOOLS"/"samtools sort $BAM -o $SORTED_BAM
$SAMTOOLS"/"samtools index $SORTED_BAM
```

    [bam_sort_core] merging from 11 files and 1 in-memory blocks...


## 3. Transcript assembly 

We are using `stringtie` with the *S .mikatae* genome as reference to perform the transcript assemby and generate the *S. mikatae* transcriptome. 


```bash
%%bash

STRINGTIE="/Users/martahuertas/Desktop/Software/stringtie"
GENOME="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/s_mikatae_fixed.gff"
BAM="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/SAM/s_mikatae_RNA_sorted.bam"
OUTPUT_GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF"

$STRINGTIE"/"stringtie $BAM -G $GENOME --rf  -o $OUTPUT_GFF"/"s_mikatae_transcriptome.gff
```

## 4. Data analysis
Separating between annotated and novel genes.
First I am creating one gff for novel genes and other one for annotated.


```bash
%%bash 

GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF"


awk '($3 == "transcript") {print}' $GFF"/s_mikatae_transcriptome.gff" | grep "reference" > s_mikatae_knowngenes.gff

awk '($3 == "transcript") {print}' s_mikatae_transcriptome.gff | grep -v "reference" > s_mikatae_novelgenes.gff
```


```bash
%%bash 

GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF"

echo "Number of annotated transcripts"
awk '($3 == "transcript") {print}' $GFF"/s_mikatae_transcriptome.gff" | grep "reference" | wc -l
```

    Number of annotated transcripts
        5921



```bash
%%bash 

GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF"

echo "Number of novel transcripts"
awk '($3 == "transcript") {print}' $GFF"/s_mikatae_transcriptome.gff" | grep -v "reference" | wc -l
```

    Number of novel transcripts
         499


Getting fasta files from gffs:


```bash
%%bash

GENOME_FA="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Imput/s_mikatae.fa"
NOVEL_GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF/s_mikatae_novelgenes.gff"
ANNOT_GFF="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF/s_mikatae_knowngenes.gff"
OUTPUT="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/"

/Users/martahuertas/Desktop/Software/bedtools-2.17.0/bin/bedtools getfasta -fi $GENOME_FA -bed $ANNOT_GFF -s -fo $OUTPUT"/"s_mikatae_known.fasta
/Users/martahuertas/Desktop/Software/bedtools-2.17.0/bin/bedtools getfasta -fi $GENOME_FA -bed $NOVEL_GFF -s -fo $OUTPUT"/"s_mikatae_novel.fasta
```

### 4.1. Identify ORFs in the sequences
Using the script given in the Campus Global `get_longestORF_multi_withrandom.pl`:


```bash
%%bash 

get_ORF="/Users/martahuertas/Desktop/PGB/Project/get_longestORF_multi_withrandom.pl"
OUTPUT_DIR="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs"


#Imput file must be written in the text file of the get_longestORF_multi_withrandom.pl
perl $get_ORF > $OUTPUT_DIR"/"novel_orfs.txt $OUTPUT_DIR"/"novel_random_orfs.txt
#perl $get_ORF > $OUTPUT_DIR"/"known_orfs.txt $OUTPUT_DIR"/"known_random_orfs.txt

```


```bash
%%bash

echo "Novel genes ORFs:"
cat "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/novel_orfs.txt" | wc -l 
```

    Novel genes ORFs:
         499



```bash
%%bash

echo "Known genes ORFs:"
cat "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/known_orfs.txt" | wc -l 
```

    Known genes ORFs:
        6131



```bash
%%bash 
echo "Novel genes maximum lengths:"
cat "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/novel_orfs.txt" | cut -f2 | sort -nr | head
```

    Novel genes maximum lengths:
    6645
    4485
    4026
    4026
    3921
    3723
    3651
    3540
    3537
    3411



```bash
%%bash 
echo "Known genes maximum lengths:"
cat "/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/known_orfs.txt" | cut -f2 | sort -nr | head
```

    Known genes maximum lengths:
    14730
    12273
    11235
    9807
    9438
    9270
    9231
    8877
    8361
    8019


### 4.2. Get sequence length
I am using `getseqlen.pl` to obtain sequence length of the novel and knonw genes.


```bash
%%bash

get_length="/Users/martahuertas/Desktop/PGB/Project"
NOVEL_F="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/s_mikatae_novel.fasta"
KNOWN_F="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/s_mikatae_known.fasta"
OUTPUT_DIR="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA"

perl $get_length"/"getseqlen.pl $NOVEL_F $OUTPUT_DIR"/"novel.fasta_length
perl $get_length"/"getseqlen.pl $KNOWN_F $OUTPUT_DIR"/"known.fasta_length
```


```bash
%%bash 

NOVEL_LEN="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/novel.fasta_length"

echo "Novel fasta lengths"
sort -k2 -nr $NOVEL_LEN | head
```

    Novel fasta lengths
    >Smik_4:1077502-1090859(-)	13357
    >Smik_4:1077502-1090859(-)	13357
    >Smik_4:1077502-1090809(-)	13307
    >Smik_1:171578-183245(+)	11667
    >Smik_4:1080998-1090809(-)	9811
    >Smik_4:1080998-1090809(-)	9811
    >Smik_13:493697-503201(+)	9504
    >Smik_1:1550-5979(+)	8859
    >Smik_12:792339-800279(+)	7940
    >Smik_12:911735-919624(+)	7889



```bash
%%bash 

KNOWN_LEN="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/known.fasta_length"

echo "Known fasta lengths"
sort -k2 -nr $KNOWN_LEN | head
```

    Known fasta lengths
    >Smik_12:328888-343618(-)	14730
    >Smik_11:525588-537861(-)	12273
    >Smik_8:279198-290433(+)	11235
    >Smik_4:1291054-1300861(+)	9807
    >Smik_12:39240-48678(-)	9438
    >Smik_2:502017-511287(-)	9270
    >Smik_15:135386-144617(+)	9231
    >Smik_12:286689-295566(-)	8877
    >Smik_2:44652-53013(-)	8361
    >Smik_7:125926-133945(+)	8019

Once I have extracted the sequence length of the novel and known genes using `getseqlen.pl`, I am going to compare both types of genes:

```{r}
known <-read.table("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/known.fasta_length",header=TRUE)
novel <-read.table("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/FASTA/novel.fasta_length",header=TRUE)

known$type='known'
novel$type='novel'

```

```{r}
summary(known)
```
      name               length     
 Length:6131        Min.   :   24  
 Class :character   1st Qu.:  630  
 Mode  :character   Median : 1140  
                    Mean   : 1379  
                    3rd Qu.: 1818  
                    Max.   :14730  
                    
```{r}
summary(novel)
```
    name               length     
 Length:499         Min.   :  207  
 Class :character   1st Qu.:  612  
 Mode  :character   Median :  978  
                    Mean   : 1580  
                    3rd Qu.: 1784  
                    Max.   :13357 
                    
                    
This summaries are already indicating that novel genes are a slightly shorter. But we will see it visually by plotting the data into a histogram. We will compare the length of the sequences with a Mann-Whitney-Wilcoxon test, which is non-parametric (it does not assume a given distribution).

```{r}
wilcox.test(known$length,novel$length)
```
  Wilcoxon rank sum test with continuity correction
    
    data:  known$length and novel$length
    W = 1580273, p-value = 0.2186
    alternative hypothesis: true location shift is not equal to 0
  
This is the first plot in the ppt:
```{r}
library(ggplot2)

p1<-ggplot(data=known,aes(x=length, fill="known")) +
geom_histogram()

p1 <- p1 + geom_histogram(data=novel, aes(x=length, fill="novel"))+
labs(title="Transcript length",
      x="Length (nt)",
      y="Count",
     fill="Dataset")+
  scale_fill_brewer(palette="Set2")+
  scale_x_continuous(limits = c(0,7500)) +
    theme_classic()
p1
```

### 4.3. Expression analysis
We will now compare the expression levels (TPM) . We need to extract the information from the gff files. We will use a parser from an R library:

For me the `refGenome`library is not working, therefore I am going to use `readGFF` from the package `rtracklayer`. 

```{r}
library("rtracklayer")

known_gff <- readGFF("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF/s_mikatae_knowngenes.gff")
novel_gff <- readGFF("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/GFF/s_mikatae_novelgenes.gff")

#transforming the TPM column in numeric values to test them statistically.
known_gff$TPM <- as.numeric(known_gff$TPM)
novel_gff$TPM <- as.numeric(novel_gff$TPM)

known_gff$t_type='known'
novel_gff$t_type='novel'

```

Getting the summary information of known expression values.
```{r}
summary(known_gff$TPM)
```
  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
    0.00    10.87    29.89   123.07    82.01 12494.49 

```{r}
summary(novel_gff$TPM)
```
  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
   2.196   11.417   31.478  489.520  217.147 6323.138 


Check if the distributions are significantly different
```{r}
wilcox.test(known_gff$TPM, novel_gff$TPM)
```


	Wilcoxon rank sum test with continuity correction

data:  known_gff$TPM and novel_gff$TPM
W = 1357827, p-value = 2.918e-05
alternative hypothesis: true location shift is not equal to 0



It was shocking that, from the mean, it can be thought that novel genes have higher expression than known. When looking at the median this values look more similar. When looking at the plot, it becomes clear that this values are generated by a high amount of known genes with lower expression than expected. 

```{r}
p2 <- ggplot(data=novel_gff, aes(x=t_type, y=log2(TPM), fill="Novel")) +
geom_violin(alpha=0.7) +
  geom_boxplot(width=0.05)

p2 <- p2 + geom_violin(data=known_gff, aes(x=t_type, y=log2(TPM), fill="Known"),alpha=0.5)+
  scale_fill_brewer(palette="Set2")+
   labs(title="Expression",
        x="Dataset",
        y="Log2(TPM)")+
  guides(fill = "none")+
  theme_classic()

p2 <- p2 +  geom_boxplot(data=known_gff, aes(x=t_type, y=log2(TPM), fill="Known"), width=0.05)


p2
```

                    
## 5. Aa composition analysis

I am using `transeq`from EMBOSS to translate the orfs. But first we need to obtain the ORFs in fasta format. For that I have modified the perl script we used to obtain the longest ORF lengths in order to get the fasta and not the length. 


```bash
%%bash 

#Getting novel and known longest ORFs and random ORFs in fasta.
LONGEST_ORF="/Users/martahuertas/Desktop/PGB/Project/get_longestORF_multi_withrandom.pl"
OUTPUT_DIR="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs"


#Imput file must be written in the text file of the get_longestORF_multi_withrandom.pl
perl $LONGEST_ORF > $OUTPUT_DIR"/"novel_orfs.fasta $OUTPUT_DIR"/"novel_random_orfs.fasta
#perl $LONGEST_ORF > $OUTPUT_DIR"/"known_orfs.fasta $OUTPUT_DIR"/"known_random_orfs.fasta


```


```bash
%%bash
#Now I am going to translate them by using transeq.
EMBOSS="/Users/martahuertas/Desktop/Software/EMBOSS-6.6.0/emboss"
NOVEL="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/novel_orfs.fasta"
NOVEL_R="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/novel_random.fasta"
KNOWN="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/known_orfs.fasta"
KNOWN_R="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/known_random.fasta"
OUTDIR="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas"

$EMBOSS"/"transeq -sequence $NOVEL_R -trim -outseq $OUTDIR"/"novel_r_orfs_aa.fasta
$EMBOSS"/"transeq -sequence $KNOWN_R -trim -outseq $OUTDIR"/"known_r_orfs_aa.fasta
$EMBOSS"/"transeq -sequence $NOVEL -trim -outseq $OUTDIR"/"novel_orfs_aa.fasta 
$EMBOSS"/"transeq -sequence $KNOWN -trim -outseq $OUTDIR"/"known_orfs_aa.fasta 
```

    Translate nucleic acid sequences
    Translate nucleic acid sequences


To study this translated orfs I am going to use some programs that are integrated with EMBOSS: `Pepinfo` and `Pepstats`, for instance.

- `Pepstats` -> Calculate properties of protein sequences such as molecular weight


```bash
%%bash

#Generating peptide statistics with pepstats

EMBOSS="/Users/martahuertas/Desktop/Software/EMBOSS-6.6.0/emboss"
NOVEL="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/novel_orfs_aa.fasta"
NOVEL_R="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/novel_r_orfs_aa.fasta"
KNOWN="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/known_orfs_aa.fasta"
KNOWN_R="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/known_r_orfs_aa.fasta"
OUTDIR="/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein"

$EMBOSS"/"pepstats -sequence $NOVEL -outfile  $OUTDIR"/"novel.pepstats
$EMBOSS"/"pepstats -sequence $NOVEL_R -outfile  $OUTDIR"/"novel_random.pepstats
$EMBOSS"/"pepstats -sequence $KNOWN -outfile  $OUTDIR"/"known.pepstats
$EMBOSS"/"pepstats -sequence $KNOWN_R -outfile  $OUTDIR"/"known_random.pepstats
```

    Calculate statistics of protein properties
    Calculate statistics of protein properties
    Calculate statistics of protein properties
    Calculate statistics of protein properties


We could be interested in **isoelectric point** (IP), **molecular weight** (MW), and **peptide length** to compare novel and known peptides. Therefore I am going to obtain that information from the pepstats files.


```python
novel=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/novel.pepstats").readlines()
novel_r=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/novel_random.pepstats").readlines()
known=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/known.pepstats").readlines()
known_r=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/known_random.pepstats").readlines()

```


```python
def getPepStatistics (pepstats):
    MW=[]
    name=[]
    PI=[]
    length=[]
    dic={}
    for line in pepstats:
        if "Molecular weight" in line:
            MW.append(float(line.split(" ")[3]))
        if "PEPSTATS" in line:
            name.append(line.split(" ")[2])
            length.append(int(line.split(" ")[6].split("\n")[0]))
        if "Isoelectric" in line:
            PI.append(float(line.split(" ")[3].split("\n")[0]))
    for i in range(0,len(name)):
        dic[name[i]] = [length[i],MW[i],PI[i]]
    return dic, length, MW, PI
        
novel_dic=getPepStatistics(novel)[0]
novel_l=getPepStatistics(novel)[1]
novel_MW=getPepStatistics(novel)[2]
novel_PI=getPepStatistics(novel)[3]
novel_r_dic=getPepStatistics(novel_r)[0]
known_dic=getPepStatistics(known)[0]
known_l=getPepStatistics(known)[1]
known_MW=getPepStatistics(known)[2]
known_PI=getPepStatistics(known)[3]
known_r_dic=getPepStatistics(known_r[0])
```

Now I am going to print the output in a tsv file:


```python
import csv 
novel_out=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/novel_aa_statistics.tsv", "w")
novel_r_out=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/novel_r_aa_statistics.tsv", "w")
known_out=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/known_aa_statistics.tsv", "w")
known_r_out=open("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/known_r_aa_statistics.tsv", "w")

def outFile (dic, output):
    writer=csv.writer(output, delimiter="\t")
    writer.writerow(["orf_name", "Length", "Mol_weight", "PI"])
    for key in dic.keys():
        writer.writerow([key, dic[key][0], dic[key][1], dic[key][2]])
            
outFile(novel_dic, novel_out)
novel_out.close()
outFile(novel_r_dic, novel_r_out)
novel_r_out.close()
outFile(known_dic, known_out)
known_out.close()
outFile(known_r_dic, known_r_out)
known_r_out.close()
```

Getting the mean length, mean molecular weight and mean IP.

Now I am going to use the pepstat output to get the 3 most present aa in each dataset. 


```python
AA=["A =", "B =", "C =", "D =", "F =", "G =", "H =", "I =", "J =", "K =", "L =", "M =", "N =", "O =","P =", "Q =", "R =", "S =", "T =", "U =", "V =", "W =","X =", "Y ", "Z ="]

aa_values_novel=[]
for line in novel:
    for aa in AA:
        if line.startswith(aa):
            aa_values_novel.append(float(line.split("		")[2]))
            
            
aa_values_known=[]
for line in known:
    for aa in AA:
        if line.startswith(aa):
            aa_values_known.append(float(line.split("		")[2]))

#With this I am obtaining a list with the positions of the most common aa in each sequence. Comparing those positions with the list of aa we can infer which are the three most common aa.            
def commonAA (aa_values):
    position_list=[]
    for n in range(0,len(aa_values),25):
        maximum=aa_values[0]
        maximum2=aa_values[0]
        maximum3=aa_values[0]
        aa_position=0
        aa_position2=0
        aa_position3=0

        j=0
        for i in range(n,n+25):
            j+=1
            if aa_values[i] > maximum:
                maximum=aa_values[i]
                aa_position=j
        position_list.append(aa_position)

        j=0
        for i in range(n,n+25):
            j+=1
            if aa_values[i] > maximum2 and aa_values[i]!=maximum:
                maximum2=aa_values[i]
                aa_position2=j
        position_list.append(aa_position2)

        j=0
        for i in range(n,n+25):
            j+=1
            if aa_values[i] > maximum3 and aa_values[i]!=maximum and aa_values[i]!=maximum2:
                maximum3=aa_values[i]
                aa_position3=j
        position_list.append(aa_position3) 
    return position_list

positionList_novel=commonAA(aa_values_novel)
positionList_known=commonAA(aa_values_known)
```

Now I am going to get the position values obtained before. Order them and obtain the number of times each position appears in the list. Then, by choosing the five with more appareances we can know the most common aa. 


```python
import itertools
#First value is the position and second the number of times it appears in the list. 
commonAA_values_novel = [(valor, sum(1 for x in ocurrencias)) for valor, ocurrencias in itertools.groupby(sorted(positionList_novel))]

commonAA_values_known = [(valor, sum(1 for x in ocurrencias)) for valor, ocurrencias in itertools.groupby(sorted(positionList_known))]
```


```python
print("The most common AA in the novel dataset are:")
print("M", commonAA_values_novel[9])
print("T", commonAA_values_novel[15])
print("J", commonAA_values_novel[7])
print("L", commonAA_values_novel[8])
```

    The most common AA in the novel dataset are:
    M (11, 282)
    T (18, 170)
    J (8, 166)
    L (10, 152)



```python
print("The most common AA in the known dataset are:")
print("M", commonAA_values_known[9])
print("T", commonAA_values_known[15])
print("L", commonAA_values_known[8])
print("J", commonAA_values_known[7])
```

    The most common AA in the known dataset are:
    M (11, 4217)
    T (18, 3052)
    L (10, 2450)
    J (8, 1500)

The analysis of this data has been performed with R as well:

```{r}
novel_aa_stats <-read.table("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/novel_aa_statistics.tsv",header=TRUE)
known_aa_stats <- read.table("/Users/martahuertas/Desktop/PGB/Project/s_mikatae/Output/ORFs/fastas/protein/pepstatistics/known_aa_statistics.tsv",header=TRUE)


novel_aa_stats$type='novel'
known_aa_stats$type='known'

summary(novel_aa_stats)

```

orf_name             Length         Mol_weight             PI             type          
 Length:485         Min.   :   2.0   Min.   :   262.4   Min.   : 3.245   Length:485        
 Class :character   1st Qu.:  37.0   1st Qu.:  4324.2   1st Qu.: 6.046   Class :character  
 Mode  :character   Median :  62.0   Median :  7011.9   Median : 8.536   Mode  :character  
                    Mean   : 163.3   Mean   : 18453.6   Mean   : 8.230                     
                    3rd Qu.: 155.0   3rd Qu.: 17561.5   3rd Qu.:10.358                     
                    Max.   :2214.0   Max.   :245082.9   Max.   :13.461 
                    
```{r}
summary(known_aa_stats)
```

orf_name             Length         Mol_weight             PI             type          
 Length:5810        Min.   :   1.0   Min.   :   149.2   Min.   : 3.245   Length:5810       
 Class :character   1st Qu.: 217.2   1st Qu.: 24670.3   1st Qu.: 5.517   Class :character  
 Mode  :character   Median : 393.0   Median : 44396.6   Median : 7.091   Mode  :character  
                    Mean   : 467.3   Mean   : 52899.4   Mean   : 7.342                     
                    3rd Qu.: 617.0   3rd Qu.: 69725.8   3rd Qu.: 9.128                     
                    Max.   :4909.0   Max.   :558500.1   Max.   :13.105                     
   
                    
Generating density plots with the novel and known length values:

```{r}
p3<-ggplot(data=novel_aa_stats,aes(x=Length, fill="Novel")) +
geom_density(alpha=0.8) 

p3 <- p3 + geom_density(data=known_aa_stats, aes(x=Length,  fill="Known"),alpha=0.5)+
  scale_fill_brewer(palette="Set2")+
  scale_x_continuous(limits = c(0,2000))+
  labs(title="Protein length",
         x="Length (aa)",
         y="Frequency",
         fill="Dataset")+
  theme_classic()

p3

```

We should test if the known orfs have statistically bigger (speaking of peptide length) peptides than the novel ones. Visually, in this case it is more or less clear that known ORFs are larger.

Novel vs. known length
```{r}
wilcox.test(novel_aa_stats$Length,known_aa_stats$Length)
```
  Wilcoxon rank sum test with continuity correction

    data:  novel_aa_stats$Length and known_aa_stats$Length
    W = 530172, p-value < 2.2e-16
    alternative hypothesis: true location shift is not equal to 0

Now I am going to statistically compare the known with the novel isoelectric points.

```{r}
p4 <- ggplot(data=novel_aa_stats,aes(y=PI, x=type, fill="Novel")) +
geom_violin(alpha=0.7) +
  geom_boxplot(width=0.05)

p4 <- p4 + geom_violin(data=known_aa_stats, aes(y=PI,x=type, fill="Known"),alpha=0.5)+
  scale_fill_brewer(palette="Set2")+
   labs(title="Isoelectric point",
        x="Dataset",
        y="Isoelectric point")+
  guides(fill = "none")+
  theme_classic()

p4 <- p4 + geom_boxplot(data=known_aa_stats, aes(y=PI,x=type, fill="Known"), width=0.05)
p4
```

Novel observed vs. known observed
```{r}
wilcox.test(novel_aa_stats$PI,known_aa_stats$PI)
```

  Wilcoxon rank sum test with continuity correction
  
    data:  novel_aa_stats$PI and known_aa_stats$PI
    W = 1715570, p-value = 1.523e-15
    alternative hypothesis: true location shift is not equal to 0
  



## 6. BLAST with novel proteins
We are going to perform a blast search comparing the novel peptides obtained by using `stringtie` and `get_longestORF_multi_withrandom`. 

To do so, we have obtained *S. cereviseae* proteome from [Ensebl](https://www.ensembl.org/Saccharomyces_cerevisiae/Info/Index). We are going to create a BLAST db and then compare the novel peptides with the annotated proteins in *S. cerevisiae* knowing that this specie has a well annotated proteome. 

As a control we will also use *Candida albicans* proteome in [Ensebl](https://www.ensembl.org/Candida_albicans/Info/Index) to compare our transcriptome with a more phylogenetically distant specie. 

```{bash}
# create the databases to compare our S. mikatae orfs
makeblastdb -in Saccharomyces_cerevisiae.R64-1-1.pep.all.fa -dbtype prot -parse_seqids -out cerevisiae.db

makeblastdb -in Candida_albicans_19f_gca_000775445.Cand_albi_19F_V2.pep.all.fa -dbtype prot -parse_seqids -out albicans.db
```
```{bash}
# run blastp and acquire the desired output
blastp -query novel_orfs_aa.fasta -db cerevisiae.db -out mikataeVScerevisiae_default.out -outfmt "6 qseqid sseqid pident evalue bitscore" -evalue 0.05

blastp -query novel_orfs_aa.fasta -db albicans.db -out mikataeVSalbicans_default.out -outfmt "6 qseqid sseqid pident evalue bitscore" -evalue 0.05
```

As can be seen in the code we used a e-value threshold of 0.05.

Then we got those proteins with or without hits with *S. cerevisiae* and from those the ones that have hit with *C. albicans*. We have two aims to perform this analysis. 

The first is to see whether the novel peptides found have homologues in other species and, therefore, are not strickly novel and should have been annotated in *S. mikatae* before. 

The second aim is to compare the length of those peptides that had no-hits, which we will consider putative novel in *S. mikatae*. Those that have hits in *S. cerevisiae* but not in *C. albicans* which will be considered novel in the *Saccharomyces* sensu stricto group. And the ones with hits in both species that will be considered as more conserved.

We expect the novels to be shorter than the ones in the sensu stricto group, which should be shorter than the ones conserved. 

```{bash}
## Parsing

# script to only retrieve the hit for every ORF with better e-value

awk 'BEGIN { C=$1; i=0; OFS="\t" }
{
  if (C==$1){
    if (i==0) {
    E=$4;
    i=1;
    }
  }
  else {
    print $1, $2, $3, $4, $5;
    i=0;
    C=$1
  }
}' mikataeVScerevisiae_default.out > mikataeVScerevisiae_parsed.out

wc -l *parsed* | tee
```

203 mikataeVSalbicans_parsed.out
260 mikataeVScerevisiae_parsed.out


```{bash}
# script to obtain a file with the ORFs that didn't produce a hit against S. cerevisiae

awk 'FNR==NR{ a[$1]; next } !($1 in a)' mikataeVScerevisiae_parsed.out orfs_sols.txt > nohit_cerevisiae.tab

# make a file with all the coordinates of every ORFs (S. mikatae)

awk '$0 ~ /^>/ { gsub(/>/, "", $1);print $1; }' novel_orfs_aa.fasta > orfs_sols.tab

# retrieve the hits only for S. cerevisiae

awk 'FNR==NR{ a[$1]; next } !($1 in a)' mikataeVSalbicans_parsed.out mikataeVScerevisiae_parsed.out > only_cerevisiae.tab

# retrieve the hits shared by S. cerevisiae and C. albicans

awk 'FNR==NR{ a[$1]; next } !($1 in a)' only_cerevisiae.txt mikataeVScerevisiae_parsed.out > hit_cere_mika_albi.tab

```

```{bash}
wc -l only_cerevisiae.tab
cat only_cerevisiae.tab | head -6
```
60 only_cerevisiae.tab

  coordinates
121529 122579
179563 180043
314519 316481
347805 351601
446379 447425

```{bash}
wc -l nohit_cerevisiae.tab
cat nohit_cerevisiae.tab | head -6
```
226 nohit_cerevisiae.tab

  coordinates
1550 5979
152906 154119
155648 158646
164814 165200
28755 30009

```{bash}
wc -l hit_cere_albi.tab
cat hit_cere_albi.tab | head -6
```
200 hit_cere_albi.tab

  coordinates
70180 70673
125170 126882
128809 134373
171578 183245
40316 42431

To sum up, the results show as follows:

| Novel (no-hit)  | Sensu stricto (hit in S. cerevisiae)  | Conserved (hit in S. cer and C. alb) |
|-----------------|---------------------------------------|--------------------------------------|
|      226        |                 60                    |                  200                 | 
 

This is the script to generate the violin plot of the lengths:

```python
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

data1 = pd.read_table("/home/maria/Desktop/PGB/blastp/graphs/novel_orfs.tab")
df_orfs = pd.DataFrame(data1)
df_orfs[['cI','cF']] = df_orfs.coordinates.str.split(expand=True,)
df_orfs['cI'] = df_orfs['cI'].astype(int)
df_orfs['cF'] = df_orfs['cF'].astype(int)
df_orfs['length'] = df_orfs['cF'] - df_orfs['cI']


data2 = pd.read_table("/home/maria/Desktop/PGB/blastp/graphs/nohit_cerevisiae.tab")
df_nohit = pd.DataFrame(data2)
df_nohit[['cI','cF']] = df_nohit.coordinates.str.split(expand=True,)
df_nohit['cI'] = df_nohit['cI'].astype(int)
df_nohit['cF'] = df_nohit['cF'].astype(int)
df_nohit['length'] = df_nohit['cF'] - df_nohit['cI']

data3 = pd.read_table("/home/maria/Desktop/PGB/blastp/graphs/only_cerevisiae.tab")
df_only = pd.DataFrame(data3)
df_only[['cI','cF']] = df_only.coordinates.str.split(expand=True,)
df_only['cI'] = df_only['cI'].astype(int)
df_only['cF'] = df_only['cF'].astype(int)
df_only['length'] = df_only['cF'] - df_only['cI']

data4 = pd.read_table("/home/maria/Desktop/PGB/blastp/graphs/hit_cere_albi.tab")
df_hitcerealbi = pd.DataFrame(data4)
df_hitcerealbi[['cI','cF']] = df_hitcerealbi.coordinates.str.split(expand=True,)
df_hitcerealbi['cI'] = df_hitcerealbi['cI'].astype(int)
df_hitcerealbi['cF'] = df_hitcerealbi['cF'].astype(int)
df_hitcerealbi['length'] = df_hitcerealbi['cF'] - df_hitcerealbi['cI']


dataplot = [df_orfs["length"], df_nohit["length"], df_only["length"], df_hitcerealbi["length"]]

headers = ["ORFs S. mikatae", "No hits with S. cerevisiae", "Hits only with S. cerevisiae", "Hits with S.cerevisiae and C. albicans"]

df = pd.concat(dataplot, axis=1, keys=headers)


plt.figure(figsize = (12,8.25))

ax = sns.violinplot(data=df, orient="v", keys=headers, palette=['coral', 'peachpuff', 'sandybrown', 'darkorange'])

plt.xticks(rotation=12.5)
ax.set_title("Length of ORFs by type")
ax.set_ylabel("Length (nº nucleotides)")

plt.savefig('blastp_output.png')
```



