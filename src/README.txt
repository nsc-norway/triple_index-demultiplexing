## Demultiplexing tool for internal dual barcodes


This tool processes two gzip'ed FASTQ files, and writes the reads in those files into
separate FASTQ files based on internal barcodes. The barcode sequences are trimmed from
the reads. Optionally also removes heterogeneity spacer sequences following the barcode
sequences.

# Usage 

     triple_index-demultiplexing/src/demultiplexer_long BARCODE_FILE SAMPLE_SHEET \
                INPUT_R1 INPUT_R2 OUTPUT_PREFIX \
                 [BARCODE_MISMATCHES_PER_READ=L1 \
                 [ALIGNMENT_MISMATCHES=BARCODE_MISMATCHES_PER_READ]]


## options

  * `BARCODE_FILE` File mapping barcode names to sequences. There are three tab-separated
    columns: Name, Barcode and Spacer. This file should contain all forward and reverse
    barcode sequences. It may contain more barcodes than are actually used. The spacer
    sequences are appended to the barcode sequences for the purpose of trimming.

Example `BARCODE_FILE`:

    f1	CCTAAACTACGG
    f2	TGCAGATCCAAC
    f4	GTGGTATGGGAG    T
    f5	ACTTTAAGGGTG    T
    f6	GAGCAACATCCT    T
    f7	TGTTGCGTTTCT    GT
    f8	ATGTCCGACCAA    GT
    f10	ACAGCCACCCAT    CGA
    f11	TGTCTCGCAAGC    CGA
    f13	GTTACGTGGTTG    ATGA
    f16	TACCGGCTTGCA    TGCGA
    f19	CACCTTACCTTA    GAGTGG
    f22	TTAACTGGAAGC    CCTGTGG
    r3	CCATCACATAGG
    r4	GTGGTATGGGAG    A
    r7	TGTTGCGTTTCT    TC
    r12	GAGGAGTAAAGC    CTA
    r15	CGTAAGATGCCT    GATA
    r16	TACCGGCTTGCA    ACTCA
    r19	CACCTTACCTTA    TTCTCT
    r22	TTAACTGGAAGC    CACTTCT


  * `SAMPLE_SHEET`: Whitespace separated file of sample name, barcode 1 and barcode 2. The
    barcode values refer to names in the BARCODE_FILE.

Example `SAMPLE_SHEET`:

    Sample01 f1 r3
    Sample02 f1 r4
    Sample03 f1 r7
    Sample04 f1 r12
    Sample05 f1 r15
    Sample06 f1 r16
    Sample07 f1 r19
    Sample08 f1 r22
    Sample09 f2 r3
    Sample10 f2 r4
    Sample11 f2 r7


  * `INPUT_R1` and `INPUT_R2`: Data files with inline barcodes and spacer sequences intact.
    Paired-end data is required, and the files must be gzip'ed FASTQ format.

  * `OUTPUT_PREFIX`: The output filename for each sample is constructed as:
    `OUTPUT_PREFIX + SAMPLE_NAME`. E.g. if the output prefix is equal to 
    `demultiplexed/`, then all files will be written to the directory demultiplexed,
    with filenames equal to the sample names (the slash at the end is required in order
    to specify a directory.

  * `BARCODE_MISMATCHES_PER_READ`: Number of mismatches to allow when assigning reads to
    samples. If the number is prefixed by "L", e.g. L1, the Levenshtein distance is used.
    If just a number is given, the Hamming distance is used -- no insertions or deletions
    are allowed. Note that the barcodes are not checked for uniqueness. If you specify a
    too large `BARCODE_MISMATCHES_PER_READ`, the assignment to samples may be ambigous,
    and the end result is undefined (the reads are usually assigned to the first sample
    in the list that matches). The default value is L1, i.e. to use Levenshtein distance
    with one edit operation allowed.

  * `ALIGNMENT_MISMATCHES`: Controls the pairwise alignment of the barcode+spacer sequence
    to the read, for the purpose of trimming. The number of edit operations allowed to match
    the read to the expected sequence. The solution with the fewest edit operations is always
    chosen, and if there are multiple solutions, the shortest match is chosen. This parameter
    is ignored if the Hamming distance is used for mismatches -- then the true length of the
    barcode+spacer is always trimmed. As the alignment and mismatch counting is done in one
    operation, the value of `ALIGNMENT_MISMATCHES` should always be greater or equal to
    `BARCODE_MISMATCHES_PER_READ`.


## output

FASTQ files are created for each sample, and a pair of files with name Undetermined, for reads
having barcodes that didn't match any of the samples. If there are no reads from a sample, no
output file is produced for that sample.

A summary is written to standard output after completion, with the following columns:

  * `SAMPLE_NAME`: As appears in the sample sheet. Undetermined is below the line.

  * `NUM_READS`: Number of read pairs assigned to this sample.

  * `PCT_READS`: The percentage of reads coming from this sample.

  * `PCT_PERFECT_BARCODE`: Percentage of the reads assigned to this sample that had no mismatches,
    in the barcode, neither in read 1 barcode or read 2 barcode.

  * `PCT_SPACER_FAIL`: Percentage of the reads assigned to this sample which failed to aligned the
    spacer sequence, and were thus untrimmed. This percentage is per end, so a read with forward
    spacer matched, and reverse spacer failed, counts as 50 %.


