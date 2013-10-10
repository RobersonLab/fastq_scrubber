fastq_scrubber
==============

## Purpose
fastq_scrubber is a horrible beast of a python program. Derived from a need to perform QC on paired-end Illumina reads. These functions are performed better and almost certainly faster with [cutadapt](http://code.google.com/p/cutadapt/), but it currently cannot handle paired-end reads simultaneously. As a result of cleaning the files separately, the read orders get out of sync. Reads will always remain in sync with this tool. Designed for Illumina sequencing in either the +64 or +33 offset, but not the log-odds scale.

## Basic QC

### Identifier sync
Since two files are processed at once, the identifiers for each read are compared. If the reads are out of sync at any point, the program exits with an error. To function properly the identifier base (unique ID, coordinates, etc) must be separated from the mate identifier and extra information by either a '#' character (CASAVA < v1.8; e.g. @MACHINE:6:1:100:1000#0/1) or a space (CASAVA >= v1.8; e.g. @MACHINE:RUN:FLOWCELLID:1:100:1000:10000 1:Y:18:TGAGGC)

### High N content
If the total N content of **both** reads is > a threshold (85% in this version), the reads are dumped to the junk files (described in Output section).

### Mismatched sequence / quality lengths
Each read is tested for equal length between the sequence and quality of read 1, and the sequence and quality in read 2. This is only tested **within** on end, so the ends can have different lengths. If the sequence and qualities are not matched length for either end, both ends are written to the junk file.

## Advanced QC
Reads passing the most basic QC steps have advanced QC performed. Steps are listed *in the order they are applied*.

### Hard clip ends
There are options (see program Usage statement) to hard clip bases from either read 1 or read 2 and either end in any combination. If you have a leading N on all of your read1 sequence, you would set "r1_ltrim=1" to get rid of it. All hard trimming functions are processed before any other QC steps. The default is no hard clipping on either end of either read.

### Terminal Ns
There should be **no** N nucleotides in Illumina sequencing. Obviously it happens, but is a result of errors or sequence so noisy that the base can't be called. These can't be aligned to the genome and aren't usable, so they are removed (if possible) from read ends.

#### 5' end
It would be a bad idea to **hard clip** a series of N bases from the 5' end of a paired-end read. While these bases can't be aligned, removing them could cause trouble with calling duplicates. Instead of hard clipping 5' N nucleotides, they remain and the quality is floored to 0 in the proper encoding.

#### 3' end
Any **sequential** Ns from the 3' end of a read are removed.

### Trim low-quality 3' bases
Bases at the 3' read end with *consecutive* qualities less than the cutoff (default 5) are removed. The cutoff can be adjusted, and quality trimming can be completely disabled. When specifying a new cutoff, use a numeric quality score with no offset. The script automatically adjusts based on input scale.

### Adapter removal
If enabled, adapter sequence will be removed from the 3' end of each read.

**Pay close attention to your *actual* adapters.** The defaults used here may have no utility for your specific application. Specify them at the command-line (comma separated, see program usage statement).

#### Removal algorithm
This is an extremely hacky strategy. To avoid globally aligning each read, the adapters are split into seed sequences. The short seeds from the adapter are used to generate regular expressions. Each read is checked for a match to the seeds. Since there are multiple seeds, they are offset. If a match is found, the position is reset from the beginning of the match to where the adapter should actually start based on the offset for the seed. Then an ungapped alignment is performed either until the end of the read or until the end of the adapter, whichever is achieved first. N characters are considered matches. After alignment if the read region has at least 85% identity with the adapter, it is removed.

### Post-QC checks
It is possible for an entire read to contain adapter contamination, Ns, or low quality sequence in some combination to cause the final clean read length to be 0. If the clean read is < a minimum size (5 bp default) then the remaining length is padded by Ns with the lowest possible quality score. The read is also flagged as "bad".

If both reads are flagged bad after QC, then they are written to the junk file. If only one is bad, then it remains and is written correctly to the clean fastq file.

## Input
Designed to work for paired-end reads; for single-end [cutadapt](http://code.google.com/p/cutadapt/) may fit your needs better.

Requires two FASTQ Illumina sequencing files with either the +64 (old) or +33 (new) quality offset. The scale of the input files is detected automatically.

**VERY IMPORTANT NOTE**
Read 1 and read 2 status is based **solely** on the order you specify them at the command line.

## Output
Four files derived from the specified base name. They are all FASTQ format.

    BASENAME_1_clean.ext
    BASENAME_2_clean.ext
    BASENAME_1_junk.ext
    BASENAME_2_junk.ext

The "clean" files contain all the post QC data in read order sync. The extension for these files depends on the *output* quality scale. If the scale is the Illumina +64 scale, the extension will be ".fq". If the output scale is +33, the output will be ".fqs" to signify it is on the Sanger / Phred quality scale.

The "junk" files contain all the failed reads (still properly paired). The *original* reads are printed here at their full length with no quality conversion if it was used. The extension depends on the input scale since the quality isn't converted, and follows the same convention the clean files use.

Run information, including options set, program version, and post QC info are printed to stdout. It is recommended to redirect that information to a log file.

## Example code

Two input files, myData_1.fq and myData_2.fq, that will be cleaned up with no adapter stripping.

    python fastq_scrubber.py myData_1.fq,myData_2.fq myCleanData 1>clean.log 2>&1
	
Two input files to be cleaned up with the few default adapters listed internally.

    python fastq_scrubber.py myData_1.fq,myData_2.fq myCleanData stripAdapter 1>clean.log 2>&1

Two input files to be cleaned up with custom adapters.

    python fastq_scrubber.py myData_1.fq,myData_2.fq myCleanData adapter=TGACAGCGAGCAGAGCCTCGA,CACATGTGTCGTCTCGGTCAA 1>clean.log 2>&1
	
Two input files to be cleaned with no adapter removal and no low-quality 3' end trimming

    python fastq_scrubber.py myData_1.fq,myData_2.fq myCleanData noRightQtrim 1>clean.log 2>&1
