# Base R DNA Sequence Analysis

## Identifying Information

- Programmer: Nita Hegde
- Language: R (base R only)
- Version: 1.0
- Submission date: 9th November, 2026
- Script name: dna_analysis.R

## Purpose

This program reads the chr1_GL383518v1_alt DNA sequence from a
FASTA file, prints specified positions, creates its reverse
complement, and counts letters in consecutive 1000-base blocks.

It also creates a data frame with A, C, G, and T counts for
each block and calculates the sum of each row.

## Files Needed

- dna_analysis.R
- chr1_GL383518v1_alt.fa
- README.md

The FASTA file is the sequence downloaded for Assignment 1.
Keep it in the same working directory as the script.

## Required Software

R is required. RStudio is optional.
No additional R packages are needed.

## How to Run

In R or RStudio, set the working directory to the folder
containing the script and FASTA file, then run:

source("dna_analysis.R")

Alternatively, from a terminal in that folder, run:

Rscript dna_analysis.R

## Output

The program prints:

- The original sequence's 10th and 758th letters.
- The reverse complement's 79th letter.
- Positions 500 through 800 of the reverse complement.
- Letter counts for the first kilobase.
- A four-column nucleotide-count data frame.
- Row sums and rows whose sums differ from 1000.

The program creates:

- kilobase_counts.rds: letter-count lists for every block.
- nucleotide_counts.csv: A, C, G, and T counts per block.
- row_sums.csv: the nucleotide-count sum for every block.

## Interpretation

Complete blocks containing only A, C, G, and T have a row
sum of 1000. The final block may be shorter. Unknown or
ambiguous bases are counted in the list but excluded from
the four-column data frame.

## Data Source

UCSC Genome Browser, human hg38 chromosome files:
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/

Input file: chr1_GL383518v1_alt.fa