# DNA-Tokenization-using-Byte-Pair-Encoding-BPE-
Large Language Models (LLMs) process text as meaningful tokens rather than raw characters. In genomics, tokenizing DNA naively by individual nucleotides (A, C, G, T) results in extremely long inputs that fail to capture biologically meaningful motifs, such as promoters, repeats, or CpG islands.   
This project implements Byte Pair Encoding (BPE) to overcome this limitation. BPE is a subword tokenization algorithm that iteratively discovers and merges the most frequently co-occurring sequence pairs into new tokens, enabling the automatic learning of frequent DNA subsequences without relying on fixed-length k-mers.

# Project Overview
The analysis is broken down into four core implementation tasks:
+ Data Preprocessing: Reads and formats Chromosome 20 from the human reference genome (hg38), stripping headers and converting it into a single continuous, uppercase DNA string.
+ Part 1: BPE Training: Implements ```train_bpe``` to learn a vocabulary of size ```K```. It iteratively counts adjacent pairs in the sequence and merges the most frequent pair into a single new token until the target vocabulary size is reached.
+ Part 2: DNA Tokenization: Implements ```tokenize_dna``` to parse new sequences using the learned BPE vocabulary, prioritizing the longest available matching tokens by sorting the vocabulary by length.
+ Parts 3 & 4: Efficiency Evaluation: Tests the tokenization efficiency across different vocabulary sizes (```K``` = 10, 50, 100, 200, 500, 1000) and sequence lengths (```M``` = 1,000, 5,000, 10,000, 50,000). It selects random subsequences for training and testing, repeating 100 times per parameter combination to compute the average tokenized sequence length.

# Setup and Data 
This project uses human chromosome 20 as the training and testing corpus.   
1. Download [```chr20.fa.gz```](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr20.fa.gz) from the UCSC Genome Browser.
2. Ensure the downloaded file remains compressed; the script reads the ```.gz``` archive directly using Python's built-in ```gzip``` module.

# Usage
The main script is ```Project_7.py```. Run it directly from your terminal:   
```
python Project_7.py
```
