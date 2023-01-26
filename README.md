# RMD signatures

## 1) window
We split the genome into 1 megabased windows. 
We calculate the trinucleotide composition of each 1-Mb window and applied a matching algorithm to reduce the differences between windows. This imperfect matching tells us how many tri-nucleotides contexts in each window we can have to reach a centain tolerance. We applied this filter (remove trinucleotide positions) to our 1 megabased windows to make them more homogeneous and allow comparison across them.

## 2) RMD counts
We parse the somatic mutations (WGS) files and calculate the number of mutations in each of the 1-Mb windows.

## 3) NMF
### 3.1) bootstraping
### 3.2
