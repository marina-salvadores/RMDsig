# RMD signatures

## 1) window
We split the genome into 1 megabased windows. 
We calculate the trinucleotide composition of each 1-Mb window and applied a matching algorithm to reduce the differences between windows. This imperfect matching tells us how many tri-nucleotides contexts in each window we can have to reach a centain tolerance. We applied this filter (remove trinucleotide positions) to our 1 megabased windows to make them more homogeneous and allow comparison across them.

## 2) RMD counts
