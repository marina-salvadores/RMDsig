# RMD signatures

## 1) window
We split the genome into 1 megabased windows. We calculate the trinucleotide composition of each 1-Mb window and applied a matching algorithm to reduce the differences between windows. This imperfect matching tells us how many tri-nucleotides contexts in each window we can have to reach a centain tolerance. We applied this filter (remove trinucleotide positions) to our 1 megabased windows to make them more homogeneous and allow comparison across them.
- input: genomic coordinates of the chrs, regions of the genome to remove (black list regions, exonic regions, etc).
- output: 1 matrix per window with the genomic coordinates of the specific regions selected to match the trinucleotide composition of every window.

## 2) RMD counts
We parse the somatic mutations (WGS) files and calculate the number of mutations in each of the 1-Mb windows.
- input: 1-Mb window coordinates (with the matching), somatic mutations for human tumors
- output: 1 matrix with the counts of mutations in each sample per each 1-Mb window

## 3) NMF
We applied an NMF to the matrix of mutations counts per samples and per windows. Each sample counts are normalized separately by chromosome arm to account for possible arm copy number amplifications.
We first applied a bootstrapping to generate a similar matrix with a few of noise. We performed the bootstrapping 100 times.
We run NMF into the 100 bootstrapped matrices for 1:30 nfactors 
We applied clustering for each factor separately for the 100 runs
We select the number of clusters and factors that is better according to the clustering silohuette index meassurement
We obtained a final subset of RMD signatures
- input: RMD counts matrices
- output: RMD signatures for the selected k and factor (in our case k=nfactor=13)

## 4) Simulated data
We performed a test to see if we can capture the variability we are interestsed in with our NMF method.
### 4.1) generated data
We generated 9 signatures (different number of windows and different strength)
We applied the signatures to a subset of samples in a random manner. We do for 9 different conditions varying the number of samples affected, the number of mutations that the signature contribute to the total mutation burden
### 4.2) NMF
we applied NMF in the same way as described in point 3 for the simulated data
