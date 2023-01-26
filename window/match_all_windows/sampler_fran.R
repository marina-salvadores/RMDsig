# ----- helper functions
library(readr)

# jitter should be between 0 (no jitter) and 1 (100% jitter)
generateJitteryRow = function(freqs, numNt, jitter) {
  trunc( freqs * runif(length(freqs), min=1-jitter, max=1+jitter) * numNt)
}
# normalize each row to sum to 1. There may be a faster way.
rowNorm =function(m) {
  t( apply(m, 1, function(x) { x/sum(x) } ) )
}
euclidean = function(a, b) {
  sqrt(sum((a - b)^2))
} 


# ----- generate simulated data (instead of this, in reality load data)

numNucleotides = 10000  # num nt available (starting length of each window)
defaultFreq = c(0.10, 0.10, 0.15, 0.15, 0.25, 0.25)  # this will have 32 trinucleotides in reality (or 64 if no strand symmetry)
jitter = 0.5  # pretty high jitter
simulData = t( replicate( 
  n=20, generateJitteryRow(defaultFreq, numNucleotides, jitter) 
  ) )


# ----- load data

#counts = simulData

all_counts = read_csv("results/windows/pentant_composition_original_windows_500kb.csv")
dim(all_counts)
head(all_counts)
hist(all_counts$nt_alig1)
table(all_counts$nt_alig1 < 500000)
all_counts = all_counts[all_counts$nt_alig1 >= 500*10^3,]

counts = as.matrix(all_counts[,3:ncol(all_counts)])
dim(counts)
counts[1:5,1:5]
table(is.na(counts))
#counts = counts[sample(1:nrow(counts), 1000),]

# ----- now the main sampling part!

# this is maximum tolerated distance (in relative frequency) in any column in any row.
# default 0.1% should be okay, if a bit stringent, for the 32 contexts (trinucleotide, strand-symmetrical).
stoppingCriterion = 0.001    
maxIter = 50000  # to prevent endless loops (in reality this can be a very high #, this works quite fast)
iter=0

while ( TRUE ) {

  iter=iter+1;
  freqs=rowNorm(counts);
  meanFreqs=colMeans(freqs, na.rm = T);

  # find the 'offending' row which is most different from the mean relFreq vector
  offender = which.max( apply(freqs, 1, function(x){ euclidean(x,meanFreqs) })  );

  diffs = meanFreqs - freqs[offender,];
  
  # in that row, find the column which is most responsible for the difference
  # however importantly we care ONLY about the negative differences in this vector!
  # i.e. those are the cases where the offending row has HIGHER freqs
  # (meaning we can correct that by removing sites... we can't add sites!!)
  correctableCol = which.min(diffs)
  
  worstCol = which.max(abs(diffs))  # this is sometimes the same as the correctable col

  # note that tolerance is expressed via worstCol not via correctableCol -- I am not sure if that is correct/optimal
  tolerance = abs(diffs[worstCol]) 
  
  # print some stats
  cat( sprintf( "Iteration: %7d, Euclidean: %.3f, Tolerance: %.3f, ShortestWin: %7d, AvgWin: %7d\n",
           iter, euclidean(freqs[offender,],meanFreqs), abs(diffs[worstCol]), min( rowSums(counts) ), round(mean( rowSums(counts) ))  ) )

  # calculate tolerance -- did we reduce the difference enough?
  if ( tolerance <= stoppingCriterion ) {
    cat( sprintf("Successfully completed optimization: tolerance reached. Poorest match at row %d col %d\n", offender, worstCol) );
    break;
  }
  
  # note this adjustment (subtraction) is too conservative, but by iterating it should converge to the right value
  subtractThis = round( diffs[correctableCol] * sum(counts[offender, ]) )
  
  if ( subtractThis == 0 ) {
    cat( sprintf("Successfully completed optimization: count reduction <=0.5. Poorest match at row %d col %d\n", offender, worstCol) );
    break;
  }
  
  # now simply decrease counts in the responsible column to get closer to the mean
  counts[offender, correctableCol] = counts[offender, correctableCol] + subtractThis
  
  if (counts[offender, correctableCol]<=0) {
    cat( sprintf("Cannot continue optimization - counts exhausted at row %d col %d\n", offender, correctableCol) );
    # note that I think we could continue... just mark this cell in the matrix as uncorrectable and move to another? not sure.
    break;
  }
  if (iter==maxIter) {
    cat( sprintf("Stopping optimization - maximum number of iterations reached.\n") );
    break;
  }  
  
} # next iteration  

new = cbind(all_counts[,1:2], counts)
new[1:5,1:5]
nts_after_matching = rowSums(new[,3:ncol(new)])
hist(nts_after_matching, breaks = 1000)
table(nts_after_matching >= 500000)
dim(new)
new = new[nts_after_matching>=500000,]
write.csv(new, "results/windows/pentant_composition_original_windows_imperfect_matching_maxtol0.004_50k.csv", row.names = F)



# print initial and final data
{ 
  cat("\n---- Initial data: ----\nCounts:\n");
  print(simulData);
  cat("\nRelFreqs:\n");
  print(rowNorm(simulData) );
  cat("\nMeanRelFreqs:\n");
  print(colMeans(rowNorm(simulData)))
  
  
  cat("\n---- Final data: ----\nCounts:\n");
  print(counts);
  cat("\nRelFreqs:\n");
  print(freqs);
  cat("\nMeanRelFreqs:\n");
  print(meanFreqs)
}

