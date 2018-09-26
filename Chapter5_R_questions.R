## Study questions for Chapter 5 computer lab session ##

# Name:
# Date:

### Note that the tutorial that will explain most of the answers to these questions is
### here: https://markravinet.github.io/Chapter5.html

### INSTRUCTIONS ###
# The easiest way to submit your work is use this R script as a template. You should show
# your code. You should also show the answers using the hash (#) symbol to comment them
# out.
#
# When you are finished, you can click, File>Save_as and alter the extension of the script
# from .R to .txt - this will save it as a text file which you can then upload to Canvas.
# We recommend you do this last of all, because it will destroy the syntax highlighting in
# your R code.

### Discussion questions

### 1. Imagine two populations. In population one, the A1 allele has the highest marginal fitness
###    whereas in population 2, the A2 allele has the highest marginal fitness. Directional selection
###    occurs in both populations over 100 generations. What happens to the frequency of the alleles in
###    both populations? How would this alter Fst between them?

### 2. Imagine the same scenario, but this time balancing selection is occurring in both populations
###    and it maintains the alleles at equal frequencies in both populations. How would this affect
###    Fst?

### R-based questions

### 3. We sample two fish populations - one in the lake and the other in a stream. We genotype them
###    at locus B. In the lake, the genotype counts are - B1B1 = 32, B1B2 = 12 and B2B2 = 6. In the
###    stream, the genotype counts are B1B1 =10, B1B2 = 16, B2B2 = 43.

### a . Calculate the average expected heterozygosity for the two populations.

### b. Calculate the expected heterozygosity for the two populations as a metapopulation.

### c. Calculate Fst between the lake and stream fish.

### 4. Using the calc_af and calc_fst functions we developed during the tutorial and the lct_freq
###    data, calculate Fst between the Han_China and the Swedish_and_Finnish_Scandinavia populations.
###    What might explain the Fst value you calculate?

### 5. Using the functions we developed in the tutorial, calculate Fst around the LCT gene between
###    Americans of European descent and also between African Americans. Plot this like we plotted the
###    the Fst between European Americans and East Asians. What is the highest value of Fst?

### NB as a hint here, allele frequencies for European Americans and African Americans are in the 4th and
### 5th columns...
