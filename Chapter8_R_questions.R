## Study questions for Chapter 7 computer lab session ##

# Name:
# Date:

### Note that the tutorial that will explain most of the answers to these questions is
### here: https://markravinet.github.io/Chapter7.html

### INSTRUCTIONS ###
# The easiest way to submit your work is use this R script as a template. You should show
# your code. You should also show the answers using the hash (#) symbol to comment them
# out.
#
# When you are finished, you can click, File>Save_as and alter the extension of the script
# from .R to .txt - this will save it as a text file which you can then upload to Canvas.
# We recommend you do this last of all, because it will destroy the syntax highlighting in
# your R code.

## 1. How far is the distance between the final 100 Kb/25 Kb jump window and the end of chr8?

## 2. Generate a set of 10 Kb sliding windows with 5000 bp jump for the sparrow chromosome 8. How many windows does this generate?
##    How far is the final window from the end of the chromosome?

## 3. Calculate the mean pairwise Fst and also the mean pairwise dxy for all of the different species comparisons
###   N.B. if you use a tidyverse solution, it may be easier to use t() to transpose and see the final result

## 4. Following on from question 2, use ggplot and geom_boxplot to visualise the distribution in these variables
##    NB: To rotate axis labels in ggplot, add the following: 
##    + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## 5. Draw a composite plot (fst, nucleotide diversity, dxy) for the house versus italian comparison.
##    Is there a similar pattern to the house/spanish data we examined in the tutorial?
##    NB you do not need to reorder the levels for the plot order, but if you do, be aware you need
##    to play with it to get it right

## 6. Plot the relationship between dxy for house vs spanish and recombination rate. Is it the same as that for Fst?
##    Explain why/why not? You will need to consult the textbook for the second part of this question. 
