## Study questions for Chapter 10 computer lab session ##

# Name:
# Date:

### Note that the tutorial that will explain most of the answers to these questions is
### here: https://markravinet.github.io/Chapter10.html

### INSTRUCTIONS ###
# The easiest way to submit your work is use this R script as a template. You should show
# your code. You should also show the answers using the hash (#) symbol to comment them
# out.
#
# When you are finished, you can click, File>Save_as and alter the extension of the script
# from .R to .txt - this will save it as a text file which you can then upload to Canvas.
# We recommend you do this last of all, because it will destroy the syntax highlighting in
# your R code.

# 1. Write a pseudo-code example of a for loop or an if-else statement. You can let your
##   imagination run wild! 

# 2. Below is some code to generate a single random number. Write an if and else statement to
##   evaluate whether it is greater than or equal to 5. You can have the statement write out whatever you
###  wish, as long as it is clear what the answer is!

rnumber <- rpois(1, lambda = 9)

# 3. Using ifelse and the starwars data, group the data in two ways. First do it based on homeworlds
### - splitting everyone who is born on Naboo from all other planets. Then do the same based on height
###  splitting the data into those with a height under 100 and those with a height over 100.
### Can you work out a way to count the numbers in these groups? 

# 4. Write a for loop that iterates over 100 numeric elements. For each iteration of the loop,
## Use your own name as named-variable within the curly brackets. Raise the number to the power
## of 10 at each iteration and then print it to the screen.

# 5. Create a for loop that prints the output of two different character vectors at once. Each vector
##  should be 5 elements long. They can be whatever you want.

# 6. Use a for loop to calculate the means of the numeric columns in the iris data

# 7. Write a for loop that will create a vector of 150,000 elements. Measure how long it
## takes and compare it to the time it takes to use seq and 1:150000


# 8. Use sapply to raise each value of 1:150 to the power of 3 and then multiply
##   by itself

# 9. Extend the greeting function so that it takes two names as arguments and greets
### both together. Call it greeting2. It must take to separate arguments - not one. 

# 10. Write a function that raises a number to the power of 10, divides it by 1 and then
###   adds it to the original number. Test this function on a single value and then a vector.
###   Finally, use apply to run the function on columns of the iris data 
###   (excluding the species name column)





