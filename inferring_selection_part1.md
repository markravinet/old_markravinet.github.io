---
title: "Inferring selection in the genome - part 1"
layout: archive
permalink: /sliding_windows/
---

There are many different ways to detect regions under divergent selection. In this tutorial, we are going to compute three of them in sliding genomic windows:
- pi, a measure of genetic variation
- Fst, a measure of genomic differentiation
- dxy, a measure of absolute divergence

Note that dxy *d*<sub>XY</sub> and pi require monomorphic sites to be present in the dataset, whereas  *F*<sub>ST</sub> is only computed on bi-allelic sites. It is thus important to filter out indels and multi-allelic sites and to keep monomorphic sites (i.e. do not apply an maf filter).

As SNP  *F*<sub>ST</sub> values are very noisy, it is also better to compute  *F*<sub>ST</sub> estimates for entire regions. Selection is expected to not only affect a single SNP and the power to detect a selective sweep is thus higher for genomic regions. How large the genomic region should be depends on the SNP density, how fast linkage disequilibrium decays, how recent the sweep is and other factors. It is thus advisable to try different window sizes. Here, we will use 100 kb windows. An alternative is to use windows of fixed number of sites instead of fixed size.

### Setting up

In order to calculate these statistics, we will use a set of custom `python` scripts written by [Simon Martin](https://simonmartinlab.org/) which you can download [here](https://github.com/simonhmartin/genomics_general). We will spend a few moments learning how to install these scripts as we will use them here and also in the practical session later today too. Thankfully, they are very easy to set up, we can download and install them using `git`. We do this like so:

```shell
# move to home directory
cd ~
# use git clone to clone the genomics_general github repo
git clone https://github.com/simonhmartin/genomics_general.git
```

This will download the scripts into a folder in your home directory. You can look inside to see all the different scripts available. However, what we really want to do is to make sure that they are accessible wherever we are on the cluster. We can do that easily too. We simply add the directory of the scripts to the `$PATH` variable like so:

```shell
# create a custom path
export PATH=$PATH:$HOME/genomics_general
```

Now if you try typing the following into the shell, you should see the help dialogue:

```shell
popgenWindows.py --help
```

However, the `$PATH` variable we created just now is temporary. If you log off, it will be gone. So we can echo it to our `.bashrc` file to make sure it is loaded each time we start a new shell session. We do that like so

```shell
echo "# create a custom path" > ~/.bashrc
echo "export PATH=$PATH:$HOME/genomics_general:$HOME/genomics_general/VCF_processing" > ~/.bashrc
```

Now the scripts are easily callable!

### Converting the input files

First, we will convert a `vcf` file into Simon Martin's `geno` file format. To do that we need to use the `parseVCF.py` script. The `geno` file is very straightforward - it is essentially a genotype matrix and easy to create.

We will convert a vcf with data from 105 individuals into a geno file that we can process downstream. First lets set up a directory to work in:

```shell
cd ~
mkdir genome_scan
cd genome_scan
```
Then we copy over our data into this directory we have just created:

```shell
cp /home/file_share/Session7/sparrows_chr8_subset.vcf.gz ./
```
Before we go further with converting the file, it is worth looking at it in a bit more detail to get some idea of what it contains:

```shell
# what individuals?
bcftools query -l sparrows_chr8_subset.vcf.gz
# how many individuals?
bcftools query -l sparrows_chr8_subset.vcf.gz | wc -l
# how many variants?
bcftools index sparrows_chr8_subset.vcf.gz
bcftools index -n sparrows_chr8_subset.vcf.gz
```

We can very easily convert the file. We do so using the `parseVCF.py` script.

```shell
parseVCF.py -i sparrows_chr8_subset.vcf.gz -o sparrows_chr8_subset.geno
```

This should take a few seconds to run and then we can have a look at the input:

```shell
head sparrows_chr8_subset.geno | cut -f 1-10 | column -t
```

Essentially this is a genotype matrix for each of our individuals. It is also a good idea to compress this file with `gzip`. This is also really straightforward.

```shell
gzip sparrows_chr8_subset.geno
```

Now we are ready to analyse the data!

### Performing the analysis

We want to calculate *d*<sub>XY</sub>, *F*<sub>ST</sub> and pi from the house sparrow and the bactrianus sparrow from this dataset. We can do this very easily, but first it will help if we create a population file that we can give to the script so that we don't need to manually input each of the sample names. We can do this fairly easily using `grep`

```shell
# generate all the sample names
bcftools query -l sparrows_chr8_subset.vcf.gz > samples
# get the bactrianus sparrow samples
grep "^[KP][Ad]" samples > bactrianus
# get the house sparrow samples
grep -v "^[KP][Ad]*" samples > house
```

However, the population file also needs a column with the population denoted in it. So we can generate this with `awk`.

```shell
awk '{print $1"\tbac"}' bactrianus > pop_file
awk '{print $1"\thouse"}' house >> pop_file
# now look at the pop_file
cat pop_file
```

With this, we are ready to perform our analyses.

```shell
popgenWindows.py -g sparrows_chr8_subset.geno.gz -o div_stats.csv \
   -f phased -w 100000 -m 10 -s 25000 \
   -p house -p bac \
   --popsFile pop_file \
   --writeFailedWindow
```
Some notes on the options here:

- `-g` the input geno file
- `-o` the output stats file
- `-f` the format of the input file - here phased simply means we have both sites (i.e. it is diploid)
- `-w` the size of the window in bp
- `-m` the minimum number of sites within a window - here set to 10 to get an output
- `-s` the size of the step for the sliding window - here it is 25 Kb.
- `-p` the population names to compare from the pops file.
- `--popsFile` the name of the pops file
- `--writeFailedWindow` write to the output even if the window does not meet filtering requirements

This will run very quickly and write out our results after printing a dialogue to the screen. We can then use `head` to examine the results. Of course it is worth remembering here that this is just a dummy dataset that we used to run the analyses quickly. In the next section, when we visualise this in `R`, we will use a real dataset!

### Visualising the results

The real analysis we will visualise in this step was conducted on 413 individuals from across the *Passer* sparrow distribution. This includes different species but here we will use a comparison of *d*<sub>XY</sub>, *F*<sub>ST</sub> and pi from house sparrows in Europe versus bactrianus sparrows from Kazahkstan and Iran.

We can get the real data like so:

```shell
cd genome_scan
cp /home/file_share/Session7/house_bac_div_stats.csv
```

Now we will turn to `R` to visualise our data. The first thing we need to do is prepare the `R` environment.

```r
# clear the r environment
rm(list = ls())
# load the tidyverse library
library(tidyverse)
```

Then with this done, we can read in the data.

```r
infile <- "./house_bac_div_stats.csv"
house_bac <- read_csv(infile)
```

We can now take a look at the `house_bac` dataframe in a bit more detail. Note that it is a `tibble` so that we can easily get a quick summary of the data. This contains estimates of pi, *d*<sub>XY</sub> and *F*<sub>ST</sub> in 100Kb windows with a 25Kb step along our chromosome of interest (chromosome 8 here).

Let's plot our statistics to get a better idea of what is going on here. We will start with *F*<sub>ST</sub>:

```r
# fst
ggplot(house_bac, aes(mid, fst_bac_house)) + geom_line()
```

Next we will look at *d*<sub>XY</sub>

```r
# dxy
ggplot(house_bac, aes(mid, dxy_bac_house)) + geom_line()
```

Clearly there is something going on at the centre of the chromosome here. We should also look at the estimates of pi to see if there is anything we can learn from that too. Since we want to plot pi from both species on the same plot, it requires a bit more leg work. First we need to use the `gather` function to reshape the data so that we can plot it.

```r
# pi
# first gather the data
pi <- house_bac %>%
  select(mid, contains("pi")) %>%
  gather(key = "species", value = "pi", -mid)
# then plot it
ggplot(pi, aes(mid, pi, colour = species)) + geom_line()
```

pi also shows a large drop in the level of diversity at the middle of the chromosome. This suggests an area of reduced recombination, most likely caused by the presence of the centromere. This can also explain the largest *F*<sub>ST</sub> peak present in the data. This suggests we should try an alternative approach.
