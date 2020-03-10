# Inferring selection in the genome - part 2

### Going further - using haplotype statistics

One alternative to using measures such as *F*<sub>ST</sub>  is to use a haplotype homozygosity statistic, as these are robust to confounding factors such as variation in recombination rate. To calculate that, we need to use **statistically phased data**.

Haplotypes are best assessed with phased dataset. Phasing basically means figuring out for heterozygous positions which of the alleles are part of the same haplotype or chromosome (e.g. the chromosome inherited from the mother). For instance if an individual is heterozygous at two SNPs with genotypes `AG` and `TC`, phasing would tell us if the allele A at SNP1 one was inherited on the same chromosome like `T` or like `C` at the second SNP. Phased genotypes would be given as `A|G` and `C|T` meaning that A and C are on the same chromosome (e.g. maternal) and `G` and `T` are on the same chromosome (e.g. paternal). In the absence of long reads that span these SNPs, we can use statistical phasing using all sequenced individuals of the same population (the more the better). There are lots of different tools for phasing and most of them also impute missing genotypes. This means that they infer missing genotypes statistically resulting in a dataset without missing data. If you have a linkage map, it is recommended to use it for making phasing more accurate. However, here we will just perform a very basic phasing with [SHAPEIT2](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gettingstarted) which does not impute genotypes and does not require recombination rate information.

#### Phasing the dataset

Before we start phasing, we will create the directory we will be working in.

```shell
# move into home directory
cd ~
# make phasing directory
mkdir phasing
# move inside
cd phasing
```

Now we are ready to download and look at the data.

```shell
cp /home/file_share/Session7/sparrows_chr8_variants_subset_unphased.vcf.gz* .
```

Next we will make an environmental variable to make our downstream analysis more straightforward.

```shell
VCF=~/phasing/sparrows_chr8_variants_subset_unphased.vcf.gz
```

Using environmental variables is a good way of running your analyses in Unix as it lets you easily rerun scripts and generalise them. We can also use Unix string manipulation to make it easier to write outfiles. For instance:

```shell
# remove the suffix
echo ${VCF%_u*}
# add the phased suffix
echo ${VCF%_u*}_phased
# now we're happy with this, assign a new variable
OUTPUT=${VCF%_u*}_phased
```

With this taken care of, we are ready to run `shapeit`. Note that here we are running the program only for a single chromosome - this is the most practical (and sensible) way to phase - it makes it much easier to parallelise.

We do this like so:

```shell
shapeit --input-vcf $VCF \
 -O $OUTPUT \
--window 0.5 -T 4 --force
```

What did we do with this command?

* `--input-vcf` - this flag allows us to read in the vcf we want to phase.
* `--window` - sets the window size within which the phasing is carried out. The default is 1 Mb, here we set it to 0.5 Mb.
* `-T` - set the number of threads to use for parallelisation - we used 4 here

Next we run the analysis. It will take a few minutes and will print updates to the screen as it runs. Take a short break while it runs... and when it is done we can have a look at the results!

#### Examining the phased data and converting it to a phased vcf

Once `shapeit` has finished running we can see it has produced two files, `sparrows_chr8_phased.haps` and `sparrows_chr8_phased.samples`. We can look at these in more detail:

```shell
head sparrows_chr8_variants_subset_phased.sample
```

The `sparrows_chr8_phased.samples` file is simply a list of samples with the ID for each sample repeated in two columns and a third column showing the proportion of missing data. Since there is no missing data in this dataset, this column is just full of zeros.

We can also look at the `.haps` output.

```shell
less sparrows_chr8_variants_subset_phased.haps
```

This is basically a matrix with the first three 5 columns identical to those in a vcf - i.e. chromosome, ID, position, reference allele, alternative allele. After this, each entry is the phased allele for each individual, where `0` is the reference allele and `1` is the alternative.

If we want to use this data downstream, particularly in the `R` package `rehh`, we need to convert it to a vcf. Luckily this is easy using `shapeit`:

```shell
shapeit -convert \
--input-haps ${OUTPUT} \
--output-vcf ${OUTPUT}.vcf
```

We can then compress and index the vcf, like so:

```shell
bgzip ${OUTPUT}.vcf
bcftools index ${OUTPUT}.vcf.gz
```

Before moving on to the next step, let's have a quick look at our phased vcf to see how it is different from a standard vcf.

```shell
bcftools view -H ${VCF} | head | cut -f 1-12
bcftools view -H ${OUTPUT}.vcf.gz | head | cut -f 1-12
```

Comparing the two, we can see tht the phased vcf only contains the genotypes - all the other information has been stripped out. Furthermore, in the phased vcf, genotypes are encoded as `0|0`, `0|1` or `1|1` instead of `0/0`, `0/1` or `1/1`. This is because in a vcf, `|` is typically used to denote that the phase of these loci are known. Thus, it is possible to read the haplotype of an individual by reading downwards across loci on either side of the `|`.

### Subsetting the vcf for a selection scans

We will be using our phased vcf for long-range haplotype statistic estimation in `rehh`. While it is possible to split the populations in the vcf apart in R, it is a bit clumbsy to do so. instead, it is easier to split the vcf using `bcftools`. To do this, we first need the sample names

```shell
# set a new variable
VCF=~/phasing/sparrows_chr8_phased.vcf.gz
# look at the sample names
bcftools query -l $VCF
# extract sample names
bcftools query -l $VCF > samples
```

As we learned earlier, this vcf contains data from two house sparrow subspecies, the European house sparrow (from France and Norway) and the wild Bactrianius sparrow (from Kazahkstan and Iran). We split the population data like so:

```shell
grep "^[8F]" samples > house
grep -v "^[8F]" samples > bac
```

Note we used the `-v` flag to do an inverse `grep` - i.e. extract the opposite of this pattern. Next we split the vcf:

```shell
bcftools view -S house -O z -o house_chr8.vcf.gz $VCF
bcftools view -S bac -O z -o bac_chr8.vcf.gz $VCF
```

What did we do here? We used the `bcftools index` command to extract the samples for each population. The `-S` flag extracts the samples listed in each file. `-O z` specifies that we want a compressed vcf. Finally `-o` tells `bcftools` where to write the output. Now all we need to do is index the vcfs.

```shell
bcftools index house_chr8.vcf.gz
bcftools index bac_chr8.vcf.gz
```

Now we're ready to read this into `R` for a selection scan analysis.

#### Getting access to a larger dataset

Remember that above, we performed our phasing on a drastically reduced dataset in order to complete our analysis in time. In order to get the real dataset, we should do the following.

```shell
cd ~
mkdir rehh
cd rehh
cp /home/file_share/Session7/[hb]*.vcf.gz .
```

These are the full vcfs for the house and bactrianus datasets that we can read into `R` to perform our genome scan with.

### Running an analysis with rehh

To compute and detect regions with extended haplotype lengths, we will use the excellent `rehh` R package. For more information see the [vignette](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html). `rehh` is well maintained, continually updated and has [a very informative tutorial](https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html) which I recommend you also check out.

To run `rehh` and perform our analyses, we need to run things in `R`. You can either download the data from [our github](https://github.com/speciationgenomics) and run it locally on your own machine, or we can use our `RStudio` server. Once `R` is running, we are ready to go!

#### Setting up the R environment

The first thing we need to do is clear our `R` environment and load the packages we need. Like so:

```r
# clear environment
rm(list = ls())
# load packages
library(rehh)
library(tidyverse)
```

### Reading in data from vcfs

A nice new addition to `rehh` is the ability to read in (and also filter) your data from a vcf. However, it is still quite tricky to split up individuals so instead we will read in a vcf for each population. We read in data using the `data2haplohh` function:

```r
# read in data for each species
# house
house_hh <- data2haplohh(hap_file = "./rehh/house_chr8.vcf.gz",
                   polarize_vcf = FALSE)
# bactrianus
bac_hh <- data2haplohh(hap_file = "./rehh/bac_chr8.vcf.gz",
                         polarize_vcf = FALSE)
```

This will take a few moments but once the commands are run, we will have read in our data. It is really important that we set `polarize_vcf` to `FALSE` because we have not used an outgroup genome to set our alleles as derived or ancestral. **NB. Some plots and functions in `rehh` will still refer to ancestral and derived even if you use this option**. Instead `rehh` will use the minor and major alleles (in terms of frequency) from our data.

Next, we will filter our data on a **minor allele frequency** or **MAF**. This is really simple in `rehh` with the `subset` function:

```r
# filter on MAF - here 0.05
house_hh_f <- subset(house_hh, min_maf = 0.05)
bac_hh_f <- subset(bac_hh, min_maf = 0.05)
```

This will remove some sites - and we'll be ready to run our haplotype scans.

### Performing a haplotype genome scan - *xpEHH*

Before we can calculate the statistic we are interested in *xpEHH* - we need to calculate *iES* statistics. Luckily this is really easy using the `rehh` function `scan_hh`.

```r
# perform scans
house_scan <- scan_hh(house_hh_f, polarized = FALSE)
bac_scan <- scan_hh(bac_hh_f, polarized = FALSE)
```

Next we will calculate *xpEHH* which is the cross-population *EHH* test. This is essentially a test for the probability that if we randomly sampled haplotypes from different populations, we would get different haplotypes. Again, `rehh` makes this simple with the `ies2xpehh` function.

```r
# perform xp-ehh
house_bac <- ies2xpehh(bac_scan, house_scan,
                       popname1 = "bactrianus", popname2 = "house",
                       include_freq = T)
```

Here we provide the names of our previous *iES* scans (`bac_scan` and `house_scan`). We can also provide the function with the names of our populations and finally, if we set `include_freq` to `TRUE`, we get the frequencies of alleles in our output, which might be useful if we want to see how selection is acting on a particular position.

Next, we can plot the *xpEHH* values, like so:

```r
# plot
ggplot(house_bac, aes(POSITION, XPEHH_bactrianus_house)) + geom_point()
```

In this plot, highly negative values suggest selection in population 2 (house in this case) whereas positive values indicate selection in population 1. Alternatively, like with *iHS*, we could plot the log *P* values.

```r
ggplot(house_bac, aes(POSITION, LOGPVALUE)) + geom_point()
```

### Examining haplotype structure around a target of selection

One other nice feature of `rehh` is that we can examine haplotype structure around SNPs we think might be under selection. Before we do that, we need to identify the SNP in our dataset with the strongest evidence of being an *xpEHH* outlier.

```r
# find the highest hit
hit <- house_bac %>% arrange(desc(LOGPVALUE)) %>% top_n(1)
# get SNP position
x <- hit$position
```

Here we also set the position of our putative selection SNP as the object `x`. This is because we need to identify where it occurs in our haplotype objects - unfortunately we cannot use the position for this. In the code below, we find the marker id for both our datasets.

```r
marker_id_h <- which(house_hh_f@positions == x)
marker_id_b <- which(bac_hh_f@positions == x)
```

Now we are ready to plot the bifurcation of haplotypes around our site of selection. We do this like so:

```r
house_furcation <- calc_furcation(house_hh_f, mrk = marker_id_h)
bac_furcation <- calc_furcation(bac_hh_f, mrk = marker_id_b)
```

We can also plot both of these to have a look at them:

```r
plot(house_furcation, xlim = c(19.18E+6, 19.22E+6))
plot(bac_furcation, xlim = c(19.18E+6, 19.22E+6))
```

Calculating the furcation pattern also makes it possible to calculate the haplotype length around our signature of selection.

```r
house_haplen <- calc_haplen(house_furcation)
bac_haplen <- calc_haplen(bac_furcation)
```

With the haplotype length calculated, we can now plot this to see how haplotype structure differs between our two populations.

```r
plot(house_haplen)
plot(bac_haplen)
```

Here we can see the blue haplotype is much larger around this target and is also more numerous in the European house sparrow.

### Writing out the data for later used

Finally, before we move on to the last tutorial, we are going to write out the data. We'll also make the column names all smaller letters, to make downstream scripting a bit easier.

Use the following code to achieve this:

```r
# write out house bactrianus xpEHH
house_bac <- tbl_df(house_bac)
colnames(house_bac) <- tolower(colnames(house_bac))
write_tsv(house_bac, "./house_bac_xpEHH.tsv")
```

In the last tutorial, we'll use `R` to identify genes that are close to our outlier SNPs.
