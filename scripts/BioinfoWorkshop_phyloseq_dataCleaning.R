###############################################################
# Microbiome Bioinformatic Analysis                           #
# Introduction to Phyloseq - Data Cleaning                    #
# 16S amplicon sequencing                                     #  
# Data: Mice - Antibiotic experience                          #
# Author: ArrietaLab - University of Calgary                  #
# Date: July 2025                                             #
# Location: IUCBC - Cordoba (Argentina)                       #
###############################################################

# Package installation and set up #### -----------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager") # a package that facilitates the download of packages from the BioConductor repository.

#BiocManager::install("phyloseq") #(McMurdie & Holmes, 2013)
##install.packages("tidyverse") #(Wickham et al., 2019)
#install.packages("dyplr") #(Wickham et al., 2014)

# Next we need to load these packages to be able to use their functions in this session. 
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")

# Extra functions for phyloseq (including ggrare for rarefaction curves)
scripts <- c("graphical_methods.R",
             "tree_methods.R",
             "plot_merged_trees.R",
             "specificity_methods.R",
             "ternary_plot.R",
             "richness.R",
             "edgePCA.R",
             "copy_number_correction.R",
             "import_frogs.R",
             "prevalence.R",
             "compute_niche.R")
urls <- paste0("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/", scripts)

for (url in urls) {
  source(url)
}

# Set the working directory to the appropriate folder (where the input files are stored). 
setwd("~/Documents/Postdoc/Argentina Workshop")


# Load and clean data #### -----------------------------------------------------
load("antibio_dada2.RData") # loads both the taxonomy and ASV tables that were generated in the DADA2 session

# Taxonomy table
head(taxa); class(taxa) # needs to be in a matrix format

# ASV table
head(seqtab.nochim); class(seqtab.nochim) # can be in a matrix or dataframe format

View(taxa)
View(seqtab.nochim)

# Load the metadata table
metadata <- read.csv2("mouse_metadata.csv", sep = ",")  # make sure this csv file is saved in your working directory
rownames(metadata) <- metadata$Sample_ID # Need to make sure the rownames match the otu table

head(metadata); class(metadata) # Needs to be in a dataframe format
View(metadata)

# Take a look at the dimensions of all 3 datasets (ASV, taxonomy, and metadata tables)
dim(taxa); dim(seqtab.nochim); dim(metadata)

# To make the phyloseq object, we need to make sure that colnames and rownames match
identical(colnames(seqtab.nochim),rownames(taxa)) 
identical(sort(rownames(seqtab.nochim)),sort(rownames(metadata))) # Use sort to order names 
setdiff(rownames(seqtab.nochim),rownames(metadata))


# Create and explore the Phyloseq object ####-----------------------------------
# Convert into phyloseq formats
otuTable  <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
taxTable  <- tax_table(taxa)
metaTable <- sample_data(metadata)

# Create a phyloseq object
ps <- phyloseq(otuTable, taxTable, metaTable)
ps # Examine object

# Rather than looking at the entire sequences, let's rename by numbering the ASVs (e.g., ASV1, ASV2, ASV3, ...)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps))) # replace sequences w/ ASV___

# Double check that the renaming of the ASVs looks correct:
taxa_names(ps) # names of each of the ASVs

# Explore the phyloseq object and package functions ####
# Here we will focus on the functions that are relevant for the downstream analysis in this tutorial. 
# For additional functions please see phyloseq documentation. 

# There are also a few functions for editing the phyloseq object that will be used throughout the tutorial. 
# Here is a basic description of the functions to give some context. 
# Please see documentation for argument details and additional functions. 
?prune_taxa # for filtering out unwanted ASVs in phyloseq object
?prune_samples # for filtering out unwanted samples in phyloseq object
?filter_taxa # filter ASVs based on across-sample abundance with a specified function
?subset_samples # subset ASVs based on the provided conditions
?subset_samples # subset samples based on the provided conditions
?tax_glom # merges ASVs at the indicated taxonomy level
?transform_sample_counts # transform ASV abundance counts in each sample with a function of choice
?psmelt # melt the phyloseq object to a dataframe 

# Take a look at the top section of each phyloseq component:
head(otu_table(ps)) # for both constructing and accessing the table ASV abundance 
head(tax_table(ps)) # for both constructing and accessing the table of taxonomic names
head(sample_data(ps)) # for both constructing and accessing the table of sample-level variables

  # We can also use "View" to look at the entire tables. For example:
  View(otu_table(ps))

# Take a look at the different taxonomy and variable names
rank_names(ps) # gives the taxonomic ranks included in the taxa table - useful to know available ranks to analyze
variable.names(sample_data(ps)) #names of the variables (in the metadata table)

  # We can also "colnames" to look at these same variables. For example:
  colnames(sample_data(ps))

# We can also extract the components from the phyloseq object, if needed. 
taxa <- as.matrix(tax_table(ps)); head(taxa) # extract the taxonomy table as a matrix
comm <- as.matrix(otu_table(ps)); head(comm) # extract the ASV table as a matrix
meta <- data.frame(sample_data(ps)); head(meta) # extract the sample metadata as a data frame


# Let's look at summarizing some of the aspects of our phyloseq object:
ntaxa(ps) # number of unique taxa across the entire dataset
nsamples(ps) # number of samples

taxa_sums(ps) # Total reads per ASV summed across all samples
summary(taxa_sums(ps)) 
hist(taxa_sums(ps))
hist(log10(taxa_sums(ps))) # Distribution of sequences

sample_sums(ps) # Total number of sequences (sequencing depth) per sample 
summary(sample_sums(ps)) 
hist(sample_sums(ps)) # distribution of sequences in samples


# Let's save these values of sequencing depth per sample to the metadata for reference. 
# NOTE: @ for S4 objects (like phyloseq) is similar to $ for data frames in that it specifies the data component to be accessed.  
ps@sam_data

# We'll create a new variable "depth" in the metadata based on the sample_sums
ps@sam_data$depth <- sample_sums(ps) # create a new column in the sample data for the total number of ASVs observed for each sample (sequencing depth)
ps@sam_data # confirm that the new variable has been added
summary(ps@sam_data$depth)


# SLIDES ####

# Pre-processing and denoising #### --------------------------------------------
# If you have sequenced your negative and positive samples, which is strongly recommended, here you can also verify them
# To evaluate the accuracy of sequencing experiments:
# 1. Check if the ASVs in positive samples are similar to your initial mock samples. Then, you can remove the positive samples from your dataset
# 2. Ensure negative samples have close to zero sequences (or significantly low number of reads compared to the other samples). 
# The reads in negative samples are potential contaminants that you can use later for decontaminating your dataset.

# 1. Remove common contaminants and undefined kingdoms

# Identify common contaminants in our dataset
subset_taxa(ps, Kingdom == "Archaea")
subset_taxa(ps, Phylum == "Cyanobacteria") # If you are analyzing soil or phyllosphere microbiome, depending on your research question, you might need to remove Cyanobacteria and/or Chloroplast 
subset_taxa(ps, Order == "Chloroplast")
subset_taxa(ps, Family == "Mitochondria") # these are often host sequences that need to be removed

  # We can also look at the overall taxonomy table and search these contaminants
  View(tax_table(ps)) # use the search box at the top and type in each contaminant name

  # Lastly, you can also look at all of the entries for a given taxonomic level to identify possible contaminants.
  # For example:
  unique(tax_table(ps)[, "Order"]) # Chloroplast is identified in the list

# Remove identified common contaminants
ps1 <- subset_taxa(ps, Order != "Chloroplast")
ps1 <- subset_taxa(ps1, Family != "Mitochondria") # notice how we changed which phyloseq object we're subsetting from


# Remove unidentified kingdoms and phyla
subset_taxa(ps1, is.na(Kingdom))
subset_taxa(ps1, is.na(Phylum))
# Question: How many taxa are unidentified at the Kingdom and Phylum levels?

  # We can also look to see if there are NAs by listing all of the entries at a given taxonomic level (as we did above):
  # For example:
  unique(tax_table(ps1)[, "Phylum"])

ps1 <- subset_taxa(ps1, !is.na(Phylum)) # note: !is.na means "is not NA"
# In our case, the same ASV is unidentified at both the Kingdom and Phylum level 
# so we don't need to repeat this line of code for each taxonomy level

100-(ntaxa(ps1)/ntaxa(ps))*100 # percentage of filtered taxa by the above threshold 
ps1 # how many ASVs are now in the dataset? 


# 2. Remove contaminants from negative controls ####
# Use the decontam package
# library(decontam); packageVersion("decontam"): https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

  # In this case, we don't have positive or negative controls in our dataset so we will
  # skip this step. However, it is recommended to include these controls when sequencing
  # your own samples. 

  # Once the negative controls have been used to remove identified contaminants, 
  # they can be removed from the phyloseq object. It is helpful to have a column in your
  # metadata to identify which sample is a true sample and which is a control. For example:
  # ps1 <- prune_samples(sample_data(ps1)$SampleType != "Control", ps1) 
  # with "SampleType" being the variable name to differentiate between samples and controls


# 3. Remove singleton and doubleton ASVs ####
ps2 <- prune_taxa(taxa_sums(ps1) > 2, ps1)
# Singleton and doubleton ASVs are ASVs that only occur once or twice in the dataset, respectively
# Question: How many ASVs were singleton or doubletons?

100-(ntaxa(ps2)/ntaxa(ps1))*100 # percentage of filtered taxa by the above threshold 

 
# 4. Determine the threshold of sequencing reads per sample to keep ####
# First, let's look at the sequencing reads in increments of 1%
quantile(sample_sums(ps2), seq(0,1,0.01))
# For example, this means that in 10% of samples, there are >21,000 sequencing reads

max(sample_sums(ps2)) # maximum number of sequencing reads per sample
min(sample_sums(ps2)) # minimum number of sequencing reads per sample

# Now, let's look at the rarefaction curves
ggrare(ps2, step=100, se=FALSE, color="Sample_ID") 
# You can adjust the step depending on how deep your sequencing is. For example, if you have
# millions of reads per sample, you can increase the step; otherwise, it will take a long time
# to generate the rarefaction curves. 

# In this case, as we increase the sample size (x-axis), we don't get any more new ASVs 
# (species richness; y-axis) after ~10,000 reads, with the exception of one sample.
# Since the minimum sequencing reads in a sample was >18,000 (above the plateau point),
# let's keep all of our samples. However, if we had samples with sequencing depths below
# the plateau point, we may want to omit those samples. 

# For example:
# ps_filtered <- subset_samples(ps2, depth >= x) # with x being your plateau point

# Some researchers may also want to rarefy their data which can be done using the following code:
# ps_rarefied <- rarefy_even_depth(ps2, rngseed = 1923, sample.size = x) 
# with x being your plateau cutoff point, and rngseed being a random number (any number will do!)


# 5. Finalize our cleaned phyloseq object
ps_clean <- ps2; ps_clean 

# Checking for the number of unique phyla and genera across all samples
unique(tax_table(ps_clean)[, "Phylum"]) # 7 present
unique(tax_table(ps_clean)[, "Genus"]) # How many genera are present? Are there NAs?

# Save denoised and decontaminated phyloseq object
saveRDS(ps_clean,"ps.rds")


# Explore the denoised and decontaminated phyloseq object #### -----------------
# We can now calculate the new sequencing depth and compare this to the original to determine the reads lost from filtering. 
ps_clean@sam_data$depth_filt <- sample_sums(ps_clean) # save number of ASVs per sample after filtering to new column
metadata_clean <- data.frame(sample_data(ps_clean)) # save sample dataframe with new and old sequencing depth columns to perform stats below
View(metadata_clean)

# Calculate % of reads lost after filtering
sdf_read_lost <- metadata_clean %>%
  mutate(read_lost = (depth - depth_filt)) %>% # new column with difference in reads pre- and post-filtering
  mutate(percent_lost = 100*(read_lost/depth)) %>% # convert reads lost to percentage for easier interpretation
  summarise(mean_lost=mean(read_lost), sd = sd(read_lost), mean_percent_lost=mean(percent_lost)) # create dataframe with stats
sdf_read_lost


# Save the workspace ####
save.image("antibio_ps.RData")

