###############################################################
# Microbiome Bioinformatic Analysis                           #
# 16S amplicon sequencing data with DADA2                     #
# Based on DADA2 tutorial by Benjamin Callahan                #
# Data: Mice - Antibiotic Experiment                          #
# Author: ArrietaLab - University of Calgary                  #
# Date: July 2025                                             #
# Location: IUCBC - Córdoba (Argentina)                       #
###############################################################

# Install & Load Packages ####

# Install and load DADA2 (Callahan et al., 2016)
# First install the BiocManager package so that you can install packages from the BioConductor repository
# To see other ways of package installation please see the Introduction to R slides
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager") 
#library(BiocManager) # library() loads the installed package

# Then use the install() function from BiocManager to install DADA2 (only required once)
#BiocManager::install("dada2")  

# Now you can load DADA2 package 
library(dada2); packageVersion("dada2") #1.30.0
# packageVersion() to verify the version of package you use for the analysis
# This is essential for reproducibility when writing scientific papers.

# For any installation issues, visit: https://benjjneb.github.io/dada2/dada-installation.html

# Install and load tidyverse (Wickham et al., 2019)
# Useful for streamlined coding
#install.packages("tidyverse") 
library(tidyverse); packageVersion("tidyverse") #2.0.0

# Note:
# You only need to install a package once on our computer. 
# To use the installed package, you need to load the library every time you start a new R/RStudio environment.
# If you update R, you need to re-install all packages

# Set Up R Environment ####
# Clear objects from the workspace (global environment). 
# This ensures only the data you load into your environment is present and helps avoid any errors.
rm(list = ls(all = TRUE))  

# Download required data from https://github.com/memoll/MicrobiomeBioinformaticsAnalysisWorkshop2025

# Set working directory to the file path where your data is located
setwd("~/Documents/Argentina_bioinfoWorkshop_July2025") 
# Check if your working directory is correctly set
getwd() 

# Set path to the folder where sequencing files are located
#path <- "~/Documents/Argentina_bioinfoWorkshop_July2025/mouse_demultiplexed" 
path <- "mouse_demultiplexed" # Creates path object (you don't need to specify entire path if it's the same as your working directory)
list.files(path, pattern = "fastq") # Lists fastq sequencing files at path location that have a specific naming convention

# List and sort forward and reverse fastq files 
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE)) # Forward sequences, indicated by R1
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE)) # Reverse sequences, indicated by R2
fnFs
fnRs

# Inspect Sequence Quality ####
plotQualityProfile(fnFs[1:2]) # Forward sequences quality for the first 2 samples
plotQualityProfile(fnRs[1:2]) # Reverse sequences quality for the first 2 samples

# Repeat for all 16 or you can check all at once
plotQualityProfile(fnFs[1:16]) # Forward sequences quality for all 16 samples
plotQualityProfile(fnRs[1:16]) # Reverse sequences quality for all 16 samples

# Overall sequence quality across all samples
plotQualityProfile(fnFs, aggregate = TRUE) # Forward sequences quality aggregated in one plot
plotQualityProfile(fnRs, aggregate = TRUE) # Forward sequences quality aggregated in one plot

#*****SLIDES####

# Filter & Trim ####
# Create a list of sample names
# basename: sets the path to an object
# strsplit: splits the elements of a character vector x into substrings according to the matches to substring split within them
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names) # make sure the file list is correct

# Get more information about a function by adding "?" in front of it to open the help pages
# ?sapply

# Assign the filenames for the filtered fastq.gz files
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) # Forward filtered sequences
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) # Reverse filtered sequences
# Use samples names for the filtered files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names 

# Filter & Trim
# ?filterAndTrim #default: compressed files
# truncLen: discard reads shorter than this; (forward trim length, reverse trim length)
  # Since the quality of the reverse reads are always worse, we trim more off of these
# maxN: max number of ambiguous bases (N) allowed; DADA2 cannot handle ambiguous bases so we keep it at 0.
# maxEE: max number of estimated errors allowed by an individual read; increase ONLY if you have low quality reads (see out2)
# truncQ: truncates reads at the first instance of a Q score less than or equal to the value specified
# multithread: using multiple CPU threads to run intensive tasks in parallel; 
  # if TRUE: uses all available CPUs, integer: number of CPUs to run the task; 
  # if FALSE (default): only uses one single CPU (more processing time)
  # On Windows, set multithread = FALSE; Windows does not support multithreading

out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, truncLen = c(240,160),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE) 
head(out)


# Play around with truncation length and expected errors
#out2 <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, truncLen = c(240,155),
 #                    maxN = 0, maxEE = c(2,5), truncQ = 2, rm.phix = TRUE,
  #                   compress = TRUE, multithread = TRUE) 
#head(out2)
# In this case, the least amount of reads are lost in the first output - we will therefore use it in the subsequent steps.

# Compare sequence quality profiles of unfiltered and filtered files
plotQualityProfile(c(fnFs[1],filtFs[1])) # Compares first forward sequence raw vs. filtered
plotQualityProfile(c(fnRs[1],filtRs[1])) # Compares first reverse sequence raw vs. filtered

# Assess the percentage of reads lost with filtering and trimming
# Function to calculate the percentage of remaining reads
percentage_reads <- function(output){
  as.data.frame(output) %>% mutate(percentage = 100*reads.out/reads.in)
} 
percentage_reads(out) 

# Tells you the average percentage of reads out across all samples
mean(percentage_reads(out)$percentage)

#*****SLIDES####

# Error Rates ####
# Learn the error rates
errF <- learnErrors(filtFs, multithread = TRUE) # Creates list with forward error rates
errR <- learnErrors(filtRs, multithread = TRUE) # Creates list with reverse error rates

# Visualize estimated error rates
plotErrors(errF, nominalQ = TRUE) # Forward error rates
plotErrors(errR, nominalQ = TRUE) # Reverse error rates

#*****SLIDES####

# Sample Inference (Denoising) ####
# Determine the number of unique ASVs per sample
dadaFs <- dada(filtFs, err = errF, pool = "pseudo", multithread = TRUE) 
dadaRs <- dada(filtRs, err = errR, pool = "pseudo", multithread = TRUE) 
dadaFs[[1]] 
dadaRs[[1]]
# pool = TRUE: pools all samples together prior to sample denoising.
# pool = FALSE: denoising is performed on each sample separately. 
# pool = "pseudo": pseudo-pooling includes 2 steps:
# 1. Denoising and processing samples individually to detect prior ASVs (to increase the algorithm's sensitivity to rare variants)
# 2. Re-processing samples individually based on prior ASVs to detect rare variants 
# Recommended when the diversity among samples is high, but it doubles the processing time.

#*****SLIDES####

# Construct Sequence Table ####

# Merge denoised forward and reverse sequences 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 12, maxMismatch = 0, returnRejects = FALSE, verbose = TRUE) # Creates list of data frames for each sample
head(mergers[[1]])
# minOverlap: Min overlap length: 12 bp (default)
  # May need to modify the previous trimming to make reads overlap
# returnRejects = TRUE, doesn't discard pairs failing the mismatch criteria and included them in the output dataframe. 
  # This offers flexibility in data handling, allowing users to tailor the merging process to their specific research needs.

# Construct sequence table (ASV table) from the merged reads
seqtab <- makeSequenceTable(mergers) # matrix
dim(seqtab) # rows (first number): number of samples, columns (second number): number of unique ASVs

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) # Tells you the number of sequences at each length
# The lengths of the merged sequences should all fall within the expected range for the length of the V4 amplicon in the 16S rRNA gene 
# 515F–806R primer pair size typically ranges from 300 to 350 bp, including primer sequences, natural length variations within the V4 region, and specific primer design considerations
# Without primers, the 16S V4 region is 250-254 bp.

# Remove chimeras
# These are single sequences generated from multiple parent sequences and are considered an artifact of sequencing technologies.
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim) 

# Percentage of non-chimeric sequences
sum(seqtab.nochim)/sum(seqtab) * 100 

# Track Reads Through Pipeline ####
# How many reads are lost after each step?

# Creates a function to tell you the frequency of each sequence
getN <- function(N){
  sum(getUniques(N))
} 
# getUniques: looks at object, gets sequence name, and tells you how many times it is observed

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)) #cbind: combine by column
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") # name columns
rownames(track) <- sample.names #name rows
track
view(track) # opens track in a new window
# We can see that most of our reads were maintained from step to step, apart from removing chimeras.

#*****SLIDES####

# Download the DADA2-formatted reference database from https://benjjneb.github.io/dada2/training.html
# For future work, check for the most recent update of the 16S database here https://www.arb-silva.de/

# Assign Taxonomy ####
# This step may take a while to run, especially when working with large datasets
# It is good practice to save your environment prior to running this step to avoid losing your work if your computer crashes or freezes

save.image("~/Documents/Argentina_bioinfoWorkshop_July2025/results/antibio_dada2.RData") # Save the workspace

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/Argentina_bioinfoWorkshop_July2025/silva_nr99_v138.2_toSpecies_trainset.fa.gz", multithread = TRUE) # Assigns taxonomy to species level
# This step may take a while to run, especially when working with large datasets

colnames(taxa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") # Assign column names
unname(head(taxa)) # Hide row names (long sequences)
dim(taxa) # Dimensions of taxa table; # unique ASVs (rows), # columns (KPCOGFS)

view(taxa) # Examine taxa table

# Save your final workspace ####
save.image("~/Documents/Argentina_bioinfoWorkshop_July2025/results/antibio_dada2.RData")
