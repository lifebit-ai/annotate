#!/usr/bin/env Rscript

### Description


############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | "help" %in% args | (length(args) == 0) | (length(args) == 1) ) {
  cat("
      The helper R Script prep_hwe.R
      Mandatory arguments:
         --predicted_ancestries=path  - Path to the predicted ancestries file, produced by infer_ancestry
				        process of Ancestry and relatedness pipeline.

         --unrelated_list=path        - Path to unrelated list file, produced by king_coefficients process
                                        of Ancestry and relatedness pipeline. In the pipeline it is called
					autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.in.id

      Optionnal arguments:
          --help                      - you are reading it

      Usage:
      
          The typical command for running the script is as follows:
    
          ./prep_hwe.R --predicted_ancestries='predicted_ancestries.tsv' \\
		       --unrelated_list='autosomes_LD_pruned_1kgp3Intersect_triangle_HWE1_5.king.cutoff.in.id'
     
     Output:
	  The script will generate four .keep files that contain unrelated participants for four populations:
	      AFR_unrelated_pop.keep
              EUR_unrelated_pop.keep
              SAS_unrelated_pop.keep
              EAS_unrelated_pop.keep
  \n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs    <- function(x) strsplit(sub("^--", "", x), "=")

argsL        <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args         <- argsL
rm(argsL)

############################## LIBRARIES SECTION #############################

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))

# ######################### VARIABLES REASSIGNMENT SECTION ###############################

# Facilitates testing and protects from wh-spaces, irregular chars

# required
predicted_ancestries      <- args$predicted_ancestries
unrelated_list            <- args$unrelated_list

cat("\n")
cat("ARGUMENTS SUMMARY")
cat("\n")
cat("predicted_ancestries : ", predicted_ancestries,   "\n",sep="")
cat("unrelated_list       : ", unrelated_list,         "\n",sep="")


# ############################### SCRIPT SECTION ###############################

dat <- fread(predicted_ancestries) %>% as_tibble()
unrel <- fread(unrelated_list) %>% as_tibble()

# Subset the ancestries table to retain only participants that have been proven as unrelated.
dat <- dat %>% filter(Sample %in% unrel$IID)

# For each population produce a .keep file with those individuals that pass population ancestry threshold.
# We need to keep only the platekeys and not values. Moreover, the .keep file needs the platekeys column
# twice, without a header, so we just keep Sample column for each population and duplicate it.

for(col in c("AFR","EUR","SAS","EAS")){dat[(dat[col]>0.8)[,1],c("Sample")] %>%
	mutate(duplicate_col=Sample) %>%
	write.table(paste0(col,"_unrelated_pop.keep"), quote = F, row.names=F, col.names=F, sep = "\t")}

