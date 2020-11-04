#!/usr/bin/env Rscript

### Description
# This script is used in make_header process to generate new headers for annotated bcf files.
# Current version is making a header used for autosomes only.


library(data.table)
library(dplyr)
library(magrittr)
library(stringr);


#Read in the header options
toAdd <- fread("all_seen_flags.txt", header = F) %>% as_tibble() %>% 
	filter(V1 != 'PASS');
        starter <- '##FILTER=<ID=';
        mid <- ',Description="';
        end <- '.">'


#Construct descriptions. We won't be all that descriptive for the filter field
#more so for the info field. Just state failures for filter
#Probably should add what the failing values are. Maybe only include these for the single instance onces
#to save on space
singlecases <- setNames(
	as.list(c(
	  'Fails depth, median depth < 10',
          'Fails GQ, median GQ <15',
          'Fails missingness, missingness > 0.05',
          'Fails completeGTRatio, compelete sites <0.5',
          'Fails AB ratio, AB ratio of (het sites passing binomial distribution with p-value < 0.01 / all het sites) < 0.25',
          'Fails phwe_eur, site is out of HWE (<10e-6) for inferred unrelated inferred eur superpopulation'
        )),
        c('depth',
          'GQ',
          'missingness',
          'completeGTRatio',
          'ABratio',
          'phwe_eur'
	)
   )

singlecases <- paste0(starter, names(singlecases), mid, singlecases, end)
#Now do all the other filter combinations
	toAdd %<>% mutate(V2 = case_when(str_count(V1, ':') == 1 ~ str_replace(V1, ':',' and '),
                               str_count(V1,':') > 1 ~ str_replace_all(V1,':',', '),
                                TRUE ~ V1))
        toAdd %<>% mutate(V2 = ifelse(str_count(V2,',') > 1, sub(".([^,]*)$", " and\\1", V2), V2))
        toAdd %<>% mutate(V3 = paste0(starter,V1, mid, 'Fails ', V2, end))
        #Now sort out the info
        infostart <- '##INFO=<ID='
        infomid <- ',Number=.,Type=Float,Description="'
        infoend <- '">'
    
infos <- setNames(
            as.list(c(
	      'Median depth (taken from the DP FORMAT field) of all samples. Used for filter flag depth.',
              'Median depth (taken from the DP FORMAT field) from samples with non-missing genotypes.',
              'Median genotype quality(taken from the GQ FORMAT field) from samples with non-missing genotypes. Used for filter flag GQ.',
              "Ratio of fully missing genotypes (GT = './.' and FORMAT/DP = 0)",
              'The ratio of complete sites/total number of samples',
              'For each het call, a binomial test is conducted for reads supporting the ref and alt alleles. AB ratio is the hets showing imbalance (p<0.01) divided by the total number of hets.',
              'The number of Mendel Errors at this site, calculated from confirmed trios',
              'HWE mid p-value in inferred unrelated inferred afr superpop.',
              'HWE mid p-value in inferred unrelated inferred amr superpop.',
              'HWE mid p-value in inferred unrelated inferred eas superpop.',
              'HWE mid p-value in inferred unrelated inferred eur superpop.',
              'HWE mid p-value in inferred unrelated inferred sas superpop.'
            )),
            c("medianDepthAll",
              "medianDepthNonMiss",
              "medianGQ",
              "missingness",
              "completeGTRatio",
              "ABratio",
              "MendelSite",
              "phwe_afr",
              "phwe_amr",
              "phwe_eas",
              "phwe_eur",
              "phwe_sas"
             )
    )


#Now build this info full string
infosout <- paste0(infostart, names(infos), infomid, infos, infoend )

#Now print this out to file, and add it to the rest of the header
names(singlecases <- NULL)
d <- c(singlecases, toAdd$V3, infosout) %>% unlist() %>% as.data.frame() 
fwrite(d, 'additional_header.txt', quote = F, col.names = F)
