#!/usr/bin/env Rscript

### SEX CHROM ANNOTATION ###

library(data.table)
library(magrittr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly = T)

argList <- c('pos',
             'start',
             'missing1',
             'missing2',
             'coverageAll',
             'coverageNonMiss',
             'quality',
             'allHets',
             'passHets',
             'mendel',
             'outdir',
             'nsamp')

# The arguments provided to this script in the order above should not be actual files,
# but rather basenames of files that have _XX _XY endings where that is applicable.
# The endings are automatically then added by the script itself and files are properly handled.

allArgs <- list()
allArgs <- lapply(1:length(argList),
                  function(i){
                    allArgs <- as.character(args[i])
                  })
names(allArgs) <- argList
#Also adding the hardy location here
allArgs$hardy <- '.'


# Funcs -------------------------------------------------------------------
startProc <- function(x){
  dat <- fread(x$start) %>%
    as_tibble() %>%
    mutate(ID = paste0(V1, ":", V2, "-", V3, "/", V4, "-", V5)) %>%
    select(ID, everything()) %>%
    setNames(c("ID","#CHROM","POS","REF","ALT", "OC_OM")) %>%
    select(-OC_OM)
  return(dat)
}

missProc <- function(dat, x){
  miss1 <- fread(x$missing1) %>%
    as_tibble() %>%
    setNames(c("ID", "miss1"))

  miss2 <- fread(x$missing2) %>% as_tibble() %>% setNames(c("ID", "miss2"))

  miss <- full_join(miss1, miss2, by="ID") %>%
    mutate(miss1 = replace_na(miss1, 0),
           miss2 = replace_na(miss2, 0)) %>%
    mutate(missingness = miss1/(miss2 + miss1)) %>%
    select(-miss1)

  dat %<>% left_join(miss, by="ID") %>%
    mutate(missingness = replace_na(missingness, 0)) #sites as NAs here should
  #be set to fully missing
}

covProc <- function(dat, x){
  covAll <- fread(x$coverageAll) %>%
    as_tibble() %>%
    setNames(c("ID", "medianCovAll"))
  covNonMiss <- fread(x$coverageNonMiss) %>%
    as_tibble() %>%
    setNames(c("ID", "medianCovNonMiss"))

  dat %<>%
    left_join(covAll, by = "ID") %>%
    left_join(covNonMiss, by = "ID") %>%
    mutate(medianCovNonMiss = replace_na(medianCovNonMiss, 0))#sites as NAs here should
  #be set to 0 coverage
}

gqProc <- function(dat, x){
  GQ <- fread(x$quality) %>%
    as_tibble() %>%
    setNames(c("ID", "medianGQ"))

  dat %<>% left_join(GQ, by="ID") %>%
    mutate(medianGQ = replace_na(medianGQ, 0))#sites as NAs here should
  #be set to 0 coverage
}

abProc <- function(dat, x){
  hetAll <- fread(x$allHets) %>%
    as_tibble() %>%
    rename('ID'='V1', 'AB_hetAll'='V2')

  hetPass <- fread(x$passHets) %>%
    as_tibble() %>%
    as_tibble() %>%
    rename('ID'='V1', 'AB_hetPass'='V2')

  #We expect fewer passes than all hets, so NAs will be introduced,
  #set these to 0
  ab <- hetAll %>%
    left_join(hetPass, by='ID') %>%
    mutate(AB_hetPass = replace_na(AB_hetPass, 0)) %>%
    mutate(AB_Ratio = AB_hetPass/AB_hetAll) %>%
    select(ID, AB_Ratio)

  dat %<>%
    left_join(ab, by = "ID") %>%
    mutate(AB_Ratio = replace_na(AB_Ratio, 99)) #scores of 99 will be changed downstream
}

completeProc <- function(dat, x){
  nsamples <- fread(x$nsamp) %>%
    as_tibble()
  dat %<>% mutate( miss2 = replace_na(miss2, 0),
                   completeSites = 1-((nsamples$V1-miss2)/nsamples$V1))
}

hweProc <- function(dat){
  hardyfiles <- paste0(c('AFR','EUR','EAS','SAS'))
  hwe <- lapply(hardyfiles,
                function(pop) {
                  f <- fread(file.path(allArgs$hardy,
                                       paste0(allArgs$pos,
                                              '_',
                                              pop,
                                              '.hwe'))) %>%
                    as_tibble()
                  names(f) <- c(names(f)[-ncol(f)], paste0('phwe_',tolower(pop)))
                  f <- f[,ncol(f)]
                })
  #Remember format of plink files switches ref and alt, aside from that the
  #vars should be in the same order
  out <- bind_cols(hwe)
  dat %<>% bind_cols(out)
  return(dat)
}

# MendelErrors ------------------------------------------------------------
#Only relevant to the female samples
mendelProc <- function(dat, x){
  mend <- paste0(gsub('_XX','',x$mendel))
  me <- fread(mend) %>%
    as_tibble()
  if(nrow(me) != nrow(dat)){
    stop('Mendel error site file length does not match backbone file')
  } else {
    dat$MendelSite <- me$N
  }
  return(dat)
}

isPAR <- function(){
  #Change FILTER if we are in a PAR region
  #PAR regions - https://m.ensembl.org/info/genome/genebuild/human_PARS.html
  #PAR1 - chromosome:GRCh38:Y:10001 - 2781479 is shared with X: 10001 - 2781479 (PAR1)
  #chromosome:GRCh38:Y:56887903 - 57217415 is shared with X: 155701383 - 156030895 (PAR2)
  #First identify the region, then split the chunk into PAR-nonPAR.
  par1 <- c(10001, 2781479)
  par2 <- c(56887903, 57217415)
  #Make a sequence of positions possible in chunk, and then pull pos overlapping
  #par regions
  positions <- strsplit(allArgs$pos, split = '_', fixed = T)[[1]][2:3] %>%
    as.numeric()

  inPAR <- seq(from = positions[1],
               to= positions[2],
               by = 1)
  inPAR <- inPAR[(inPAR %between% c(par1[1], par1[2])) |
                   inPAR %between% c(par2[1], par2[2])]

  return(inPAR)
}

XX_PAR_FILTER <- function(dat, XX){
  #SET XX TO TRUE/FALSE accordingly
  #These filters are used for all females, and for male in PAR regions
  if(XX){
    dat %<>%
      mutate(FILTER =
               ifelse((missingness <= 0.05 &
                         medianCovAll >= 10 &
                         medianGQ >= 15 &
                         completeSites >= 0.5 &
                         AB_Ratio >= 0.25 &
                         phwe_eur >= 10e-6),
                      "PASS", 'NA'),
             FILTER = ifelse(missingness > 0.05, paste0(FILTER, ':missingness'), FILTER),
             FILTER = ifelse(medianCovAll < 10, paste0(FILTER, ':depth'), FILTER),
             FILTER = ifelse(AB_Ratio < 0.25, paste0(FILTER, ':ABratio'), FILTER),
             FILTER = ifelse(completeSites < 0.5, paste0(FILTER, ':completeGTRatio'), FILTER),
             FILTER = ifelse(medianGQ < 15, paste0(FILTER, ':GQ'), FILTER),
             FILTER = ifelse(phwe_eur < 10e-6, paste0(FILTER, ':phwe_eur'), FILTER)) %>%
      mutate(FILTER =
               ifelse(grepl('^NA:', FILTER), str_sub(FILTER, 4), FILTER))
  } else {
    cat('Annotating PAR region...\n')
    dat %<>%
      mutate(FILTER = #NOTE WE DON'T INCLUDE AB RATIO FOR MALES
               ifelse((missingness <= 0.05 &
                         medianCovAll >= 10 &
                         medianGQ >= 15 &
                         completeSites >= 0.5),
                      "PASS_m", 'NA'),
             FILTER = ifelse(missingness > 0.05, paste0(FILTER, ':missingness_m'), FILTER),
             FILTER = ifelse(medianCovAll < 10, paste0(FILTER, ':depth_m'), FILTER),
             FILTER = ifelse(completeSites < 0.5, paste0(FILTER, ':completeGTRatio_m'), FILTER),
             FILTER = ifelse(medianGQ < 15, paste0(FILTER, ':GQ_m'), FILTER)) %>%
      mutate(FILTER =
               ifelse(grepl('^NA:', FILTER), str_sub(FILTER, 4), FILTER))
  }
}

annotateXX <- function(){
  xxArgs <- lapply(allArgs, function(x) paste0(x, '_XX'))
  dat <- startProc(xxArgs)
  cat('Start file done...\n')
  dat %<>% missProc(xxArgs)
  cat('Missingness done...\n')
  dat %<>% covProc(xxArgs)
  cat('Coverage done...\n')
  dat %<>% gqProc(xxArgs)
  cat('GQ done...\n')
  dat %<>% abProc(xxArgs)
  cat('AB ratio done...\n')
  dat %<>% completeProc(xxArgs)
  dat %<>% hweProc()
  cat('HWE...\n')
  dat %<>% mendelProc(xxArgs)
  dat %<>% XX_PAR_FILTER(XX=T)
}

annotateXY <- function(){
  xyArgs <- lapply(allArgs, function(x) paste0(x, '_XY'))
  dat <- startProc(xyArgs)
  cat('Start file done...\n')
  dat %<>% missProc(xyArgs)
  cat('Missingness done...\n')
  dat %<>% covProc(xyArgs)
  cat('Coverage done...\n')
  dat %<>% gqProc(xyArgs)
  cat('GQ done...\n')
  dat %<>% completeProc(xyArgs)
  #Pull regions that are in PAR
  par <- isPAR()
  if(length(par) > 0){
    dat %<>% mutate(index = 1:nrow(dat))
    datPAR <- dat %>% filter(POS %in% par)
    datPAR %<>% XX_PAR_FILTER(XX = F)
    dat %<>% filter(!POS %in% par)
  }
  #Non-par bits
  dat %<>%
    mutate(FILTER =
             ifelse((missingness <= 0.05 &
                       medianCovAll >= 5 &
                       medianGQ >= 15 &
                       completeSites >= 0.5),
                    "PASS_m", 'NA'),
           FILTER = ifelse(missingness > 0.05, paste0(FILTER, ':missingness_m'), FILTER),
           FILTER = ifelse(medianCovAll < 5, paste0(FILTER, ':depth_m'), FILTER),
           FILTER = ifelse(completeSites < 0.5, paste0(FILTER, ':completeGTRatio_m'), FILTER),
           FILTER = ifelse(medianGQ < 15, paste0(FILTER, ':GQ_m'), FILTER)) %>%
    mutate(FILTER =
             ifelse(grepl('^NA:', FILTER), str_sub(FILTER, 4), FILTER))

  if(length(par) > 0){
    dat <- bind_rows(datPAR, dat) %>%
      arrange(index) %>%
      select(-index)
  }
  names(dat) <- c(names(dat[1:5]), paste0(names(dat)[6:ncol(dat)], '_m'))
  return(dat)
}

combine_sets <- function(XX, XY){
  XX %>% full_join(XY, by=(c('ID'='ID',
                                   '#CHROM'='#CHROM',
                                   'POS'='POS',
                                   'REF'='REF',
                                   'ALT'='ALT')))
}


final_clean <- function(dat){
  dat %<>%
    mutate_at(.vars = vars(-ID,
                           -`#CHROM`,
                           -POS,
                           -REF,
                           -ALT,
                           -FILTER,
                           -FILTER_m),
              .funs = list( ~ gsub("(\\.?<![0-9])0+", "",
                                   round(x = ., digits = 3),
                                   perl = TRUE)))  %>%
    mutate_at(.vars = vars(),
              .funs = list(~ replace_na(., '.')))
  dat %<>% mutate(AB_Ratio = ifelse(AB_Ratio == 99, '.', AB_Ratio))
  dat %<>% select(-ID, everything())
  dat %<>% select(-'miss2', -'miss2_m')
  dat %<>% rename('medianDepthAll'='medianCovAll',
                  'medianDepthNonMiss'='medianCovNonMiss',
                  'completeGTRatio'='completeSites',
                  'medianDepthAll_m'='medianCovAll_m',
                  'medianDepthNonMiss_m'='medianCovNonMiss_m',
                  'completeGTRatio_m'='completeSites_m',
                  'ABratio'='AB_Ratio')
  return(dat)
}



# Run ---------------------------------------------------------------------
datXX<- annotateXX()

datXY <- annotateXY()

#Combine data
datOut <- combine_sets(datXX, datXY) %>%
  final_clean()
datOut %>%
  write.table(paste0(allArgs$outdir,'/BCFtools_site_metrics_',allArgs$pos, '.txt'),
              quote=F,
              row.names = F,
              sep = '\t')


