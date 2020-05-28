suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lmtest)
library(gridExtra)
library(gtools)
library(scales)
library(here)
library(binom)


options(stringsAsFactors = FALSE)
options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)
x<-as.numeric(args[1])

#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))


Mval<-function(beta) log2(beta/(1-beta))
organoid_Mval = apply(organoid_beta, 1, Mval) # need mvalues for combat
organoid_Mval = as.data.frame(organoid_Mval)
organoid_Mval = t(organoid_Mval)

# variable mvalues only
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
ref_range_dnam<-sapply(1:nrow(organoid_Mval), function(x) Variation(organoid_Mval[x,]))
organoid_beta_VeryVariable<-organoid_beta[which(ref_range_dnam>=2.75),]#  51545

print(paste("There are ",nrow(organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))


converted<-as.numeric(round(organoid_beta_VeryVariable[,x]*1000,0))
counts<-rep(1000, length(converted))
res = em(converted, counts, .419, .31, .27, 0.01, .1, .1, .90, .03, .5, .05)
passage<-paste("passage: ",epic.organoid$passage.or.rescope.no_numeric[x],"\nIndividual: ", epic.organoid$case.no[x],"\nPrior I: ",round(res$prior_I,2), sep="")
draw_fit_params_gg(converted, counts, res,passage)

ggsave(paste(here("figs"),"/tmp_",x,"_EM_iter.pdf", sep=""), width=6, height=5)
