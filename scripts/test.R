#'---
#'title: Exploration of CpGs effected by passage
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' #### Load Libraries
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(scales))
suppressMessages(library(here))
suppressMessages(library(binom))


options(stringsAsFactors = FALSE)
options(scipen = 999)


#' #### Load Functions
source(here("scripts","00_pretty_plots.R"))
source(here("scripts","00_EM_array_uniform_background_maximise_betabinom.R"))


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))


#' ##### Beta value distribution change with passage

# beta plot all CpGs
Beta_melted<- melt(organoid_beta)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
Beta_Plot$passage.or.rescope.no_numeric.factor <- factor(Beta_Plot$passage.or.rescope.no_numeric, levels = c(16,14,11,10,8,7,6,4,3,2,1))

ggplot(Beta_Plot, aes(Beta,color=passage.or.rescope.no_numeric.factor))+
  geom_density(size=0.5)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col,name="Passage\nNumber")+th+theme(legend.text = element_text(size=7),
                                                                      legend.title = element_text(size=10),
                                                                      legend.key.size = unit(0.7,"line"))


ggsave(here("figs","Passage_all_CpGs.pdf"),width = 3.75, height = 2.5)
ggsave(here("figs/jpeg","Passage_all_CpGs.jpeg"), width = 3.75, height = 2.5)

