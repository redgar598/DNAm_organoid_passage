load(file=here("data/validation/DNAm","validation_betas_normalized.RData"))
load(here("data","beta_organoids.RData"))



#' ### representative CpGs
epic.organoid_minimal<-epic.organoid[,c(2, 14, 17)]
colnames(epic.organoid_minimal)[1]<-"Assay.Name"
epic.organoid_minimal$cohort<-"Cohort 1 Organoids"

validation_epic.organoid$Sample_ID<-paste(validation_epic.organoid$individual, validation_epic.organoid$Segment, validation_epic.organoid$condition)
validation_epic.organoid_minimal<-validation_epic.organoid[,c(10, 19, 12)]
colnames(validation_epic.organoid_minimal)[1]<-"Assay.Name"
validation_epic.organoid_minimal$cohort<-"Validation Organoids"
colnames(validation_epic.organoid_minimal)[2:3]<-c("sample_ID","passage.or.rescope.no_numeric")

sample_info_both<-rbind(validation_epic.organoid_minimal,epic.organoid_minimal)
sample_info_both$passage.or.rescope.no_numeric<-as.numeric(as.character(sample_info_both$passage.or.rescope.no_numeric))

validation_organoid_beta<-validation_organoid_beta[which(rownames(validation_organoid_beta)%in%rownames(organoid_beta)),]
organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%rownames(validation_organoid_beta)),]

validation_organoid_beta<-validation_organoid_beta[match(rownames(organoid_beta), rownames(validation_organoid_beta)),]
identical(rownames(validation_organoid_beta), rownames(organoid_beta))

betas_both<-cbind(validation_organoid_beta,organoid_beta)

identical(colnames(betas_both),sample_info_both$Assay.Name )


#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(betas_both))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

## PC vs PC plot
Loadings$Assay.Name<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sample_info_both, by="Assay.Name")

ggplot(Loadings_meta, aes(PC1, PC2, fill=cohort))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (9%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("red", "blue"))


ggplot(Loadings_meta, aes(PC2, PC3, fill=cohort))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (9%)")+ylab("PC3 (5%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("red", "blue"))


ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (9%)")+ylab("PC3 (5%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col)





#' ## Overall Variance Across most Variable CpGs with Passage
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
Mval<-function(beta) log2(beta/(1-beta))

organoid_mval = apply(betas_both, 1, Mval)
organoid_mval = as.data.frame(organoid_mval)
organoid_mval = t(organoid_mval)

ref_range_dnam<-sapply(1:nrow(organoid_mval), function(x) Variation(organoid_mval[x,]))
validation_organoid_beta_VeryVariable<-betas_both[which(ref_range_dnam>=2.75),]

print(paste("There are ",nrow(validation_organoid_beta_VeryVariable), " variable CpGs (10th-90th quantile range in M value >2.75)",sep=""))

#' dim(validation_organoid_beta_VeryVariable<-validation_organoid_beta[rev(order(ref_range_dnam)),])
#' #' Include the same number of vairable CpGs as organoids varible. So take the 71384 most variable
#' dim(validation_organoid_beta_VeryVariable<-validation_organoid_beta_VeryVariable[1:71384 ,])


# beta plot variable CpGs
Beta_melted<- melt(validation_organoid_beta_VeryVariable)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sample_info_both, by.x="ID", by.y="Assay.Name")

ggplot(Beta_Plot, aes(Beta,  color=as.factor(passage.or.rescope.no_numeric)))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")

ggplot(Beta_Plot, aes(Beta,  color=as.factor(passage.or.rescope.no_numeric)))+
  geom_density(size=0.75)+theme_bw()+xlab("DNAm Beta Value")+ylab("Density")+
  scale_color_manual(values=pass_col, name="Passage\nNumber")+facet_wrap(~cohort)
