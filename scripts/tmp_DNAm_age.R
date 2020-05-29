dnamage<-read.csv("../ibd/output/clock/ibd_beta_organoid_funnorm.output.csv")
dnamage$SampleID<-gsub("X","",dnamage$SampleID)


#' #### Load Normalized Data
load(here("data","beta_organoids.RData"))


dnamage_80<-dnamage[which(dnamage$SampleID%in%epic.organoid$array.id),]

sampleinfo<-merge(dnamage_80, epic.organoid, by.x="SampleID", by.y="array.id")



#overall correlation plot
ggplot(sampleinfo, aes(Age, DNAmAge))+geom_abline(color="lightgrey")+geom_point(size=2, color="black", fill="grey", shape=21)+
  stat_smooth(method="lm", se=F, color="black")+th+theme_bw()+xlab("Chronological Age")+
  annotate("text", x=20, y=5, label = paste("r = ", round(cor(sampleinfo$Age, sampleinfo$DNAmAge), 2)), size = 5)+
  xlim(0,25)+ylim(0,25)

ggsave(here("figs/jpeg","tmp_organoid_Age_cor.jpeg"),width = 5, height = 5)


ggplot(sampleinfo, aes(sample.site, AgeAccelerationResidual, fill=sample.site))+geom_boxplot(color="black", outlier.shape=NA)+
  geom_point(shape=21,position=position_jitter(w=0.15), size=2)+theme_bw()+scale_fill_manual(values=c("#a6d96a","cornflowerblue"), name="Sample\nSite")+th+xlab("")+ylab("Age Acceleration Residual ")
ggsave(here("figs/jpeg","tmp_organoid_Age_acc.jpeg"),width = 5, height = 5)


ggplot(sampleinfo, aes(sample.site, AgeAccelerationResidual, fill=sample.site))+
  geom_boxplot(color="black", outlier.shape=NA)+
  geom_line(aes(group=case.no), color="grey70")+
  geom_point(shape=21, size=2)+
  theme_bw()+
  scale_fill_manual(values=c("#a6d96a","cornflowerblue"), name="Sample\nSite")+th+xlab("")+ylab("Age Acceleration Residual ")
ggsave(here("figs/jpeg","tmp_organoid_Age_acc_paired.jpeg"),width = 5, height = 5)


summary(aov(sampleinfo$AgeAccelerationResidual~sampleinfo$sample.site))

## paired T test

sampleinfo$case.pass<-paste(sampleinfo$case.no, sampleinfo$passage.or.rescope.no_numeric, sep="_")

paired_ageacc<-do.call(rbind,lapply(as.character(unique(sampleinfo$case.pass)), function(x){
  matched<-sampleinfo[which(sampleinfo$case.pass==x),]
  if(length(unique(matched$sample.site))<2){}else{
    data.frame(case.no=unique(matched$case.no),
               TI=matched[which(matched$sample.site=="TI"),"AgeAccelerationResidual"],
               SC=matched[which(matched$sample.site=="SC"),"AgeAccelerationResidual"])}}))

t.test(paired_ageacc$TI, paired_ageacc$SC, paired = TRUE)


ggplot(sampleinfo, aes(passage.or.rescope.no_numeric,AgeAccelerationResidual, 
                          fill=as.factor(passage.or.rescope.no_numeric)))+
  geom_line(aes(passage.or.rescope.no_numeric,AgeAccelerationResidual, group=sample_ID), color="lightgrey")+
  geom_point(shape=21, size=3)+scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+th+theme_bw()

ggsave(here("figs/jpeg","tmp_organoid_Age_acc_passage.jpeg"),width = 7, height = 5)
