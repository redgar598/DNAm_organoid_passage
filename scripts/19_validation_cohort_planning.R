library(here)
options(stringsAsFactors = FALSE)

samples<-read.csv(here("data/validation_dataset/Final samples for submission 2021.csv"))

manifest<-read.csv(here("data/validation_dataset/Methylation_infinium_manifest_template.csv"), skip=12, header=T)
manifest$sentrix<-sapply(1:nrow(manifest), function(x) strsplit(manifest$Well[x], ":")[[1]][2] )
manifest<-manifest[-97,]

manifest_1<-manifest
manifest_2<-manifest
manifest_2$Plate.ID<-"Plate 2"
manifest_2$sentrix<-paste(manifest_2$sentrix, "_2",sep="")
manifest_3<-manifest[1:16,]
manifest_3$Plate.ID<-"Plate 3"
manifest_3$sentrix<-paste(manifest_3$sentrix, "_3",sep="")

manifest<-rbind(manifest_1, manifest_2, manifest_3)

table(manifest$Plate.ID)
table(manifest$sentrix)

# will need 3 plates, 2 full at 96 samples, on with only 16 samples (2 arrays)

## will place the replicates and randomized the rest
# 2 reps on plate 3, 3 reps on plate 1 and 2
set.seed(5)
manifest_reps<-rbind(manifest_1[sample(1:96,3),],manifest_2[sample(1:96,3),], manifest_3[which(manifest_3$sentrix=="1_3"),][sample(1:8,1),],manifest_3[which(manifest_3$sentrix=="2_3"),][sample(1:8,1),])
manifest_reps
manifest_reps$Name<-sample(samples$Sample.Name[which(samples$Biobank.Rachel.Replicates=="Replicates")])
manifest_reps


set.seed(5)
manifest_samples<-rbind(manifest_1[-sample(1:96,3),],manifest_2[-sample(1:96,3),], manifest_3[which(manifest_3$sentrix=="1_3"),][-sample(1:8,1),],manifest_3[which(manifest_3$sentrix=="2_3"),][-sample(1:8,1),])

set.seed(20)
manifest_samples$Name<-sample(samples$Sample.Name[which(samples$Biobank.Rachel.Replicates!="Replicates")])

manifest_samples<-rbind(manifest_samples, manifest_reps)


#check balanace of all things
manifest_samples_merge<-merge(manifest_samples, samples, by.x="Name",by.y="Sample.Name")


## plate balanced?
balance_plate<-lapply(c(11:13,15:19,30), function(x) table(manifest_samples_merge[,x],manifest_samples_merge$Plate.ID))
names(balance_plate)<-colnames(manifest_samples_merge)[c(11:13,15:19,30)]
balance_plate

## sentrix balanced?
balance_sentrix<-lapply(c(11:13,15:19,30), function(x) table(manifest_samples_merge[,x],manifest_samples_merge$sentrix))
names(balance_sentrix)<-colnames(manifest_samples_merge)[c(11:13,15:19,30)]
balance_sentrix

manifest_samples_merge$sentrix_numeric<-as.numeric(sapply(1:nrow(manifest_samples_merge), function(x) strsplit(manifest_samples_merge$Well[x], ":")[[1]][2] ))

manifest_samples_merge<-manifest_samples_merge[order(manifest_samples_merge$Plate.ID, manifest_samples_merge$sentrix_numeric,manifest_samples_merge$Well),]


#  They need gender column? in manifest? and concentration and volume and total ng?
## Komal will handel but copying and pasting sample names into fncy format sample sheets

write.csv(file="data/validation_dataset/sample_order_randomized.csv",manifest_samples_merge)

dim(samples)
length(unique(intersect(samples$Sample.Name, manifest_samples_merge$Name)))
