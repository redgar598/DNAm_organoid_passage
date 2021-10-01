#' ## DNAm passage genes
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))

hypo_diffexp_original<-sleuth_significant_UD[which(sleuth_significant_UD$ext_gene%in%diff_genes_db_hypovalidation_original$V1),]
hyper_diffexp_original<-sleuth_significant_UD[which(sleuth_significant_UD$ext_gene%in%diff_genes_db_hypervalidation_original$V1),]

head(hypo_diffexp_original[order(hypo_diffexp_original$pval),])
head(hyper_diffexp_original[order(hyper_diffexp_original$pval),])

length(unique(c(hyper_diffexp_original$ext_gene, hypo_diffexp_original$ext_gene)))


diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))

sleuth_sig_DNAm<-low_not_high_diff[which(low_not_high_diff$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]



low_not_high_treatment<-sleuth_significant_low[which(!(sleuth_significant_low$ext_gene%in%sleuth_significant_high$ext_gene)),]

#' ## DNAm passage genes original 80 and validation
diff_genes_db_hypovalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypomethylation_UTUD.txt"))
diff_genes_db_hypervalidation_original<-read.table(file=here("data/validation/DNAm/","validation_original_genes_hypermethylation_UTUD.txt"))
sleuth_sig_DNAm<-low_not_high_treatment[which(low_not_high_treatment$ext_gene%in%c(diff_genes_db_hypovalidation_original$V1, diff_genes_db_hypervalidation_original$V1)),]
