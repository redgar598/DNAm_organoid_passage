

#rmarkdown::render("scripts/01_organoids_preprocessing.R", output_dir = "output")
#rmarkdown::render("scripts/04_passage_CpGs.R", output_dir = "output")

#rmarkdown::render("scripts/09_enrichment_heteroskedasticity_CpG.R", output_dir = "output")
#rmarkdown::render("scripts/10_ensembl_genes_EPIC_annotation.R", output_dir = "output")
#rmarkdown::render("scripts/11_EPIC_gene_annotation.R", output_dir = "output")

rmarkdown::render("scripts/13_prep_For_GO_ORA.R", output_dir = "output")
rmarkdown::render("scripts/14_post_GO_ORA.R", output_dir = "output")

rmarkdown::render("scripts/15_cancer_beta_simplified.R", output_dir = "output")
#rmarkdown::render("scripts/16_MTAB4957_passage.R", output_dir = "output")
rmarkdown::render("scripts/17_GSE141256_organoids.R", output_dir = "output")
#rmarkdown::render("scripts/18_GSE142213_organoids.R", output_dir = "output")


