suppressPackageStartupMessages({
  library(nychanesmicrobiome)
  library(dplyr)
  })

NYC_HANES <- loadQiimeData(sas7bdat::read.sas7bdat("/home/km/work/GSEA_paper/DSAScaleError/nyc-hanes-datasets-and-resources-analytic-data-sets-sas-file.sas7bdat")) %>% annotateFactors(.)
# Extract alternative smokers who also smoke cigarettes and remove them
alt_smokers_cigarettes <- data.frame(sample_data(NYC_HANES)) %>% dplyr::filter(smokingstatus == 'Alternative smoker') %>% dplyr::filter(CIGARETTES == 'Yes') %>% dplyr::select(Burklab_ID) %>% t
NYC_HANES <- prune_samples(!(sample_names(NYC_HANES) %in% alt_smokers_cigarettes), NYC_HANES)
sample_data(NYC_HANES)$smokingstatus <- relevel(sample_data(NYC_HANES)$smokingstatus,"Never smoker")

# Subset alternative smokers
NYC_HANES.alt <- subset_samples(NYC_HANES, smokingstatus %in% c("Never smoker", "Alternative smoker"))
levels(sample_data(NYC_HANES.alt)$smokingstatus) <- c('Never','Alternative')

biosis.tsv <- system.file("extdata","biosis.tsv", package="nychanesmicrobiome", mustWork = TRUE)
biosis <- read.csv(biosis.tsv,header = TRUE, sep = '\t')

NYC_HANES.genus <- tax_glom(NYC_HANES, taxrank = "Genus")
NYC_HANES.genus.relab <-
  transform_sample_counts(NYC_HANES.genus, function(x)
    x / sum(x))
NYC_HANES.relab <-
  transform_sample_counts(NYC_HANES, function(x)
    x / sum(x))

NYC_HANES.genus.relab <- merge_taxa(NYC_HANES.genus,
                                    taxa_names(filter_taxa(NYC_HANES.genus.relab, function(x)
                                      mean(x) < 2e-4, TRUE)))
tax_table(NYC_HANES.genus.relab)[is.na(tax_table(NYC_HANES.genus.relab)[, "Genus"]), "Genus"] <-
  "Other"
NYC_HANES.genus.relab <-
  tax_glom(transform_sample_counts(NYC_HANES.genus.relab, function(x)
    x / sum(x)),
    taxrank = "Genus")

NYC_HANES.phylum <- tax_glom(NYC_HANES, taxrank = "Phylum")
NYC_HANES.phylum.relab <-
  transform_sample_counts(NYC_HANES.phylum, function(x)
    x / sum(x))

full_eset <- ExpressionSet(assayData=as.matrix(otu_table(NYC_HANES)), phenoData = new("AnnotatedDataFrame", sample_data(NYC_HANES)))
full_eset$CATCOTININE <- ifelse(full_eset$COTININE < 3,'low','high')
full_annotated_ds <- tax_table(NYC_HANES) %>% as.data.frame %>% bind_cols(data.frame(OTU=rownames(tax_table(NYC_HANES)))) %>% filter(Domain != "Unassigned") %>% left_join(biosis, by = c("Genus"="X1"))

alt_full_eset <- full_eset[,(full_eset$smokingstatus %in% c("Never smoker","Alternative smoker"))]
alt_full_eset$smokingstatus <- droplevels(alt_full_eset$smokingstatus)
levels(alt_full_eset$smokingstatus) <- c('Never','Alternative')
alt_full_eset$GROUP <- ifelse(alt_full_eset$smokingstatus == "Alternative",1,0)

EnrichmentBrowser::configEBrowser("OUTDIR.DEFAULT","./results")

sns.eset <- full_eset[,full_eset$smokingstatus %in%  c("Cigarette","Never smoker")]
NYC_HANES.sns <- NYC_HANES %>% subset_samples(., smokingstatus %in% c("Cigarette","Never smoker"))
sns.eset$GROUP <- ifelse(sns.eset$smokingstatus == "Cigarette", 1, 0)

design = model.matrix(~smokingstatus, data=data.frame(sample_data(NYC_HANES.sns)))

sns_ds <- full_annotated_ds %>% filter(OTU %in% taxa_names(NYC_HANES.sns))
sns_sAero <- (sns_ds %>% filter(X2=="Aerobic"))$OTU %>% as.character
sns_sAnae <- (sns_ds %>% filter(X2=="Anaerobic"))$OTU %>% as.character
sns_sFana <- (sns_ds %>% filter(X2=="F Anaerobic" | X2=="Aero / Facultative Anaerobic"))$OTU %>% as.character
sns_my.gs <- list(aero=sns_sAero, anae=sns_sAnae, fana=sns_sFana)

sns_ds <- full_annotated_ds %>% filter(OTU %in% taxa_names(NYC_HANES.sns))
sns_sAero <- (sns_ds %>% filter(X2=="Aerobic"))$OTU %>% as.character
sns_sAnae <- (sns_ds %>% filter(X2=="Anaerobic"))$OTU %>% as.character
sns_sFana <- (sns_ds %>% filter(X2=="F Anaerobic" | X2=="Aero / Facultative Anaerobic"))$OTU %>% as.character
sns_my.gs <- list(aero=sns_sAero, anae=sns_sAnae, fana=sns_sFana)
design = model.matrix(~smokingstatus, data=data.frame(sample_data(NYC_HANES.sns)))
x <- otu_table(NYC_HANES)[,colnames(otu_table(NYC_HANES))%in%row.names(design)]
saveRDS(x, "~/data/output/smoking_vs_non_smoking.RDS")
saveRDS(design[,2], "~/data/output/smoking_vs_non_smoking_categories.RDS")
saveRDS(sns_my.gs, "~/data/output/smoking_vs_non_smoking_microbe_sets.RDS")
