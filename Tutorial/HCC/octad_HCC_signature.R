library('octad')

#1. load phenotype data
phenoDF <- octad.db::get_ExperimentHub_data("EH7274")
head(phenoDF)[1:10] #visualize phenotype table
HCC_primary=subset(phenoDF,cancer=='liver hepatocellular carcinoma'& sample.type == 'primary') #select data
case_id=HCC_primary$sample.id #select cases

#2. Compute reference tissue.
#computing reference tissue
HCC_adjacent=subset(phenoDF,cancer=='liver hepatocellular carcinoma'&sample.type == 'adjacent'&data.source == 'TCGA') #select data
control_id=HCC_adjacent$sample.id #select cases

#3. Compute DE
#compute differential expression, filter result file
res=diffExp(case_id,control_id,source='octad.whole',output=T,n_topGenes=10000,file='octad.counts.and.tpm.h5')
res=subset(res,abs(log2FoldChange)>1&padj<0.05)

#4. Prepare into GPS format
res$Value <- res$log2FoldChange
res$GeneSymbol <- res$Symbol
res <- res[, c("GeneSymbol", "Value")]
write.csv(res, "DZSIG__HCC.csv", row.names = FALSE)
