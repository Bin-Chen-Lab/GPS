unique_drug_profile <- read.csv('unique_drug_profile.csv', stringsAsFactors = F)
unique_drug_value <- read.csv('unique_drug_value.csv', stringsAsFactors = F)

## map gene names
fn = paste0('Predictabilities.csv')
gene_quality_table = read.csv(file=fn, header = T, sep=',', stringsAsFactors = F)   # ge name as L5
go_features_l4 = read.csv(file='../data/output/go_fingerprints_l4.csv', header = T, stringsAsFactors = F)
genes_go_l4 = go_features_l4$gene #not change
genes_low_dose = gene_quality_table$Gene #change to l4 genes
genes_unique_drug_value = colnames(unique_drug_value) #change to l4 genes

"HLA.DMA" %in% genes_low_dose
"HLA.DMA" %in% genes_go_l4
"HLA.DRA" %in% genes_low_dose
"HLA.DRA" %in% genes_go_l4

## important: gene quality, chage HLA-DMA, HLA-DRA to HLA.DMA, HLA.DRA in meta. because both gene quality and go_fingerprints_l4 has . instead of -
gene_meta = read.csv('Landmark_gene_ID_name_3versions.csv', stringsAsFactors = F)
gene_meta[gene_meta=='HLA-DMA']='HLA.DMA'
gene_meta[gene_meta=='HLA-DRA']='HLA.DRA'

all.equal(gene_meta$gene_symbol_2016, gene_meta$gene_symbol_2017)
sum(genes_low_dose %in% gene_meta$gene_symbol_2020) #low dose new gene name
sum(genes_unique_drug_value %in% gene_meta$gene_symbol_2020) #low dose new gene name
sum(genes_go_l4 %in% gene_meta$gene_symbol_2017) # l4 corresponds to 2016/2017

# map low dose to l4 gene names
gene_name_new = rep(NA, nrow(gene_quality_table))
for (i in 1:nrow(gene_quality_table)){
  gn = gene_quality_table$Gene[i]
  idx = which(gene_meta$gene_symbol_2020==gn)
  gene_name_new[i] = gene_meta$gene_symbol_2017[idx]
}
sum(gene_name_new %in% go_features_l4$gene)==length(gene_name_new)
all.equal(gene_name_new, genes_low_dose)
gene_quality_table$Gene=gene_name_new # finish1

gene_name_new = rep(NA, ncol(unique_drug_value))
for (i in 1:ncol(unique_drug_value)){
  gn = genes_unique_drug_value[i]
  idx = which(gene_meta$gene_symbol_2020==gn)
  gene_name_new[i] = gene_meta$gene_symbol_2017[idx]
}
sum(gene_name_new %in% go_features_l4$gene)==length(gene_name_new)
all.equal(gene_name_new, genes_unique_drug_value)
colnames(unique_drug_value)=gene_name_new # finish2

sum(colnames(unique_drug_value) %in% go_features_l4$gene)
sum(gene_quality_table$Gene %in% go_features_l4$gene)
## end of map gene names ####

# select gene based on p and r
p_cut = 0.05
r_cut = 0.3
hist(gene_quality_table$P_MSE)
tmp = which(gene_quality_table$P_MSE < p_cut) # & gene_quality_table$Average_PearsonR > r_cut
if(length(tmp) < 2){
  stop("Too few predictable genes!")
}
selected_genes = gene_quality_table$Gene[tmp]

# reorder genes (drug_value and selected genes affected)
if(length(selected_genes) == 1){
  unique_drug_value = as.data.frame(unique_drug_value[, which(colnames(unique_drug_value) %in% selected_genes)])
  colnames(unique_drug_value) = selected_genes
}else{
unique_drug_value = unique_drug_value[, which(colnames(unique_drug_value) %in% selected_genes)]
tmp = order(colnames(unique_drug_value))
unique_drug_value = unique_drug_value[, tmp]
}

tmp = order(selected_genes)
selected_genes = selected_genes[tmp]
all.equal(colnames(unique_drug_value), selected_genes)

fn = paste0('selected_genes.csv')
write.csv(selected_genes, file=fn, row.names = F, quote = F)

# reorder drugs (drug_profile and drug_value affected)
tmp = order(unique_drug_profile$SMILES)
unique_drug_profile = unique_drug_profile[tmp,]

if(length(selected_genes) == 1){
  unique_drug_value = as.data.frame(unique_drug_value[tmp,])
  colnames(unique_drug_value) = selected_genes
}else{
  unique_drug_value = unique_drug_value[tmp,]
}

# train test split, sample partial high profiles as test set
profile_quality_threshold = 0.7
tmp = which(unique_drug_profile$cor_mean > profile_quality_threshold) 
if(length(tmp) <2){
stop('Too few high-quality profiles!')
}
number_of_HQ_DP = length(tmp)

if (number_of_HQ_DP < 200){
  stop("Not enough samples for training!")
}

n_test = round(0.1 * number_of_HQ_DP)

set.seed(1)
profile_test = sample(tmp, n_test, replace=FALSE)
profile_HQ_small = setdiff(tmp, profile_test)

unique_drug_profile_test = unique_drug_profile[profile_test, ]
if(length(selected_genes) == 1){
  unique_drug_value_test = as.data.frame(unique_drug_value[profile_test, ])
  colnames(unique_drug_value_test) = selected_genes
}else{
  unique_drug_value_test = unique_drug_value[profile_test, ]
}

unique_drug_profile_train = unique_drug_profile[-profile_test, ]

if(length(selected_genes) == 1){
  unique_drug_value_train = as.data.frame(unique_drug_value[-profile_test, ])
  colnames(unique_drug_value_train) = selected_genes
}else{
  unique_drug_value_train = unique_drug_value[-profile_test, ]
}

# flatten drugs-genes matrix to (drug, gene) vector, go_feature is a dictionary, no need to match 
is.unsorted(colnames(unique_drug_value))
unique_drug_value_test_v = c(t(unique_drug_value_test))
unique_drug_profile_test_drug = as.character(rep(unique_drug_profile_test$SMILES, each=ncol(unique_drug_value)))
unique_drug_profile_test_quality = rep(unique_drug_profile_test$cor_mean, each=ncol(unique_drug_value))
unique_drug_profile_test_gene = as.character(rep(selected_genes, nrow(unique_drug_profile_test)))

unique_drug_value_train_v = c(t(unique_drug_value_train))
unique_drug_profile_train_drug = as.character(rep(unique_drug_profile_train$SMILES, each=ncol(unique_drug_value)))
unique_drug_profile_train_quality = rep(unique_drug_profile_train$cor_mean, each=ncol(unique_drug_value))
unique_drug_profile_train_gene = as.character(rep(selected_genes, nrow(unique_drug_profile_train)))

# discretenize y value
# 3-class
labels = c(0, 1, 2)
unique_drug_value_test_v_discrete = labels[as.numeric(cut(unique_drug_value_test_v, c(-99999, -1.5, 1.5, 99999)))]
unique_drug_value_train_v_discrete = labels[as.numeric(cut(unique_drug_value_train_v, c(-99999, -1.5, 1.5, 99999)))]

unique_test = as.data.frame(cbind(unique_drug_profile_test_drug, unique_drug_profile_test_gene, unique_drug_value_test_v, unique_drug_value_test_v_discrete, unique_drug_profile_test_quality))
unique_train = as.data.frame(cbind(unique_drug_profile_train_drug, unique_drug_profile_train_gene, unique_drug_value_train_v, unique_drug_value_train_v_discrete, unique_drug_profile_train_quality))

colnames(unique_test) = c('smiles', 'gene', 'target', 'label', 'quality')
colnames(unique_train) = c('smiles', 'gene', 'target', 'label', 'quality')

unique_test$quality = as.numeric(as.character(unique_test$quality))
unique_train$quality = as.numeric(as.character(unique_train$quality))

write.csv(unique_train, file='data_train.csv', row.names = F, quote = F)
write.csv(unique_test, file='data_test.csv', row.names = F, quote = F)

#--------------------------------------------------------------------
