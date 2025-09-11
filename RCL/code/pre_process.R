load('LINCS_NEW_LOW_DOSE_L4_MCF7.RData')

drug_profile = as.data.frame(zscore[, c('SMILES')])
colnames(drug_profile)='SMILES'
drug_profile$counts = rep(1, nrow(drug_profile))
drug_value = zscore[, 2:ncol(zscore)]  # drug_profile and drug_value correspondence

unique_drug_profile = aggregate(counts ~ SMILES, data = drug_profile, sum)

cor_mean = rep(NA, nrow(unique_drug_profile))
unique_drug_value = matrix(ncol = ncol(drug_value), nrow = nrow(unique_drug_profile))
for (i in 1:nrow(unique_drug_profile)){
  s = unique_drug_profile$SMILES[i]
  
  df = drug_value[which(drug_profile$SMILES == s), ]
  df = as.matrix(df)
  
  n_rep = nrow(df)
  if (n_rep > 1){
    m = colMeans(df)
    v = rep(NA, n_rep)
    for (j in 1:n_rep){
      v[j] = cor(df[j, ], m)
    }
    unique_drug_value[i, ] = t(v)%*%df/sum(v)
    cor_mean[i] = mean(v)
  }else{
    cor_mean[i]=0
    unique_drug_value[i, ] = as.vector(df)
  }
  if (i%%1000==1){print(i)}
}
hist(cor_mean)
unique_drug_profile$cor_mean = cor_mean
colnames(unique_drug_value) = colnames(drug_value)

# remove 1 replicates
tmp = which(unique_drug_profile$counts==1)
if (length(tmp) != 0){
  unique_drug_profile = unique_drug_profile[-tmp, ]
  unique_drug_value = unique_drug_value[-tmp, ]
}

write.csv(unique_drug_profile, 'unique_drug_profile.csv', row.names = F)
write.csv(unique_drug_value, 'unique_drug_value.csv', row.names = F)


