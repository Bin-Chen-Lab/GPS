library('octad')

data = read.csv('res_MUC5B+.csv', row.names = 1)
  
sRGES = runsRGES(data, max_gene_size=100, permutations=10000)
  
sRGES = sRGES[sRGES$sRGES < -0.2,]
sRGES = sRGES[!grepl("BRD-", sRGES$pert_iname),]

write.csv('octad_MUC5B+.csv', srges_nom)
