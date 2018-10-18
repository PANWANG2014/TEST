######################################################################################################

   ### Part I. FUNCTIONS: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#the following three functions are used to match the names of samples, mirnas and mrnas
#By using these three functions, we can get the matched mrna expression data and mirna expression data
name_func = function(name_base, mirna_base, mrna_base){
  x = match(name_base[, "miRNA"], mirna_base[1, ])
  y = match(name_base[, "mRNA"], mrna_base[1, ])
  return(list(x,y))
}
mirna_matrix = function(name_base, mirna_base, mirna_name, mirna){  
  mirna_name1 = rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i] = unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use = mirna_base[mirna_name1%in%mirna, mirna_name]
  mirna_name1 = mirna_name1[mirna_name1%in%mirna]
  return(list(mirna_use, mirna_name1))
}
mrna_matrix = function(name_base, mrna_base, mrna_name, mrna){
  mrna_exist = rep(0,nrow(mrna_base))
  mrna_name_all = matrix(data = 0, ncol = 2, nrow = nrow(mrna_base), byrow = TRUE)
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 1] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    mrna_name_all[i, 2] = unlist(strsplit(mrna_base[i, 1], "\\|"))[2]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name1 = mrna_name1[mrna_exist]
  mrna_name_sp = mrna_name_all[mrna_exist, ]
  mrna_fullname = mrna_base[mrna_exist, 1]
  return(list(mrna_use, mrna_name1, mrna_name_sp, mrna_fullname))
}
#mirna_mrna_data_unselected will select the RNAs that are expressed in at least one sample
#mirna_mrna_data is used for targetscan, which will select the RNAs that are expressed in a certain percentage of samples
mirna_mrna_data_unselected = function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_sp, mrna_fullname){
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0
  mirna_sgn = seq(1, nrow(mirna_use), 1)
  mrna_sgn = seq(1, nrow(mrna_use), 1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) == ncol(mirna_use)){
      mirna_sgn[i] = 0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) == ncol(mrna_use)){
      mrna_sgn[i] = 0
    }
  }  
  mirna_use1 = mirna_use[mirna_sgn, ]
  mrna_use1 = mrna_use[mrna_sgn, ]
  mirna_name1 = mirna_name[mirna_sgn]
  mrna_name1 = mrna_name[mrna_sgn]
  mrna_name_sp1 = mrna_name_sp[mrna_sgn, ]
  mrna_fullname1 = mrna_fullname[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_name_sp1, mrna_fullname1))
}
mirna_mrna_data = function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_sp, mrna_fullname, cutoff){
  mirna_use[is.na(mirna_use)] = 0
  mrna_use[is.na(mrna_use)] = 0
  mirna_use[mirna_use < 0] = 0
  mrna_use[mrna_use < 0] = 0
  mirna_sgn = seq(1, nrow(mirna_use), 1)
  mrna_sgn = seq(1, nrow(mrna_use), 1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,] == 0) > ncol(mirna_use) * cutoff){
      mirna_sgn[i] = 0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,] == 0) > ncol(mrna_use) * cutoff){
      mrna_sgn[i] = 0
    }
  }  
  mirna_use1 = mirna_use[mirna_sgn,]
  mrna_use1 = mrna_use[mrna_sgn,]
  mirna_name1 = mirna_name[mirna_sgn]
  mrna_name1 = mrna_name[mrna_sgn]
  mrna_name_sp1 = mrna_name_sp[mrna_sgn,]
  mrna_fullname1 = mrna_fullname[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_name_sp1, mrna_fullname1))
}
#MMI_location is used to get the position of the validation set in the predicted matrix
MMI_location = function(Vset, mirna_name, mrna_name){
  Vset_loc = c(0, nrow(Vset))
  temp1 = as.numeric(mrna_name[, 2])
  temp2 = match(Vset[, 1], temp1)
  temp3 = match(Vset[, 2], mirna_name)
  for(i in 1 : nrow(Vset)){
    if(is.na(temp2[i]) == FALSE && is.na(temp3[i]) == FALSE){
      Vset_loc[i] = length(mirna_name) * (temp2[i] - 1) + temp3[i]
    }
  }  
  return(Vset_loc)
}
#mirna_mrna_loc is used to reshape the wMRE (or qMRE) matrix to fit the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, wMRE, mrna_fullname){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  wMRE = t(wMRE)
  wMRE_use = wMRE[mirna_loc, mrna_loc]
  rownames(wMRE_use) = mirna_name
  colnames(wMRE_use) = mrna_fullname
  return(list(mirna_loc, mrna_loc, wMRE_use))
}
#sum_mirtigo is used to calculate the normalization coefficient of mirtigo algorithm
sum_mirtigo = function(mirna, mrna, wMRE){
  res = rep(0, ncol(mirna))
  for(i in 1:ncol(mirna)){
    mirna_use1 = as.numeric(mirna[, i])
    mrna_use1 = as.numeric(mrna[, i])
    temp = mirna_use1 %*% t(mrna_use1)
    temp2 = wMRE*temp
    res[i] = sum(temp2)
  }
  return(res)
}
#mirtigo_enrich is used to get the enrichment result of mirtigo algorithm
mirtigo_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, wMRE, sumup){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * wMRE / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
  }
  thres_num = floor(length(which(pro_matrix_agg != 0))/100)
  order_agg = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  value_agg = sort(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_agg/nrow(wMRE))
  mirna_list = order_agg - (mrna_list - 1) * nrow(wMRE)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_agg[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_agg[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) #this is for mirna total probability
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) #this is for mrna total probability
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] #this is for mirna names
  mrna_ranklist = mrna_name[mrna_rank] #this is for mrna names
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}
#promise_enrich is used to get the enrichment result of promise conserved mrna
promise_enrich = function(mirna_use1, mrna_use1, mirna_name, mrna_name, wMRE){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_matrix_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  wMRE_new = t(wMRE)
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(wMRE_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(wMRE_new)
    rs = roleswitch(x, z, wMRE_new)
    pro_matrix_temp = t(rs$p.x)
    pro_matrix_mrna = pro_matrix_mrna + pro_matrix_temp
  }
  thres_num = floor(length(which(pro_matrix_mrna != 0))/100)
  order_mrna = order(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
  value_mrna = sort(pro_matrix_mrna, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_mrna/nrow(wMRE))
  mirna_list = order_mrna - (mrna_list - 1) * nrow(wMRE)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_mrna[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_mrna[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) #this is for mirna total probability
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) #this is for mrna total probability
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] #this is for mirna names
  mrna_ranklist = mrna_name[mrna_rank] #this is for mrna names
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum))
}
#targetscan_enrich is used to get the enrichment result of targetscan
targetscan_enrich = function(wMRE, mirna_name, mrna_name){
  thres_num = floor(length(which(wMRE != 0))/100)
  order_tar = order(wMRE, decreasing = TRUE)[1 : thres_num]
  value_tar = sort(wMRE, decreasing = TRUE)[1 : thres_num]
  mrna_list = ceiling(order_tar/nrow(wMRE))
  mirna_list = order_tar - (mrna_list - 1) * nrow(wMRE)
  mirna_unique = unique(mirna_list)
  mrna_unique = unique(mrna_list)
  mirna_total_sum = rep(0, length(mirna_unique))
  mrna_total_sum = rep(0, length(mrna_unique))
  for(k in 1:length(mirna_unique)){
    mirna_total_sum[k] = sum(value_tar[which(mirna_list == mirna_unique[k])])
  }
  for(k in 1:length(mrna_unique)){
    mrna_total_sum[k] = sum(value_tar[which(mrna_list == mrna_unique[k])])
  }
  mirna_ranksum = sort(mirna_total_sum, decreasing = TRUE) #this is for mirna total probability
  mirna_rank = mirna_unique[order(mirna_total_sum, decreasing = TRUE)]
  mrna_ranksum = sort(mrna_total_sum, decreasing = TRUE) #this is for mrna total probability
  mrna_rank = mrna_unique[order(mrna_total_sum, decreasing = TRUE)]
  mirna_ranklist = mirna_name[mirna_rank] #this is for mirna names
  mrna_ranklist = mrna_name[mrna_rank] #this is for mrna names
  return(list(mirna_ranklist, mrna_ranklist, mirna_ranksum, mrna_ranksum)) 
}

####################################################################################################

   ### Part II. INPUT DATA: INPUT FOR POPULATION-LEVEL EVALUATION ###

#the following package is used for promise method and expression-based methods
library(Roleswitch)
library(glmnet)
#first read in the files we need to use, which contains:
#input for mirtigo algorithm: mirna expression, mrna expression, CWCS_matrix, mirna list, mrna list
#input for promise algorithm: mirna expression, mrna expression, qMRE_matrix, mirna list, mrna list
#input for comparison: from oncomirs, miRsurvival and Cancer gene_2017
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
qMRE_conserved_matrix = as.matrix(read.table("conserved_qMRE.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("wMRE_all.txt", head = TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
mirna_valid1 = as.matrix(read.table("Oncomirs.txt", head = TRUE, sep = "\t"))
mirna_valid2 = as.matrix(read.table("miRNA biomarkers.txt", head = TRUE, sep = "\t"))
mrna_valid = as.matrix(read.table("Cancer genes.txt", head = TRUE, sep = "\t"))
name_cancer = as.matrix(read.table("Sample-to-sample.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("Input miRNA expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("Input mRNA expression.txt",head = FALSE, sep = "\t"))
rm(CWCS_matrix1)
gc()

####################################################################################################

   ### Part III. MAIN PROGRAM: MAIN FUNCTIONS ###

#the following code is used to do the preparation.
x = name_func(name_cancer, mirna_cancer, mrna_cancer)
z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], mirna)
z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], mrna)
z3 = mirna_mrna_data_unselected(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]])
z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
z6 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_conserved_matrix, z3[[6]])
z7 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, CWCS_matrix, z3[[6]])
z7_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, CWCS_matrix, z4[[6]])

#the following code is used to get the enriched mirna or mrna in order
z8 = sum_mirtigo(z3[[1]], z3[[2]], z7[[3]])
z9 = mirtigo_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z7[[3]], z8)
z10 = promise_enrich(z3[[1]], z3[[2]], z3[[3]], z3[[4]], z6[[3]])
z11 = targetscan_enrich(z7_tar[[3]], z4[[3]], z4[[4]])

#mirna_num and mrna_num is the overlap of validated miRNA and mRNA in the matrix
mirna_num1 = sum(mirna_valid1 %in% z3[[3]])
mirna_num2 = sum(mirna_valid2 %in% z3[[3]])
mrna_num = sum(mrna_valid[,1] %in% z3[[4]])
mirna_num1_tar = sum(mirna_valid1 %in% z4[[3]])
mirna_num2_tar = sum(mirna_valid2 %in% z4[[3]])
mrna_num_tar = sum(mrna_valid[,1] %in% z4[[4]])

#the following code is used to get the final output, mirna_fc1, mirna_fc2 and mrna_fc contains the detailed information of enrichment result
mirna_fc1 = matrix(data = 0, nrow = 100, ncol = 9, byrow = TRUE)
mirna_fc2 = matrix(data = 0, nrow = 100, ncol = 9, byrow = TRUE)
mrna_fc = matrix(data = 0, nrow = 100, ncol = 9, byrow = TRUE)

colnames(mirna_fc1) = c("mirtigo_number", "mirtigo_fc", "mirtigo_pvalue", "promise_number", "promise_fc", "promise_pvalue", 
                        "targetscan_number", "targetscan_fc", "targetscan_pvalue")
colnames(mirna_fc2) = c("mirtigo_number", "mirtigo_fc", "mirtigo_pvalue", "promise_number", "promise_fc", "promise_pvalue", 
                        "targetscan_number", "targetscan_fc", "targetscan_pvalue")
colnames(mrna_fc) = c("mirtigo_number", "mirtigo_fc", "mirtigo_pvalue", "promise_number", "promise_fc", "promise_pvalue", 
                      "targetscan_number", "targetscan_fc", "targetscan_pvalue")

a1 = matrix(data = 0, ncol = 2, nrow = 2, byrow = TRUE)
a2 = matrix(data = 0, ncol = 2, nrow = 2, byrow = TRUE)
a3 = matrix(data = 0, ncol = 2, nrow = 2, byrow = TRUE)

for(i in 1 : 100){
  mirna_fc1[i, 1] = sum(z9[[1]][1:i] %in% mirna_valid1)
  mirna_fc1[i, 2] = (sum(z9[[1]][1:i] %in% mirna_valid1)/i)/(mirna_num1/length(z3[[3]]))
  a1[2, 1] = mirna_fc1[i, 1]
  a1[2, 2] = i - a1[2, 1]
  a1[1, 1] = mirna_num1 - a1[2, 1]
  a1[1, 2] = length(z3[[3]]) - i - a1[1, 1]
  mirna_fc1[i, 3] = fisher.test(a1, alternative = "less")$p.value
  
  mirna_fc1[i, 4] = sum(z10[[1]][1:i] %in% mirna_valid1)
  mirna_fc1[i, 5] = (sum(z10[[1]][1:i] %in% mirna_valid1)/i)/(mirna_num1/length(z3[[3]]))
  a2[2, 1] = mirna_fc1[i, 4]
  a2[2, 2] = i - a2[2, 1]
  a2[1, 1] = mirna_num1 - a2[2, 1]
  a2[1, 2] = length(z3[[3]]) - i - a2[1, 1]
  mirna_fc1[i, 6] = fisher.test(a2, alternative = "less")$p.value
  
  mirna_fc1[i, 7] = sum(z11[[1]][1:i] %in% mirna_valid1)
  mirna_fc1[i, 8] = (sum(z11[[1]][1:i] %in% mirna_valid1)/i)/(mirna_num1_tar/length(z4[[3]]))
  a3[2, 1] = mirna_fc1[i, 7]
  a3[2, 2] = i - a3[2, 1]
  a3[1, 1] = mirna_num1_tar - a3[2, 1]
  a3[1, 2] = length(z4[[3]]) - i - a3[1, 1]
  mirna_fc1[i, 9] = fisher.test(a3, alternative = "less")$p.value
}
for(i in 1 : 100){
  mirna_fc2[i, 1] = sum(z9[[1]][1:i] %in% mirna_valid2)
  mirna_fc2[i, 2] = (sum(z9[[1]][1:i] %in% mirna_valid2)/i)/(mirna_num2/length(z3[[3]]))
  a1[2, 1] = mirna_fc2[i, 1]
  a1[2, 2] = i - a1[2, 1]
  a1[1, 1] = mirna_num2 - a1[2, 1]
  a1[1, 2] = length(z3[[3]]) - i - a1[1, 1]
  mirna_fc2[i, 3] = fisher.test(a1, alternative = "less")$p.value

  mirna_fc2[i, 4] = sum(z10[[1]][1:i] %in% mirna_valid2)
  mirna_fc2[i, 5] = (sum(z10[[1]][1:i] %in% mirna_valid2)/i)/(mirna_num2/length(z3[[3]]))
  a2[2, 1] = mirna_fc2[i, 4]
  a2[2, 2] = i - a2[2, 1]
  a2[1, 1] = mirna_num2 - a2[2, 1]
  a2[1, 2] = length(z3[[3]]) - i - a2[1, 1]
  mirna_fc2[i, 6] = fisher.test(a2, alternative = "less")$p.value

  mirna_fc2[i, 7] = sum(z11[[1]][1:i] %in% mirna_valid2)
  mirna_fc2[i, 8] = (sum(z11[[1]][1:i] %in% mirna_valid2)/i)/(mirna_num2_tar/length(z4[[3]]))
  a3[2, 1] = mirna_fc2[i, 7]
  a3[2, 2] = i - a3[2, 1]
  a3[1, 1] = mirna_num2_tar - a3[2, 1]
  a3[1, 2] = length(z4[[3]]) - i - a3[1, 1]
  mirna_fc2[i, 9] = fisher.test(a3, alternative = "less")$p.value
}
for(i in 1 : 100){
  mrna_fc[i, 1] = sum(z9[[2]][1:i] %in% mrna_valid[,1])
  mrna_fc[i, 2] = (sum(z9[[2]][1:i] %in% mrna_valid[,1])/i)/(mrna_num/length(z3[[4]]))
  a1[2, 1] = mrna_fc[i, 1]
  a1[2, 2] = i - a1[2, 1]
  a1[1, 1] = mrna_num - a1[2, 1]
  a1[1, 2] = length(z3[[4]]) - i - a1[1, 1]
  mrna_fc[i, 3] = fisher.test(a1, alternative = "less")$p.value

  mrna_fc[i, 4] = sum(z10[[2]][1:i] %in% mrna_valid[,1])
  mrna_fc[i, 5] = (sum(z10[[2]][1:i] %in% mrna_valid[,1])/i)/(mrna_num/length(z3[[4]]))
  a2[2, 1] = mrna_fc[i, 4]
  a2[2, 2] = i - a2[2, 1]
  a2[1, 1] = mrna_num - a2[2, 1]
  a2[1, 2] = length(z3[[4]]) - i - a2[1, 1]
  mrna_fc[i, 6] = fisher.test(a2, alternative = "less")$p.value

  mrna_fc[i, 7] = sum(z11[[2]][1:i] %in% mrna_valid[,1])
  mrna_fc[i, 8] = (sum(z11[[2]][1:i] %in% mrna_valid[,1])/i)/(mrna_num_tar/length(z4[[4]]))
  a3[2, 1] = mrna_fc[i, 7]
  a3[2, 2] = i - a3[2, 1]
  a3[1, 1] = mrna_num_tar - a3[2, 1]
  a3[1, 2] = length(z4[[4]]) - i - a3[1, 1]
  mrna_fc[i, 9] = fisher.test(a3, alternative = "less")$p.value
}

write.table(mirna_fc1,file = "mirna fold change-oncomirs.txt", quote = FALSE, sep = "\t")
write.table(mirna_fc2,file = "mirna fold change-miRNA biomarkers.txt", quote = FALSE, sep = "\t")
write.table(mrna_fc,file = "mrna fold change-Cancer genes.txt", quote = FALSE, sep = "\t")
