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

#mirtigo_NCI is specially designed to get the result of NCI-60 data, since each cell-line only have very few samples
mirtigo_NCI = function(mirna_use1, mrna_use1, wMRE, sumup, thres_num, rep_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  cel_num = ncol(mirna_use)/rep_num
  sort_agg = matrix(data = 0, ncol = cel_num, nrow = thres_num, byrow = TRUE)
  for(m in 1 : cel_num){
    pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
    for(k in 1 : rep_num){
      j = rep_num * (m - 1) + k
      pro_matrix_temp = (mirna_use[, j] %*% t(mrna_use[, j])) * wMRE/sumup[j]
      pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
    }
    sort_agg[, m] = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  }
  return(sort_agg)
}
#promise_NCI is used to calculate the result of promise for NCI-60 data
promise_NCI = function(mirna_use1, mrna_use1, wMRE, thres_num, rep_num){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  cel_num = ncol(mirna_use)/rep_num
  
  sort_agg_mrna = matrix(data = 0, ncol = cel_num, nrow = thres_num, byrow = TRUE)
  sort_agg_full = matrix(data = 0, ncol = cel_num, nrow = thres_num, byrow = TRUE)
  wMRE_new = t(wMRE)
  for(m in 1 : cel_num){
    pro_agg_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
    pro_agg_full = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
    for(k in 1 : rep_num){
      j = rep_num * (m - 1) + k
      x = matrix(mrna_use[, j])
      rownames(x) = rownames(wMRE_new)
      z = matrix(mirna_use[, j])
      rownames(z) = colnames(wMRE_new)
      rs = roleswitch(x, z, wMRE_new)
      pro_matrix_mrna_temp = t(rs$p.x)
      pro_matrix_full_temp = t(rs$p.xz)
      pro_agg_mrna = pro_agg_mrna + pro_matrix_mrna_temp
      pro_agg_full = pro_agg_full + pro_matrix_full_temp
    }
    sort_agg_mrna[, m] = order(pro_agg_mrna, decreasing = TRUE)[1 : thres_num]
    sort_agg_full[, m] = order(pro_agg_full, decreasing = TRUE)[1 : thres_num]
  }
  return(list(sort_agg_mrna, sort_agg_full))
}
#targetscan_rank is used to get the result of targetscan for NCI-60 data
targetscan_NCI = function(wMRE, thres_num){
  target_rank = order(wMRE, decreasing = TRUE)[1 : thres_num]
  return(target_rank)
}

#result_output_NCI is used to get the comparison result of mirtigo and promise with the validation set
result_output_NCI = function(re_matrix, Vset){
  comp_res = matrix(data = 0,ncol = ncol(re_matrix), nrow = 5000, byrow = TRUE)
  comp_out = matrix(data = 0,nrow = ncol(re_matrix), ncol = 50, byrow = TRUE)
  for(i in 1:ncol(re_matrix)){
    comp_res[, i] = re_matrix[, i] %in% Vset
    for(j in 1:50){
      comp_out[i,j] = sum(comp_res[1:(100*j), i])/(100 * j)
    }
  }
  return(list(comp_res, comp_out))
}
#result_output_tar is used to get the comparison result of targetscan with the validation set
result_output_tar = function(rank_matrix, Vset){
  Vset = as.numeric(Vset)
  rank = as.numeric(rank_matrix)
  rank_result = rank %in% Vset
  rank_str = rep(0, 50)
  rank_str_percent = rep(0, 50)
  for(i in 1:50){
    temp1 = 100 * i
    rank_str[i] = sum(rank_result[1:temp1])
    rank_str_percent[i] = rank_str[i]/temp1
  }
  return(list(rank_str, rank_str_percent))
}

####################################################################################################

   ### Part II. INPUT DATA: INPUT FOR POPULATION-LEVEL EVALUATION ###

#the following package is used for promise method
library(Roleswitch)
#first read in the files we need to use, which contains:
#input for mirtigo algorithm: mirna expression, mrna expression, CWCS_matrix, mirna list, mrna list
#input for promise algorithm: mirna expression, mrna expression, qMRE_matrix, mirna list, mrna list
#input for validation: Vset
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
qMRE_matrix = as.matrix(read.table("qMRE_all.txt", head = TRUE, sep = "\t"))
qMRE_conserved_matrix = as.matrix(read.table("conserved_qMRE.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("wMRE_all.txt", head=TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
Vset = read.table("V1.txt", head = TRUE, sep = "\t")
name_cancer = as.matrix(read.table("Sample-to-sample.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("Input miRNA expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("Input mRNA expression.txt", head = FALSE, sep = "\t"))
rm(CWCS_matrix1)
gc()

####################################################################################################

   ### Part III. MAIN PROGRAM: MAIN FUNCTIONS ###

#the following code is used to do the preparation.
thres_num = 5000
rep_num = 2 # for cell lines with 4 duplicates, this parameter should be set 4
x = name_func(name_cancer, mirna_cancer, mrna_cancer)
z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], mirna)
z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], mrna)
z3 = mirna_mrna_data_unselected(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]])
z4 = mirna_mrna_data(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]], z2[[4]], 0.8)
z5 = MMI_location(Vset, z3[[3]], z3[[5]]) #this is used to locate the validation set
z5_tar = MMI_location(Vset, z4[[3]], z4[[5]]) #this is used to locate the validation set of selected targetscan
z6 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_matrix, z3[[6]])
z6_con = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, qMRE_conserved_matrix, z3[[6]])
z7 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, CWCS_matrix, z3[[6]])
z6_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, qMRE_matrix, z4[[6]])
z7_tar = mirna_mrna_loc(z4[[3]], z4[[4]], mirna, mrna, CWCS_matrix, z4[[6]])

#the following code is used to calculate the top rankers from mirtigo, promise, and targetscan
z8 = sum_mirtigo(z3[[1]], z3[[2]], z7[[3]])
z9 = mirtigo_NCI(z3[[1]], z3[[2]], z7[[3]], z8,  thres_num, rep_num) #get the rankers of mirtigo
z10 = promise_NCI(z3[[1]], z3[[2]], z6[[3]], thres_num, rep_num) #get the rankers of promise_qMRE
z11 = promise_NCI(z3[[1]], z3[[2]], z6_con[[3]], thres_num, rep_num) #get the rankers of promise_qMRE_conserve
z12 = targetscan_NCI(z6_tar[[3]], thres_num) #get the rankers of targetscan_qMRE
z13 = targetscan_NCI(z7_tar[[3]], thres_num) #get the rankers of targetscan_wMRE

z20 = result_output_NCI(z9, z5) #z20 is the result of mirtigo
z21 = result_output_NCI(z10[[1]], z5) #z21 is the result of promise mrna
z22 = result_output_NCI(z10[[2]], z5) #z22 is the result of promise full
z23 = result_output_NCI(z11[[1]], z5) #z23 is the result of promise_conserve mrna
z24 = result_output_NCI(z11[[2]], z5) #z24 is the result of promise_conserve full
z25 = result_output_tar(z12, z5_tar) #z25 is the result of targetscan_qMRE
z26 = result_output_tar(z13, z5_tar) #z26 is the result of targetscan_wMRE

write.table(z20[[2]],file="comparison result_mirtigo.txt",quote=FALSE,sep="\t")
write.table(z21[[2]],file="comparison result_promise mrna.txt",quote=FALSE,sep="\t")
write.table(z22[[2]],file="comparison result_promise full.txt",quote=FALSE,sep="\t")
write.table(z23[[2]],file="comparison result_promise_conserve mrna.txt",quote=FALSE,sep="\t")
write.table(z24[[2]],file="comparison result_promise_conserve full.txt",quote=FALSE,sep="\t")
write.table(z25[[2]],file="comparison result_targetscan qMRE.txt",quote=FALSE,sep="\t")
write.table(z26[[2]],file="comparison result_targetscan wMRE.txt",quote=FALSE,sep="\t")
