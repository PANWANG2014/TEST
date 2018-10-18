######################################################################################################

   ### Part I. FUNCTIONS: FUNCTIONS THAT ARE USED IN DATA PROCESSING AND MODELING ###

#The following three functions are used to match the names of samples, mirnas and mrnas
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
  mrna_exist = rep(0, nrow(mrna_base))
  mrna_name1 = rep(0, nrow(mrna_base))
  for(i in 1 : nrow(mrna_base)){
    mrna_name1[i] = unlist(strsplit(mrna_base[i, 1], "\\|"))[1]
    if(mrna_name1[i] %in% mrna == 1){
      mrna_exist[i] = i
    }
  }
  mrna_use = mrna_base[mrna_exist, mrna_name]
  mrna_name2 = mrna_name1[mrna_exist]
  mrna_fullname = mrna_base[mrna_exist, 1]
  return(list(mrna_use, mrna_name2, mrna_fullname))
}

#mirna_mrna_data_unselected will select the RNAs that are expressed in at least one sample
mirna_mrna_data_unselected = function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_fullname){
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
  mrna_fullname1 = mrna_fullname[mrna_sgn]
  rownames(mirna_use1) = mirna_name1
  rownames(mrna_use1) = mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_fullname1))
}
#mirna_mrna_loc is used to get the targetscan matrix that fits the expression data
mirna_mrna_loc = function(mirna_name, mrna_name, mirna, mrna, wMRE, mrna_fullname){
  mirna_loc = match(mirna_name, mirna)
  mrna_loc = match(mrna_name, mrna)
  wMRE = t(wMRE)
  wMRE_use = wMRE[mirna_loc, mrna_loc]
  rownames(wMRE_use) = mirna_name
  colnames(wMRE_use) = mrna_fullname
  return(list(mirna_loc, mrna_loc, wMRE_use))
}

#the following function sum_mirtigo and mirtigo_detail are the core functions of the mirtigo algorithm.
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
#mirtigo_basic is used to get the detailed information of the top rankers by mirtigo algorithm, both cancer level and sample level
#Here, we offer two sets of results.
#One is default, which is the first percent pairs, the other is the number of top rankers which is set by the users
mirtigo_detail = function(mirna_use1, mrna_use1, wMRE, sumup, ID, mirna_name, mrna_name, thres_num){
  numm = floor(length(wMRE[wMRE != 0]))/100
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  #pro_matrix_temp saves the probability matrix of each single sample, while pro_matrix_agg is the aggregated result
  pro_matrix_temp = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  pro_matrix_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  
  detail_output1 = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = thres_num, byrow = TRUE)
  detail_output2 = matrix(data = 0, ncol = (4 * ncol(mirna_use)), nrow = numm, byrow = TRUE)
  temp_name1 = c()
  temp_name2 = c()
  num_max = 0
  
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * wMRE / sumup[m]
    pro_matrix_agg = pro_matrix_agg + pro_matrix_temp
    sort_one_rank1 = order(pro_matrix_temp, decreasing = TRUE)[1 : thres_num]
    sort_one_prob1 = sort(pro_matrix_temp, decreasing = TRUE)[1 : thres_num]
    temp1 = ceiling(sort_one_rank1/length(mirna_name))
    temp2 = sort_one_rank1 - (temp1 - 1) * length(mirna_name)
    for(i in 1 : thres_num){
      detail_output1[i, (4 * m - 3)] = i
      detail_output1[i, (4 * m - 2)] = mrna_name[temp1[i]]
      detail_output1[i, (4 * m - 1)] = mirna_name[temp2[i]]
      detail_output1[i, (4 * m)] = sort_one_prob1[i]
    }
    temp_name1 = c(temp_name1, ID[m], "mrna name", "mirna name", "mirtigo score")
    
    num1 = floor(length(pro_matrix_temp[pro_matrix_temp != 0]))/100
    sort_one_rank2 = order(pro_matrix_temp, decreasing = TRUE)[1 : num1]
    sort_one_prob2 = sort(pro_matrix_temp, decreasing = TRUE)[1 : num1]
    temp3 = ceiling(sort_one_rank2/length(mirna_name))
    temp4 = sort_one_rank2 - (temp3 - 1) * length(mirna_name)
    for(i in 1 : num1){
      detail_output2[i, (4 * m - 3)] = i
      detail_output2[i, (4 * m - 2)] = mrna_name[temp3[i]]
      detail_output2[i, (4 * m - 1)] = mirna_name[temp4[i]]
      detail_output2[i, (4 * m)] = sort_one_prob2[i]
    }
    temp_name2 = c(temp_name2, ID[m], "mrna name", "mirna name", "mirtigo score")
    num_max = max(num_max, num1)
  }
  colnames(detail_output1) = temp_name1
  colnames(detail_output2) = temp_name2
  detail_output3 = detail_output2[1 : num_max, ]
  
  num2 = floor(length(pro_matrix_agg[pro_matrix_agg != 0]))/100
  out_matrix1 = matrix(data = 0, ncol = 4, nrow = thres_num, byrow = TRUE)
  out_matrix2 = matrix(data = 0, ncol = 4, nrow = num2, byrow = TRUE)
  
  temp_rank1 = order(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  temp_prob1 = sort(pro_matrix_agg, decreasing = TRUE)[1 : thres_num]
  temp_rank2 = order(pro_matrix_agg, decreasing = TRUE)[1 : num2]
  temp_prob2 = sort(pro_matrix_agg, decreasing = TRUE)[1 : num2]

  mrna1 = ceiling(temp_rank1/length(mirna_name))
  mirna1 = temp_rank1 - (mrna1 - 1) * length(mirna_name)
  
  mrna2 = ceiling(temp_rank2/length(mirna_name))
  mirna2 = temp_rank2 - (mrna2 - 1) * length(mirna_name)
  
  for(i in 1:nrow(out_matrix1)){
    out_matrix1[i, 1] = i
    out_matrix1[i, 2] = mrna_name[mrna1[i]]
    out_matrix1[i, 3] = mirna_name[mirna1[i]]
    out_matrix1[i, 4] = temp_prob1[i]
  }
  for(i in 1:nrow(out_matrix2)){
    out_matrix2[i, 1] = i
    out_matrix2[i, 2] = mrna_name[mrna2[i]]
    out_matrix2[i, 3] = mirna_name[mirna2[i]]
    out_matrix2[i, 4] = temp_prob2[i]
  }
  colnames(out_matrix1) = c("rank", "mrna name", "mirna name", "mirtigo score")
  colnames(out_matrix2) = c("rank", "mrna name", "mirna name", "mirtigo score")
  
  return(list(out_matrix1, out_matrix2, detail_output1, detail_output3))
}

#####################################################################################################

   ### Part II. INPUT DATA: RAW DATA FOR THE MODEL ###

#The input data contains two parts.
#The first part is to provide sequence information
#Here, CWCS_matrix is got from Targetscan website, which is a matrix of score of each pair of mRNA and miRNA
#Also, mrna and mirna are the list that can be predicted by the model, which is also got from Targetscan website
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("wMRE_all.txt", head = TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
rm(CWCS_matrix1)
gc()

#The second part is to provide expression data
#mirna_cancer and mrna_cancer are the mirna and mrna expression data, whose format should be strictly satisfied
#name_cancer provides the relationship of sample names in mirna expression and sample names in mrna expression
name_cancer = as.matrix(read.table("Sample-to-sample file.txt", head = TRUE, sep = "\t"))
mirna_cancer = as.matrix(read.table("Input miRNA expression.txt", head = FALSE, sep = "\t"))
mrna_cancer = as.matrix(read.table("Input mRNA expression.txt", head = FALSE, sep = "\t"))
#thres_num is the number of top rankers that users are interested in.
thres_num = 5000
#####################################################################################################

   ### Part III. MAIN PROGRAM: MAIN FUNCTION ###

#Our algorithm provides both population-level result and sample-level result
x = name_func(name_cancer, mirna_cancer, mrna_cancer)
z1 = mirna_matrix(name_cancer, mirna_cancer, x[[1]], mirna)
z2 = mrna_matrix(name_cancer, mrna_cancer, x[[2]], mrna)
z3 = mirna_mrna_data_unselected(z1[[1]], z2[[1]], z1[[2]], z2[[2]], z2[[3]])
z7 = mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, CWCS_matrix, z3[[5]])
z8 = sum_mirtigo(z3[[1]], z3[[2]], z7[[3]])
z9 = mirtigo_detail(z3[[1]], z3[[2]], z7[[3]], z8, name_cancer[,1], z3[[3]], z3[[5]], thres_num)

write.table(z9[[1]],file = paste("top", thres_num, "pairs of population-level miRTIGO algorithm.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z9[[2]],file = "top 1 percent pairs of population-level miRTIGO algorithm.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z9[[3]],file = paste("top", thres_num, "pairs of sample-level miRTIGO algorithm.txt", sep = " "), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z9[[4]],file = "top 1 percent pairs of sample-level miRTIGO algorithm.txt", quote = FALSE, sep = "\t", row.names = FALSE)
