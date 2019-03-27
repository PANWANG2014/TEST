library(caret)
#the first part is to load the data, including mirna expression, mrna expression, and wMRE matrix
#basic_information contains some information about patients and samples
basic_information = as.matrix(read.table("basic information.txt", head = TRUE, sep = "\t"))
#This is the total mirna and mrna list and wMRE matrix
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
wMRE_matrix_cwcs1 = as.matrix(read.table("wMRE_all.txt", head = TRUE, sep = "\t"))
wMRE_matrix_cwcs = abs(wMRE_matrix_cwcs1)
#This is used to get the CICAMS matrix
mirna_CICAMS_T = as.matrix(read.table("CICAMS-ESCC-miRNA-agilent.txt", head = FALSE, sep = "\t"))
mrna_CICAMS_T = as.matrix(read.table("CICAMS-ESCC-mRNA-agilent.txt", head = FALSE, sep = "\t"))
name_CICAMS_T = as.matrix(read.table("ESCC_T_1to1.txt", head = TRUE, sep = "\t"))
#This is used to get the TCGA matrix
mirna_CICAMS_E = as.matrix(read.table("CICAMS-E-miRNA-agilent.txt", head = FALSE, sep = "\t"))
mrna_CICAMS_E = as.matrix(read.table("CICAMS-E-mRNA-agilent.txt", head = FALSE, sep = "\t"))
name_CICAMS_E = as.matrix(read.table("ESCC_E_1to1.txt", head = TRUE, sep = "\t"))
#z3_E is expression matrix for E file
x_E = name_func(name_CICAMS_E, mirna_CICAMS_E, mrna_CICAMS_E)
z1_E = mirna_matrix(name_CICAMS_E, mirna_CICAMS_E, x_E[[1]], mirna)
z2_E = mrna_matrix(name_CICAMS_E, mrna_CICAMS_E, x_E[[2]], mrna)
z3_E = mirna_mrna_data_unselected(z1_E[[1]], z2_E[[1]], z1_E[[2]], z2_E[[2]], z2_E[[3]], z2_E[[4]])
#z3_T is expression matrix for T file
x_T = name_func(name_CICAMS_T, mirna_CICAMS_T, mrna_CICAMS_T)
z1_T = mirna_matrix(name_CICAMS_T, mirna_CICAMS_T, x_T[[1]], mirna)
z2_T = mrna_matrix(name_CICAMS_T, mrna_CICAMS_T, x_T[[2]], mrna)
z3_T = mirna_mrna_data_unselected(z1_T[[1]], z2_T[[1]], z1_T[[2]], z2_T[[2]], z2_T[[3]], z2_T[[4]])

#first match the order of expression data with basic_information, and split the training data from test data
times = 20
part = cate_part(basic_information, times) #this partition is based on the order of basic_information
compare = match_part(name_CICAMS_E, basic_information)#name_CICAMS[compare[1],] is basic_information[1,]
test_record = matrix(unlist(part[[2]]), nrow = length(part[[2]][[1]]))
t = seq(1, length(part[[1]][[1]]), 1)
part_strategy = matrix(data = 0, ncol = length(part[[1]][[1]]), nrow = times)
#the following iteration is used to split the data
for(time in 1:times){
  s = sample(t, size = ceiling(0.1 * ncol(part_strategy)), replace = FALSE)
  part_strategy[time, s] = 1
  c = s
  for(i in 2:7){
    s = sample(t[-c], size = ceiling(0.1 * ncol(part_strategy)), replace = FALSE)
    c = c(c,s)
    part_strategy[time, s] = i
  }
  for(i in 8:9){
    s = sample(t[-c], size = (ceiling(0.1 * ncol(part_strategy)) - 1), replace = FALSE)
    c = c(c,s)
    part_strategy[time, s] = i
  }
  s = t[-c]
  part_strategy[time, s] = 10
}
top_feature = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000)
#feature_number is used to save the features selected in each cross validation
feature_number_mrna = matrix(data = 0, ncol = 10, nrow = times)
#cross is used to record the accuracy among iterations when we choose different numbers of top_features
#we can use these numbers to calculate the various measurements
cross_sep_T = matrix(data = 0, ncol = times, nrow = length(top_feature))
cross_nor_T = matrix(data = 0, ncol = times, nrow = length(top_feature))
cross_sep_E = matrix(data = 0, ncol = times, nrow = length(top_feature))
cross_nor_E = matrix(data = 0, ncol = times, nrow = length(top_feature))
total_val = matrix(data = 0, ncol = times, nrow = length(top_feature))
feature_detail_mrna = list()
for(i in 1:times){
  feature_detail_mrna[[i]] = matrix(data = 0, ncol = 10, nrow = 1000)
}
#feature_vali_number is used to record the number of features in each time of validation step
feature_vali_number = matrix(data = 0, ncol = times, nrow = length(top_feature))
accu_sep_T = matrix(data = 0, ncol = times, nrow = length(top_feature))
accu_nor_T = matrix(data = 0, ncol = times, nrow = length(top_feature))
accu_sep_E = matrix(data = 0, ncol = times, nrow = length(top_feature))
accu_nor_E = matrix(data = 0, ncol = times, nrow = length(top_feature))
total_val_vali = matrix(data = 0, ncol = times, nrow = length(top_feature))

#the following part is used to do the cross validation part:
for(time in 1:times){
  #first split the training data and validation data
  mrna_E = z3_E[[2]][,compare[part[[1]][[time]]]]
  mrna_T = z3_T[[2]][,compare[part[[1]][[time]]]]
  vali_mrna_E = z3_E[[2]][,compare[part[[2]][[time]]]]
  vali_mrna_T = z3_T[[2]][,compare[part[[2]][[time]]]]
  
  #the following code is used to do the cross validation and select the features from mirna
  for(k in 1:10){
    temp2 = which(part_strategy[time,] == k) #for test
    temp1 = t[-temp2] #for training
    #the following code is used to find the features selected by t test and expression mean difference
    pos_mrna = feature_select(mrna_E, mrna_T, temp1) #pos_mrna is the ranked features
    feature_number_mrna[time, k] = length(pos_mrna)
    #feature_detail_mirna and feature_detail_mrna save the top features of mirna or mrna
    feature_detail_mrna[[time]][, k] = pos_mrna[1:1000]
    #the following is to get the parameter for naive bayes, including normal version and seperate version
    para_sep_mrna = para_sep(mrna_E, mrna_T, temp1, pos_mrna[1:1000])
    para_E_nor_mrna = para_nor(mrna_E, temp1, pos_mrna[1:1000])
    para_T_nor_mrna = para_nor(mrna_T, temp1, pos_mrna[1:1000])
    
    #the following is to finish the cross validation part for different parameter
    for(j in 1:length(top_feature)){
      prob_E_sep = test_sep(para_sep_mrna, mrna_E, temp2, pos_mrna[1:top_feature[j]])
      prob_T_sep = test_sep(para_sep_mrna, mrna_T, temp2, pos_mrna[1:top_feature[j]])
      prob_E_nor = test_nor(para_E_nor_mrna, para_T_nor_mrna, mrna_E, temp2, pos_mrna[1:top_feature[j]])
      prob_T_nor = test_nor(para_E_nor_mrna, para_T_nor_mrna, mrna_T, temp2, pos_mrna[1:top_feature[j]])  
      cross_sep_E[j, time] = cross_sep_E[j, time] + sum(prob_E_sep[1,] >= prob_E_sep[2,])
      cross_sep_T[j, time] = cross_sep_T[j, time] + sum(prob_T_sep[1,] < prob_T_sep[2,])
      cross_nor_E[j, time] = cross_nor_E[j, time] + sum(prob_E_nor[1,] >= prob_E_nor[2,])
      cross_nor_T[j, time] = cross_nor_T[j, time] + sum(prob_T_nor[1,] < prob_T_nor[2,])
      total_val[j, time] = total_val[j, time] + length(temp2)
    }
  }
  
  #the following part is used to do the internal independent validation
  temp3 = seq(1, ncol(mrna_E), 1)
  temp4 = seq(1, ncol(vali_mrna_E), 1)
  for(j in 1:length(top_feature)){
    #get the mrna features
    fea_mrna = feature_detail_mrna[[time]][1:top_feature[j], 1]
    for(k in 1:10){
      fea_mrna = intersect(fea_mrna, feature_detail_mrna[[time]][1:top_feature[j], k])
    }
    #train the parameter
    vali_sep_mrna = para_sep(mrna_E, mrna_T, temp3, fea_mrna)
    vali_E_nor_mrna = para_nor(mrna_E, temp3, fea_mrna)
    vali_T_nor_mrna = para_nor(mrna_T, temp3, fea_mrna)
    #get the probability
    prob_E_sep_vali = test_sep(vali_sep_mrna, vali_mrna_E, temp4, fea_mrna)
    prob_T_sep_vali = test_sep(vali_sep_mrna, vali_mrna_T, temp4, fea_mrna)
    prob_E_nor_vali = test_nor(vali_E_nor_mrna, vali_T_nor_mrna, vali_mrna_E, temp4, fea_mrna)
    prob_T_nor_vali = test_nor(vali_E_nor_mrna, vali_T_nor_mrna, vali_mrna_T, temp4, fea_mrna)
    #record the accuracy
    accu_sep_E[j, time] = accu_sep_E[j, time] + sum(prob_E_sep_vali[1,] >= prob_E_sep_vali[2,])
    accu_sep_T[j, time] = accu_sep_T[j, time] + sum(prob_T_sep_vali[1,] < prob_T_sep_vali[2,])
    accu_nor_E[j, time] = accu_nor_E[j, time] + sum(prob_E_nor_vali[1,] >= prob_E_nor_vali[2,])
    accu_nor_T[j, time] = accu_nor_T[j, time] + sum(prob_T_nor_vali[1,] < prob_T_nor_vali[2,])
    total_val_vali[j, time] = total_val_vali[j, time] + ncol(vali_mrna_E)
  }
  print(time)
}


summ_result = matrix(data = 0, ncol = 22, nrow = length(top_feature))
colnames(summ_result) = c("chosen feature number", "intersect feature number", "union feature number", "jaccard index",
                          "ACC using normal in CV", "ACC using sep in CV", "ACC using normal in IV", "ACC using sep in IV",
                          "PPV using normal in CV", "PPV using sep in CV", "PPV using normal in IV", "PPV using sep in IV",
                          "TPR using normal in CV", "TPR using sep in CV", "TPR using normal in IV", "TPR using sep in IV",
                          "TNR using normal in CV", "TNR using sep in CV", "TNR using normal in IV", "TNR using sep in IV",
                          "F score using normal in IV", "F score using sep in IV")

temp_inter = list()
temp_uni = list()
for(j in 1:length(top_feature)){
  temp_inter[[j]] =  feature_detail_mrna[[1]][1:top_feature[j],1]
  temp = feature_detail_mrna[[1]][1:top_feature[j], 1]
  for(k in 1:10){
    temp = intersect(temp, feature_detail_mrna[[time]][1:top_feature[j], k])
  }
  temp_uni[[j]] = temp
  for(time in 1:times){
    temp = feature_detail_mrna[[time]][1:top_feature[j], 1]
    for(k in 1:10){
      temp = intersect(temp, feature_detail_mrna[[time]][1:top_feature[j], k])
    }
    temp_inter[[j]] = intersect(temp_inter[[j]], temp)
    temp_uni[[j]] = union(temp_uni[[j]], temp)
  }
}
jaccard_mrna = rep(0, length(top_feature))
for(i in 1:length(top_feature)){
  summ_result[i, 2] = length(temp_inter[[i]])
  summ_result[i, 3] = length(temp_uni[[i]])
  summ_result[i, 4] = length(temp_inter[[i]])/length(temp_uni[[i]])
}
summ_result[, 1] = top_feature
#ACC
summ_result[, 5] = (rowSums(cross_nor_E) + rowSums(cross_nor_T))/(rowSums(total_val) + rowSums(total_val))
summ_result[, 6] = (rowSums(cross_sep_E) + rowSums(cross_sep_T))/(rowSums(total_val) + rowSums(total_val))
summ_result[, 7] = (rowSums(accu_nor_E) + rowSums(accu_nor_T))/(rowSums(total_val_vali) + rowSums(total_val_vali))
summ_result[, 8] = (rowSums(accu_sep_E) + rowSums(accu_sep_T))/(rowSums(total_val_vali) + rowSums(total_val_vali))
#PPV (precision)
summ_result[, 9] = rowSums(cross_nor_T)/(rowSums(total_val) - rowSums(cross_nor_E) + rowSums(cross_nor_T))
summ_result[, 10] = rowSums(cross_sep_T)/(rowSums(total_val) - rowSums(cross_sep_E) + rowSums(cross_sep_T))
summ_result[, 11] = rowSums(accu_nor_T)/(rowSums(total_val_vali) - rowSums(accu_nor_E) + rowSums(accu_nor_T))
summ_result[, 12] = rowSums(accu_sep_T)/(rowSums(total_val_vali) - rowSums(accu_sep_E) + rowSums(accu_sep_T))
#TPR (recall)
summ_result[, 13] = rowSums(cross_nor_T)/rowSums(total_val)
summ_result[, 14] = rowSums(cross_sep_T)/rowSums(total_val)
summ_result[, 15] = rowSums(accu_nor_T)/rowSums(total_val_vali)
summ_result[, 16] = rowSums(accu_sep_T)/rowSums(total_val_vali)
#TNR
summ_result[, 17] = rowSums(cross_nor_E)/rowSums(total_val)
summ_result[, 18] = rowSums(cross_sep_E)/rowSums(total_val)
summ_result[, 19] = rowSums(accu_nor_E)/rowSums(total_val_vali)
summ_result[, 20] = rowSums(accu_sep_E)/rowSums(total_val_vali)
#F score
summ_result[, 21] = 2 * summ_result[, 11] * summ_result[, 15] / (summ_result[, 11] + summ_result[, 15])
summ_result[, 22] = 2 * summ_result[, 12] * summ_result[, 16] / (summ_result[, 12] + summ_result[, 16])


#the following code is used to get the information of the feature
feature_detail = feature_detail_mrna[[1]]
for(time in 2:times){
  feature_detail = cbind(feature_detail, feature_detail_mrna[[time]])
}
feature_name = matrix(data = 0, ncol = ncol(feature_detail), nrow = nrow(feature_detail))
feature_ID = matrix(data = 0, ncol = ncol(feature_detail), nrow = nrow(feature_detail))
for(i in 1:ncol(feature_name)){
  feature_name[, i] = z3_E[[5]][feature_detail[, i], 1]
  feature_ID[, i] = z3_E[[5]][feature_detail[, i], 2]
}

write.table(summ_result, file = "summary result.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(feature_name, file = "mrna name.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(feature_ID, file = "mrna ID.txt", quote = FALSE, sep = "\t", row.names = FALSE)



