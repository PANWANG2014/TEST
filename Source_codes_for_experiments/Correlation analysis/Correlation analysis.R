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
#MMI_location is used to match the validation set with the expression to get the location
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

#mirtigo_corr is used to get the top 1 percent pairs of mirtigo and the intersection pairs of mirtigo
mirtigo_corr = function(mirna_use1, mrna_use1, wMRE, sumup){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_agg = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)

  for(m in 1:ncol(mirna_use)){
    pro_temp = (mirna_use[, m] %*% t(mrna_use[, m])) * wMRE / sumup[m]
    pro_agg = pro_agg + pro_temp
    num = floor(length(which(pro_temp != 0))/100)
    temp = order(pro_temp, decreasing = TRUE)[1 : num]
  }
  thres_num = floor(length(which(pro_agg != 0))/100)
  sort_agg = order(pro_agg, decreasing = TRUE)[1 : thres_num]
  return(sort_agg)
}
#promise_corr is used to get the top 1 percent pairs of promise conserved mrna
promise_corr = function(mirna_use1, mrna_use1, wMRE){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  pro_agg_mrna = matrix(data = 0, ncol = nrow(mrna_use), nrow = nrow(mirna_use), byrow = TRUE)
  
  wMRE_new = t(wMRE)
  for(m in 1:ncol(mirna_use)){
    x = matrix(mrna_use[, m])
    rownames(x) = rownames(wMRE_new)
    z = matrix(mirna_use[, m])
    rownames(z) = colnames(wMRE_new)
    rs = roleswitch(x, z, wMRE_new)
    pro_temp = t(rs$p.x)
    pro_agg_mrna = pro_agg_mrna + pro_temp
  }
  thres_num = floor(length(which(pro_agg_mrna != 0)) / 100)
  sort_agg_mrna = order(pro_agg_mrna, decreasing = TRUE)[1 : thres_num]
  return(sort_agg_mrna)
}
#targetscan_corr is used to get the top 1 percent pairs of targetscan wMRE
targetscan_corr = function(wMRE){
  thres_num = floor(length(which(wMRE != 0))/100)
  target_rank = order(wMRE, decreasing = TRUE)[1 : thres_num]
  return(target_rank)
}
#correlation_detail is used to get the detailed information of the correlation result
correlation_detail = function(pos_vec, mirna_use1, mrna_use1, mirna_name, mrna_name){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  cor_matrix = matrix(data = 0, ncol = 5, nrow = length(pos_vec), byrow = TRUE)
  colnames(cor_matrix) = c("rank", "mirna name", "mrna name", "mrna ID", "Pearson correlation")
  
  mrna_vec = ceiling(pos_vec/length(mirna_name))
  mirna_vec = pos_vec - (mrna_vec - 1) * length(mirna_name)
  for(i in 1:length(pos_vec)){
    cor_matrix[i, 1] = i
    cor_matrix[i, 2] = mirna_name[mirna_vec[i]]
    cor_matrix[i, 3] = mrna_name[mrna_vec[i], 1]
    cor_matrix[i, 4] = mrna_name[mrna_vec[i], 2]
    cor_matrix[i, 5] = cor.test(mrna_use[mrna_vec[i], ], mirna_use[mirna_vec[i], ], method = "pearson")$estimate
  }
  cor_med = median(as.numeric(cor_matrix[, 5]))
  return(list(cor_matrix, cor_med))
}
#base_correlation is used to calculate the baseline correlation distribution and median
base_correlation = function(mirna_use1, mrna_use1){
  mirna_use = matrix(as.numeric(mirna_use1), nrow = nrow(mirna_use1))
  mrna_use = matrix(as.numeric(mrna_use1), nrow = nrow(mrna_use1))
  corMat = cor(t(mirna_use), t(mrna_use))
  med = median(corMat)
  return(list(corMat, med))
}

####################################################################################################

### Part II. INPUT DATA: INPUT FOR POPULATION-LEVEL EVALUATION ###

#the following package is used for promise method
library(Roleswitch)
#first read in the files we need to use, which contains:
#input for mirtigo algorithm: mirna expression, mrna expression, CWCS_matrix, mirna list, mrna list
#input for promise algorithm: qMRE_matrix
mrna = as.matrix(read.table("mrna_list.txt", head = TRUE, sep = "\t"))
mirna = as.matrix(read.table("mirna_list.txt", head = TRUE, sep = "\t"))
qMRE_conserved_matrix = as.matrix(read.table("conserved_qMRE.txt", head = TRUE, sep = "\t"))
CWCS_matrix1 = as.matrix(read.table("wMRE_all.txt", head=TRUE, sep = "\t"))
CWCS_matrix = abs(CWCS_matrix1)
name_cancer = as.matrix(read.table("Sample-to-sample file.txt", head = TRUE, sep = "\t"))
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

z8 = sum_mirtigo(z3[[1]], z3[[2]], z7[[3]])
z9 = mirtigo_corr(z3[[1]], z3[[2]], z7[[3]], z8) #get the result of mirtigo
z10 = promise_corr(z3[[1]], z3[[2]], z6[[3]]) #get the result of promise conserved qMRE
z11 = targetscan_corr(z7_tar[[3]])

z13 = correlation_detail(z9, z3[[1]], z3[[2]], z3[[3]], z3[[5]]) #this is for mirtigo
z14 = correlation_detail(z10, z3[[1]], z3[[2]], z3[[3]], z3[[5]]) #this is for promise conversed mrna
z15 = correlation_detail(z11, z4[[1]], z4[[2]], z4[[3]], z4[[5]]) #this is for targetscan
base_cor = base_correlation(z3[[1]], z3[[2]])

median_cor = matrix(data = 0, ncol = 2, nrow = 4, byrow = TRUE)
colnames(median_cor) = c("methods", "median")
median_cor[, 2] = c(z13[[2]], z14[[2]], z15[[2]], base_cor[[2]])
median_cor[, 1] = c("mirtigo", "promise conserved", "targetscan CWCS", "base line")

write.table(z13[[1]],file="correlation detail inforamtion of mirtigo.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z14[[1]],file="correlation detail inforamtion of promise conserved mrna.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(z15[[1]],file="correlation detail inforamtion of targetscan.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(median_cor,file="correlation median.txt", quote = FALSE, sep = "\t", row.names = FALSE)

a1 = cbind(as.numeric(z13[[1]][,5]), rep(1, length(z13[[1]][,5])))
a2 = cbind(as.numeric(z14[[1]][,5]), rep(1, length(z14[[1]][,5])))
a3 = cbind(as.numeric(z15[[1]][,5]), rep(1, length(z15[[1]][,5])))
a00 = as.vector(base_cor[[1]])
a0 = cbind(as.numeric(a00), rep(0, length(a00)))
b1 = rbind(a1, a0)
b2 = rbind(a2, a0)
b3 = rbind(a3, a0)
c1 = wilcox.test(b1[, 1] ~ b1[, 2], data = b1)$p.value
c2 = wilcox.test(b2[, 1] ~ b2[, 2], data = b2)$p.value
c3 = wilcox.test(b3[, 1] ~ b3[, 2], data = b3)$p.value

pdf(file = "Rplot.pdf", width = 10, height = 7.5)
par(mar=c(5, 5, 2, 2) + 0.1)
plot(density(as.numeric(z14[[1]][,5])), col=rgb(r = 56, g = 108, b = 176, maxColorValue = 255), xlim=c(-0.5, 0.5), ylim=c(0,6), lwd=2, xlab="", ylab="", main="", cex.axis=1.5) #this is for promise
lines(density(as.numeric(z15[[1]][,5])), col=rgb(r = 217, g = 95, b = 2, maxColorValue = 255),lwd=2) #this is for targetscan
lines(density(unlist(base_cor)), col="black", lty=2, lwd=2)#this is for base
lines(density(as.numeric(z13[[1]][,5])), col=rgb(r = 240, g = 2, b = 127, maxColorValue = 255), lwd=2) #this is for mirtigo
legend("topleft",legend = c(paste("miRTIGO (n=",nrow(z13[[1]]),", shift=",round(z13[[2]] - base_cor[[2]], 5),", p=", signif(c1, 3),")", sep = ""),
                            paste("ProMISe (n=",nrow(z14[[1]]),", shift=",round(z14[[2]] - base_cor[[2]], 5),", p=",signif(c2, 3),")", sep = ""), 
                            paste("TargetScan (n=",nrow(z15[[1]]),", shift=",round(z15[[2]] - base_cor[[2]], 5),", p=",signif(c3, 3),")", sep = ""), 
                            "all the miRNA-mRNA correlations"),
       col = c(rgb(r=240, g=2, b=127, maxColorValue = 255), 
               rgb(r=56, g=108, b=176, maxColorValue = 255),
               rgb(r=217, g=95, b=2, maxColorValue = 255), 
               "black"), ncol=1, lwd = c(3,3,3,3), lty=c(1,1,1,2), cex=1.5, bty="n")
mtext("Density", side = 2, line = 3, cex = 2)
mtext("Pearson correlation coefficient", side = 1, line = 3.5, cex = 2)
dev.off()
