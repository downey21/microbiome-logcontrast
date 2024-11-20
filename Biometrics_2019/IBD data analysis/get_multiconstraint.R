#set correct directory
getwd()

x_compo <- read.csv(paste(getwd(),"/IBD data analysis/IBD.csv",sep=""), header=T)
x_compo = x_compo[2:98]

ind = (x_compo>0)*1
ind_d = colSums(ind)
pick = ind_d>= 111*0.2
x_or = x_compo[,pick]
p = dim(x_or)[2]

bac_list  = colnames(x_or)
length(bac_list)

bac_genus_list = c()
for (i in 1:length(bac_list)){
  bac_str = bac_list[i]
  bac_l=unlist(strsplit(bac_str,"__"))[2]
  bac_genus= unlist(strsplit(bac_l,"_"))[1]
  bac_genus_list = c(bac_genus_list,bac_genus)
}

constr_list = as.vector(table(bac_genus_list)>1)*1

constr = matrix(0,p+1,sum(constr_list))

for (i in 1:sum(constr_list)){
  id = which(constr_list==1)[i]
  bac_name = names(table(bac_genus_list))[id]
  constr[2:(p+1),i] = (bac_genus_list == bac_name)*1
}


write.table(constr, file = paste(getwd(),"/IBD data analysis/multi_constr.csv",sep=""),row.names=FALSE, col.names=FALSE, sep=",")

