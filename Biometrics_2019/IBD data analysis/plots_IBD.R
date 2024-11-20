library("boot", lib.loc="~/R/win-library/3.2")
library("ggplot2", lib.loc="~/R/win-library/3.2")
library("pastecs", lib.loc="~/R/win-library/3.2")
library("gridExtra", lib.loc="~/R/win-library/3.2")
#library("refund", lib.loc="~/R/win-library/3.2")
library("scales", lib.loc="~/R/win-library/3.2")
library("gplots", lib.loc="~/R/win-library/3.2")
library("ROCR", lib.loc="~/R/win-library/3.2")
library(plotROC)

getwd()

y = c(rep(0,8), rep(1,28))

###IBD variable selection plot
CI = read.csv(paste(getwd(),"/IBD data analysis/results/IBD_result_CI.csv",sep=""), header=FALSE)

p = dim(CI)[1]/2
x = as.data.frame(matrix(0,p*2,4))
x[1:p,1] = seq(1,p);
x[(p+1):(2*p),1] = seq(1,p);
x[1:(p*2),2] = CI[1:(p*2),1];
x[1:(p*2),3] = CI[1:(p*2),2];
x[1:p,4] = "De-biased"
x[(p+1):(p*2),4] = "Lasso"
names(x)=c("bacteria","Parameter","se", "Group")
x[1,3]=0;

de = x[1:p,];

plot <- ggplot(x, aes(x = bacteria, y =Parameter, shape = Group)) + geom_point() + scale_shape_manual(values = c(1,4))
plot<- plot + geom_errorbar(data = de, aes(ymin = Parameter - se, ymax = Parameter + se), width = 0.2)  
plot <- plot + xlab("Bacteria")+ ylim(-0.4,0.3)
plot <-plot+theme(axis.title.x=element_text(size=8))
plot <-  plot +annotate("text", x=2, y=-0.25, label="copri",color="red",size=4)
plot <-plot+annotate("text", x=12, y=-0.32, label="bromii",color="red",size=4)
plot <-plot+annotate("text", x=34, y=-0.3, label="leptum",color="red",size=4)
plot <-plot+annotate("text", x=64, y=0.23, label="coli",color="red",size=4)
plot <-plot+annotate("text", x=70, y=0.28, label="gnavus",color="red",size=4)
plot



ggsave(paste(getwd(),"/IBD data analysis/plots/IBD_vs_new.pdf",sep=""), width = 20, height = 10, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/IBD_vs_new.png",sep=""), width = 20, height = 10, units = "cm")


####
##score plot
score <- read.csv(paste(getwd(),"/IBD data analysis/results/score_prob.csv",sep=""), header=FALSE)

names(score) <- c("score", "y")
score$prob <- exp(score$score)/(1+exp(score$score))

p1 <-ggplot(score) + geom_line(aes(x=score,y=prob))
#p1<-ggplot(score, aes(x = score, y=y,color =factor(y)))+geom_point()
p1 <- p1 + geom_point(aes(x = score, y=y,color =factor(y)))
p1 <- p1+ labs(color="Y",size=30) +scale_color_discrete(labels=c("Control", "Case"))
p1
ggsave(paste(getwd(),"/IBD data analysis/plots/score_plot.pdf",sep=""), width = 20, height = 20, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/score_plot.png",sep=""), width = 20, height = 10, units = "cm")


#score plot using only five
score1 <- read.csv(paste(getwd(),"/IBD data analysis/results/score_prob_5.csv",sep=""), header=FALSE)
names(score1) <- c("score", "y")
m1 <- glm(score1$y~score1$score,family=binomial)

score1$prob <- m1$fitted.values


p1 <-ggplot(score1) + geom_line(aes(x=score,y=prob))
p1 <- p1 + geom_point(aes(x = score, y=y,color =factor(y)))
p1 <- p1+ labs(color="Y",size=30) +scale_color_discrete(labels=c("Control", "Case"))
p1
ggsave(paste(getwd(),"/IBD data analysis/plots/score_plot_5_nointer.pdf",sep=""), width = 20, height = 20, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/score_plot_5_nointer.png",sep=""), width = 20, height = 10, units = "cm")


# individual plot
porp<- read.csv(paste(getwd(),"/IBD data analysis/results/x_select.csv",sep=""), header=FALSE)
names(porp) <-c("v1","v2","v3","v4","v5")

#all log plot
porp_resize3 <- matrix(0,111*5,3)
porp_resize3[,1] <- cbind(porp$v1,porp$v2,porp$v3,porp$v4,porp$v5)
porp_resize3[,2] <- c(rep(1,111),rep(2,111),rep(3,111),rep(4,111),rep(5,111))
porp_resize3[,3] <- c(rep(0,26),rep(1,85),rep(0,26),rep(1,85),rep(0,26),rep(1,85),rep(0,26),rep(1,85),rep(0,26),rep(1,85))

porp_resize3 = as.data.frame(porp_resize3)
names(porp_resize3) = c("abundance", "group", "treatment")

p2 <- ggplot(porp_resize3)+geom_boxplot(aes(x = factor(group), y= log(abundance), fill = factor(treatment)), width=0.3)
p2 <- p2 + labs(fill=NULL,size=30) +scale_fill_discrete(labels=c("Control", "Case"))
p2 <- p2 + xlab("group")+ylab("log relative abundance")+ theme(legend.title = element_text(face = "bold"))
p2 <- p2 + scale_x_discrete(labels=c("Prevotella_copri","Ruminococcus_bromii","Clostridium_leptum","Escherichia_coli","Ruminococcus_gnavus"))
p2

ggsave(paste(getwd(),"/IBD data analysis/plots/abund_log_all.png",sep=""), width = 30, height = 15, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/abund_log_all.pdf",sep=""), width = 30, height = 15, units = "cm")


# individual testing:

logprop <-read.csv(paste(getwd(),"/IBD data analysis/plots/log_proportion.csv"), header=FALSE)

pvalue = rep(0,77)
for (i in 1:77){
  case = logprop[1:26,i]
  control = logprop[27:111,i]
  test = wilcox.test(case, control)
  pvalue[i] =test$p.value
}
padjust = p.adjust(pvalue)

pv = as.data.frame(matrix(0,77,3))
pv$V1 = pvalue
pv$V2 = seq(1,77)
pv$V3[c(1,11,33,63,69)]=1
names(pv) = c("pvalue", "taxa","group")

p<-ggplot(pv)+geom_point(aes(x= taxa, y = -log10(pvalue), color= factor(group)))
p<-p + scale_color_manual(values=c("black","red"),labels=c("others", "Selected variables"))
p <- p +annotate("text", x=1, y=-log10(pvalue[1])+0.3, label="copri")
p <-p+annotate("text", x=11, y=-log10(pvalue[11])-0.3, label="bromii")
p <-p+annotate("text", x=33, y=-log10(pvalue[33])+0.3, label="leptum")
p <-p+annotate("text", x=63, y=-log10(pvalue[63])+0.3, label="coli")
p <-p+annotate("text", x=69, y=-log10(pvalue[69])+0.5, label="gnavus")
p <- p + xlab("taxa number") + labs(color=NULL)
p

ggsave(paste(getwd(),"/IBD data analysis/plots/pvalue.png"), width = 15, height = 15, units = "cm")

#adjust
pv_adj= as.data.frame(matrix(0,77,3))
pv_adj$V1 = padjust
pv_adj$V2 = seq(1,77)
pv_adj$V3[c(1,11,33,63,69)]=1
names(pv_adj) = c("pvalue", "taxa","group")

p<-ggplot(pv_adj)+geom_point(aes(x= taxa, y = -log10(pvalue), color= factor(group)))
p<-p + scale_color_manual(values=c("black","red"),labels=c("others", "Selected variables"))
p <- p +annotate("text", x=1, y=-log10(padjust[1])+0.3, label="copri")
p <-p+annotate("text", x=11, y=-log10(padjust[11])-0.3, label="bromii")
p <-p+annotate("text", x=33, y=-log10(padjust[33])+0.3, label="leptum")
p <-p+annotate("text", x=63, y=-log10(padjust[63])+0.3, label="coli")
p <-p+annotate("text", x=69, y=-log10(padjust[69])+0.5, label="gnavus")
p <- p + xlab("taxa number") + labs(color=NULL)
p
ggsave(paste(getwd(),"/IBD data analysis/plots/pvalue_adj.png"), width = 15, height = 15, units = "cm")


#stability plot
stability <- read.csv(paste(getwd(),"/IBD data analysis/plots/Stability.csv"), header=FALSE)

stab <- matrix(0,21*77,4)
stab[,1] <- rep(seq(0,0.5,0.025),77)
for (i in 1:71){
  stab[((i-1)*21+1):(i*21),2] = 1-stability[,i]
  if (i==1){
  stab[((i-1)*21+1):(i*21),4] = 1
  }else if(i==11){
    stab[((i-1)*21+1):(i*21),4] = 2
  }else if(i==33){
    stab[((i-1)*21+1):(i*21),4] = 3
  }else if(i==63){
    stab[((i-1)*21+1):(i*21),4] = 4
  }else if(i==69){
    stab[((i-1)*21+1):(i*21),4] = 5
  }
}  
stab[,3] <- rep(seq(1,77),each=21)

stab <- as.data.frame(stab)
names(stab) <- c("x","prob","group","ind")

p3 <-ggplot(stab, aes(x=x, y = prob,group=factor(group),color = factor(ind)))+geom_line(size=1)+scale_x_reverse()
p3 <-p3 + scale_color_manual(values=c("grey","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),labels=c("others", "copri","bromii","leptum","coli","gnavus")) + xlab("lambda")+ylab("selecton probability")+labs(color=NULL,size=30)
p3

ggsave(paste(getwd(),"/IBD data analysis/plots/stability.png"), width = 15, height = 15, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/stability.pdf"), width = 15, height = 15, units = "cm")



#stability plot for not using constraints
stability <- read.csv(paste(getwd(),"/IBD data analysis/plots/stability_noc.csv"), header=FALSE)
stab <- matrix(0,21*77,4)
stab[,1] <- rep(seq(0,0.5,0.025),77)
for (i in 1:71){
  stab[((i-1)*21+1):(i*21),2] = 1-stability[,i]
  if (i==1){
    stab[((i-1)*21+1):(i*21),4] = 1
  }else if(i==11){
    stab[((i-1)*21+1):(i*21),4] = 2
  }else if(i==33){
    stab[((i-1)*21+1):(i*21),4] = 3
  }else if(i==63){
    stab[((i-1)*21+1):(i*21),4] = 4
  }else if(i==69){
    stab[((i-1)*21+1):(i*21),4] = 5
  }
}  
stab[,3] <- rep(seq(1,77),each=21)

stab <- as.data.frame(stab)
names(stab) <- c("x","prob","group","ind")

p3 <-ggplot(stab, aes(x=x, y = prob,group=factor(group),color = factor(ind)))+geom_line(size=1)+scale_x_reverse()
p3 <-p3 + scale_color_manual(values=c("grey","#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"),labels=c("others", "copri","bromii","leptum","coli","gnavus")) + xlab("lambda")+ylab("selecton probability")+labs(color=NULL,size=30)
p3

ggsave(paste(getwd(),"/IBD data analysis/plots/stability_noc.png"), width = 15, height = 15, units = "cm")
ggsave(paste(getwd(),"/IBD data analysis/plots/stability_noc.pdf"), width = 15, height = 15, units = "cm")


## ROC
## all variable
prob <- read.csv(paste(getwd(),"/IBD data analysis/ROC_result.csv",sep=""), header=FALSE)
y = c(rep(0,8), rep(1,28))
nrep = dim(prob)[1]

AUC = matrix(0, nrep,4)
for (i in 1: nrep){
  pre<-prediction(t(prob[i,1:36]), y)
  auc1 <-performance(pre, measure = "auc")
  #plot(roc)
  pre2<-prediction(t(prob[i,37:72]),y)
  auc2<-performance(pre2, measure = "auc")
  #plot(roc2)
  pre3<-prediction(t(prob[i,73:108]),y)
  auc3<-performance(pre3, measure = "auc")
  #
  pre4<-prediction(t(prob[i,109:144]),y)
  auc4 <-performance(pre4, measure = "auc")
  AUC[i,] = as.vector(c(auc1@y.values[[1]],auc2@y.values[[1]],auc3@y.values[[1]],auc4@y.values[[1]]))
}

colMeans(AUC)
apply(AUC, 2, sd)

## 5 variable
prob <- read.csv(paste(getwd(),"/IBD data analysis/ROC_result_select.csv",sep=""), header=FALSE)

prob1 = matrix(0,50,144)
for (i in 1:50){
  prob1[i,] = prob[(144*(i-1)+1):(144*i),1]
}
prob = as.data.frame(prob1)

y = c(rep(0,8), rep(1,28))
nrep = dim(prob)[1]


AUC = matrix(0, nrep,4)
for (i in 1: nrep){
  pre<-prediction(t(prob[i,1:36]), y)
  auc1 <-performance(pre, measure = "auc")
  #plot(roc)
  pre2<-prediction(t(prob[i,37:72]),y)
  auc2<-performance(pre2, measure = "auc")
  #plot(roc2)
  pre3<-prediction(t(prob[i,73:108]),y)
  auc3<-performance(pre3, measure = "auc")
  #
  pre4<-prediction(t(prob[i,109:144]),y)
  auc4 <-performance(pre4, measure = "auc")
  AUC[i,] = as.vector(c(auc1@y.values[[1]],auc2@y.values[[1]],auc3@y.values[[1]],auc4@y.values[[1]]))
}

colMeans(AUC)
apply(AUC, 2, sd)

#ROC for random forest
prob <- read.csv(paste(getwd(),"/IBD data analysis/ROC_result_rf.csv",sep=""), header=FALSE)
y = c(rep(0,8), rep(1,28))
nrep = dim(prob)[1]

AUC = matrix(0, nrep,1)

for (i in 1: nrep){
  pre<-prediction(t(prob[i,]), y)
  auc1 <-performance(pre, measure = "auc")
  #plot(roc)
  AUC[i,] = as.vector(auc1@y.values[[1]])
}

mean(AUC)
sd(AUC)

