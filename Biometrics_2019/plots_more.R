library("ggplot2", lib.loc="~/R/win-library/3.2")
library("gridExtra", lib.loc="~/R/win-library/3.2")
library("scales", lib.loc="~/R/win-library/3.2")

#
#need to set working directory first
getwd()

#box plot of the coverage probability for each variables
x <- read.csv(paste(getwd(),"/result/coverage_sim50.csv",sep=""), header=FALSE)
p <- dim(x)[1]
N <- dim(x)[2]
data = as.data.frame(matrix(0,p*N,2))
for (i in 1:N){
  data[(p*(i-1)+1):(p*i),1] = x[,i]
  data[(p*(i-1)+1):(p*i),2] = i
}
names(data) = c("P", "cate")

x$cat = rep("others",p)
x$cat[c(2,3,4,6,12,14,17)]="none zero"

#50
plot_true <-ggplot(x,aes(x=x[,1],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,2],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,3],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,4],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")


savepath="CI_n50_p50.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

#100
plot_true <-ggplot(x,aes(x=x[,5],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,6],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,7],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,8],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")

savepath="CI_n100_p50.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

#200
plot_true <-ggplot(x,aes(x=x[,9],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,10],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,11],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,12],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")


savepath="CI_n200_p50.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()
#500
plot_true <-ggplot(x,aes(x=x[,13],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,14],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="ONe")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,15],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,16],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")

savepath="CI_n500_p50.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

# p=100
x <- read.csv(paste(getwd(),"/result/coverage_sim100.csv",sep=""), header=FALSE)
p <- dim(x)[1]
N <- dim(x)[2]
data = as.data.frame(matrix(0,p*N,2))
for (i in 1:N){
  data[(p*(i-1)+1):(p*i),1] = x[,i]
  data[(p*(i-1)+1):(p*i),2] = i
}
names(data) = c("P", "cate")

x$cat = rep("others",p)
x$cat[c(2,3,4,6,12,14,17)]="none zero"

#50
plot_true <-ggplot(x,aes(x=x[,1],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,2],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,3],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,4],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")


savepath="CI_n50_p100.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

#100
plot_true <-ggplot(x,aes(x=x[,5],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,6],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,7],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,8],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")

savepath="CI_n100_p100.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

#200
plot_true <-ggplot(x,aes(x=x[,9],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,10],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="One")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,11],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,12],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")


savepath="CI_n200_p100.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()
#500
plot_true <-ggplot(x,aes(x=x[,13],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="True")
plot_true <- plot_true + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_one <-ggplot(x,aes(x=x[,14],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="ONe")
plot_one <- plot_one + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_no <-ggplot(x,aes(x=x[,15],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="No")
plot_no <- plot_no + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")
plot_wrong <-ggplot(x,aes(x=x[,16],fill = factor(cat)))+geom_histogram(binwidth=0.01,colour="black")+xlim(0.75,1)+labs(title="Wrong")
plot_wrong <- plot_wrong + xlab("Coverage probability") + scale_fill_discrete(name = "Group") + geom_vline(aes(xintercept = 0.95), colour="black", linetype = "dashed")

savepath="CI_n500_p100.pdf"
pdf(savepath,width=10, height = 10)
grid.arrange(plot_true,plot_one,plot_no,plot_wrong,ncol=1, nrow=4)  
dev.off()

