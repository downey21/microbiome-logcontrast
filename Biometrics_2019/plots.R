library("boot", lib.loc="~/R/win-library/3.2")
library("Rcpp")
library("ggplot2", lib.loc="~/R/win-library/3.2")
library("pastecs", lib.loc="~/R/win-library/3.2")
library("gridExtra", lib.loc="~/R/win-library/3.2")
library("scales", lib.loc="~/R/win-library/3.2")

#this part of the code is used to generate the plots 

#need to set working directory first
getwd()

#coverage probability for simulation based on generated data
#need to change the filename from 50 to 100 manually. 
#boxplots for coverage probability
x <- read.csv(paste(getwd(),"/result/coverage_sim50.csv",sep=""), header=FALSE)
p <- dim(x)[1]
N <- dim(x)[2]
data = as.data.frame(matrix(0,p*N,2))
for (i in 1:N){
  data[(p*(i-1)+1):(p*i),1] = x[,i]
  data[(p*(i-1)+1):(p*i),2] = i
}
names(data) = c("P", "cate")

#coverage probability plot
plot <- ggplot(data, aes(x = factor(cate), y = P)) +geom_boxplot(position = position_dodge(width = 0), width = 0.5)+stat_summary(fun.y = "mean", geom = "point", size = 2, color = "red")
#set y range
#plot <- plot + ylim(0,1)
#add dashed line
plot <- plot + geom_hline(yintercept = 0.95, linetype = "dashed",size =0.5) + geom_vline(xintercept = 4.5, linetype = "dashed",size =0.5) + geom_vline(xintercept = 8.5, linetype = "dashed",size =0.5) + geom_vline(xintercept = 12.5, linetype = "dashed",size =0.5) 
#add new labels on x axis
plot <- plot + scale_x_discrete(labels = c("True", "One", "No", "Wrong","True", "One", "No", "Wrong","True", "One", "No", "Wrong","True", "One", "No", "Wrong"))
#hide title on x axis
plot <- plot +ylab("Coverage Probability")+xlab("Constraints")
#+xlab("n =50                                      n=100                                      n=200                                      n=500")
#add text
#plot <- plot + annotate("text", x=2.5, y =0, label = "n =50")
# theme
plot <- plot + theme_bw()
#hide lines
plot <- plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#rotate 90 
plot <- plot + theme(axis.text.x = element_text(angle =90, hjust =1, vjust = 0.5))
plot <- plot+coord_cartesian(ylim = c(0.5,1))

plot


## Length of CI
x <- read.csv(paste(getwd(),"/result/length_sim50.csv",sep=""), header=FALSE)
p <- dim(x)[1]
N <- dim(x)[2]
data = as.data.frame(matrix(0,p*N,2))
for (i in 1:N){
  data[(p*(i-1)+1):(p*i),1] = x[,i]
  data[(p*(i-1)+1):(p*i),2] = i
}
names(data) = c("P", "cate")

#length of CI plot
plot <- ggplot(data, aes(x = factor(cate), y = P)) +geom_boxplot(position = position_dodge(width = 0), width = 0.5)+stat_summary(fun.y = "mean", geom = "point", size = 2, color = "red")
#set y range
#plot <- plot + ylim(0,1)
#add dashed line
plot <- plot +  geom_vline(xintercept = 4.5, linetype = "dashed",size =0.5) + geom_vline(xintercept = 8.5, linetype = "dashed",size =0.5) + geom_vline(xintercept = 12.5, linetype = "dashed",size =0.5) 
#add new labels on x axis
plot <- plot + scale_x_discrete(labels = c("True", "One", "No", "Wrong","True", "One", "No", "Wrong","True", "One", "No", "Wrong","True", "One", "No", "Wrong"))
#hide title on x axis
plot <- plot +ylab("Expected Length of CI")+xlab("Constraints")
#+xlab("n =50                                      n=100                                      n=200                                      n=500")
#add text
#plot <- plot + annotate("text", x=2.5, y =0, label = "n =50")
# theme
plot <- plot + theme_bw()
#hide lines
plot <- plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#rotate 90 
plot <- plot + theme(axis.text.x = element_text(angle =90, hjust =1, vjust = 0.5))
plot <- plot+coord_cartesian(ylim = c(0,4))
plot
