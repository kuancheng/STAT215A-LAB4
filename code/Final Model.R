library(RColorBrewer)
library(grDevices)
library(R.matlab)
library(ggplot2)
library(plyr)
library(class) ##pkg for knn
library(pROC) ##multiclass roc
library(kknn) ##weighted knn
library(e1071)
library(caret) ##confusion matrix
library(lattice)
library(spdep)
load('embryo_data.Rda')
load('embryo_imgs.Rda')
source('./utils.R')

############AT this stage, read in Final_data.csv which contains all the infomration and fiji features

data<-read.csv("Final_data.csv",header=T) ##train+test+fiji features
data$label_pred<-0 ##build a fake label for later testing
id <- sum(num.sp[1:150]) ##after load embryo_data.Rda, you will have this num.sp
train<-data[1:id,]
test<-data[(id+1):dim(data[1])[1],]

############we will test the performance of KNN on Test set using spatial features, k=151 is selected as best parameter
x.train<-train[1:3] ##feel free to add more fiji features, but it will only make classification worse
y.train<-train[,4]
x.test<-test[1:3]
test$label_pred<-knn(x.train, x.test,cl=y.train,k=151) ##classify on test set based on knn


##########test the performance of this classifier using KNN
confusionMatrix(test$label,test$label_pred)


################visualize predicted labels 
data2<-rbind(train,test)

###jj ranges from 151 to 170, one can simply pick one to plot and compare predicted image with sp-labeled image
jj<-159
emb<-data2[(sum(num.sp[1:(jj-1)])+1):(sum(num.sp[1:jj])),]
emb$sp_id<-sort(unique(as.numeric(embryo.seg[,,jj])))[-1] ##this will be sp id
embryo.seg_new<-embryo.seg[,,jj] 

###match predicted labels
for (s in (1:length(emb$sp_id))){
  embryo.seg_new<-ifelse(embryo.seg_new==emb$sp_id[s],emb$label_pred[s],embryo.seg_new)
}

embryo.seg_new2<-matrix(as.numeric(embryo.seg_new),ncol=600)
seg.colors <- c(brewer.pal(12, 'Set3'))

par(mfrow=c(2, 1))
image(flipImg(embryo.labels[,,jj]), col=seg.colors)
image(flipImg(embryo.seg_new2), col=seg.colors)
