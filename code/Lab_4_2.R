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
setwd("D:/Courses/UC Berkeley/Fall 2015/Statistical Models Theory and Application/Lab/Lab4/lab")
load('embryo_data.Rda')
load('embryo_imgs.Rda')
source('./utils.R')

############AT this stage, read in finaldata.csv which contains all the infomration and fiji features

data<-read.csv("finaldata.csv",header=T) ##train+test+fiji features
data$label_pred<-0 ##build a fake label for later testing
id <- sum(num.sp[1:150]) ##after load embryo_data.Rda, you will have this num.sp
train<-data[1:id,]
test<-data[(id+1):dim(data[1])[1],]

############Test spatial autocorrelation

JC_1<-rep(0,150)
JC_2<-rep(0,150)
JC_3<-rep(0,150)
###compute Joint counts statistics for each embryo
for (j in 1:150){
  ##obtain x,y, and label for each embryo
  if (j<150){
    emb_j<-data[(sum(num.sp[1:j])+1):(sum(num.sp[1:(j+1)])),]
    obj<-matrix(cbind(emb_j$x,emb_j$y),ncol=2)
    nb<-dnearneigh(obj,20,200)
    jointest<-joincount.test(as.factor(emb_j$label),nb2listw(nb))
    JC_1[j]<-jointest[[1]]$p.value ##class 1
    JC_2[j]<-jointest[[2]]$p.value ##class 2
    JC_3[j]<-jointest[[3]]$p.value ##class 3
  }
  if (j==150){
    emb_j<-data[(sum(num.sp[1:149])+1):(sum(num.sp[1:150])),]
    obj<-matrix(cbind(emb_j$x,emb_j$y),ncol=2)
    nb<-dnearneigh(obj,20,200)
    jointest<-joincount.test(as.factor(emb_j$label),nb2listw(nb))
    JC_1[j]<-jointest[[1]]$p.value ##class 1
    JC_2[j]<-jointest[[2]]$p.value ##class 2
    JC_3[j]<-jointest[[3]]$p.value ##class 3
  }
}
JC<-data.frame(cbind(JC_1,JC_2,JC_3))
JC$id<-seq(1:150)

###plot p-values for three classes
ggplot(JC, aes(x=id, y = value, color = variable)) + 
  geom_point(aes(y = JC_1, col = "y1")) + 
  geom_point(aes(y = JC_2, col = "y2")) +
  geom_point(aes(y = JC_3, col = "y3"))
############Tune individual models using cross validation

###method 1: random forest
##check python code
##random forest results suggest spatial features are the most important ones,
##this makes us to focus on spatial dependnce, for instance, k-nearest neighbor

##################Method 2: Tuning k-nearest neighbor, solely based on spatial information:coordinates and distance
###use CV to find k

set.seed(100)
klist <- seq(1,200,5) # we test values of k from 1 to 200
nfolds <- 6 # we make 6-fold cross-validation, since 6 can be divided by 32664 train without reminder

##prepare training set
x.train=train[c("x","y","distance")]  ###one can add more features, but since adding fiji features could make knn worse, we only include three spatial features.
y.train=train[4]

# Perform nfolds-cross validation of kNN, for the values of k in klist
# Number of instances
n.train <- nrow(x.train)

# Matrix to store predictions
p.cv <- matrix(NA, nfolds, length(klist))
p.cv_class1 <- matrix(NA, nfolds, length(klist))
p.cv_class2 <- matrix(NA, nfolds, length(klist))
p.cv_class3 <- matrix(NA, nfolds, length(klist))
# Prepare the folds
s <- split(sample(n.train),rep(1:nfolds,length=n.train))

# Cross-validation for KNN based on spatial information
for (i in seq(nfolds)){
  for (j in (1:length(klist))) {
    k<-klist[j]
    knn_fit <- knn(train=x.train[-s[[i]],], test=x.train[s[[i]],],cl=y.train[-s[[i]],],k=k)
    conf<-confusionMatrix(y.train[s[[i]],],knn_fit)
    p.cv[i,j]<-conf$overall[1]
    p.cv_class1[i,j]<-conf$byClass[1] # sensitivity 
    p.cv_class2[i,j]<-conf$byClass[2]
    p.cv_class3[i,j]<-conf$byClass[3]    
  }
}

CV_error<-colMeans(p.cv) 
CV_error_class1<-colMeans(p.cv_class1)
CV_error_class2<-colMeans(p.cv_class2)
CV_error_class3<-colMeans(p.cv_class3)
knn_CV<-data.frame(cbind(CV_error,CV_error_class1,CV_error_class2,CV_error_class3))
knn_CV$k<-seq(1,200,5)

##best k = 151 with overall prediction accuracy 86.8%

##this graph shows out prediction on class 1 is better than class 2 and 3
ggplot(knn_CV, aes(x=k, y = value, color = variable)) + 
  geom_line(aes(y = CV_error, col = "overall")) + 
  geom_line(aes(y = CV_error_class1, col = "class 1")) +
  geom_line(aes(y = CV_error_class2, col = "class 2")) +
  geom_line(aes(y = CV_error_class3, col = "class 3")) + ylab("Prediction accuracy rate")


###########################Method 3: Tuning weighted k-nearest neighbor
###use CV to find k
set.seed(100)
# Matrix to store predictions
p.cv.G <- matrix(NA, nfolds, length(klist))
p.cv.O <- matrix(NA, nfolds, length(klist))
p.cv.E <- matrix(NA, nfolds, length(klist))

wx.train=train[c("x","y","distance","label")]
# Cross-validation for weighted KNN
for (i in seq(nfolds)){
  for (j in (1:length(klist))) {
    k<-klist[j]
    wknn_fit_G <- kknn(as.factor(label)~., wx.train[-s[[i]],], wx.train[s[[i]],], na.action = na.omit(),
                     k = k, distance = 2, kernel = "gaussian", ykernel = NULL, scale=TRUE)
    wknn_fit_O <- kknn(as.factor(label)~., wx.train[-s[[i]],], wx.train[s[[i]],], na.action = na.omit(),
                       k = k, distance = 2, kernel = "optimal", ykernel = NULL, scale=TRUE)
    wknn_fit_E <- kknn(as.factor(label)~., wx.train[-s[[i]],], wx.train[s[[i]],], na.action = na.omit(),
                       k = k, distance = 2, kernel = "epanechnikov", ykernel = NULL, scale=TRUE)   
    p.cv.G[i,j]<-confusionMatrix(y.train[s[[i]],],wknn_fit_G$fitted.values)$overall[1]
    p.cv.O[i,j]<-confusionMatrix(y.train[s[[i]],],wknn_fit_O$fitted.values)$overall[1]
    p.cv.E[i,j]<-confusionMatrix(y.train[s[[i]],],wknn_fit_E$fitted.values)$overall[1]
  }
}

CV_error_G<-colMeans(p.cv.G) 
CV_error_O<-colMeans(p.cv.O)
CV_error_E<-colMeans(p.cv.E)
wknn_CV<-data.frame(cbind(CV_error_G,CV_error_O,CV_error_E))
wknn_CV$k<-seq(1,200,5)

##best k= 151 with overall prediction accuracy 86.7%

ggplot(wknn_CV, aes(x=k, y = value, color = variable)) + 
  geom_line(aes(y = CV_error_G, col = "Gaussain kernel")) + 
  geom_line(aes(y = CV_error_O, col = "Optimal kernel")) +
  geom_line(aes(y = CV_error_E, col = "Epanechnikov")) + ylab("Prediction accuracy rate")

####
fitted<-list()
for (i in seq(nfolds)){
  probab<- kknn(as.factor(label)~., wx.train[-s[[i]],], wx.train[s[[i]],], na.action = na.omit(),
                              k = 151, distance = 2, kernel = "epanechnikov", ykernel = NULL, scale=TRUE)$prob
  fitted[[i]]<-cbind(probab,wx.train[s[[i]],]$label)
}

output_wknn<-data.frame(do.call("rbind", fitted))
names(output_wknn)[names(output_wknn)=="X1"] <- "prob.class1"
names(output_wknn)[names(output_wknn)=="X2"] <- "prob.class2"
names(output_wknn)[names(output_wknn)=="X3"] <- "prob.class3"
names(output_wknn)[names(output_wknn)=="V4"] <- "label"
write.csv(output_wknn,"output_wknn.csv") ##to be called by python to construct ROC curve

#############Method 4: QDA

dist<-dist(matrix(cbind(train$x,train$y),ncol=2))
mda_pred <- qda(as.factor(label)~x+y+distance,data=train)
posteriors<-predict(mda_pred,test,type="post")
confusionMatrix(softmax(posteriors), test$label)
confusionMatrix(posteriors, test$label)


#############################################visulize the predicted labels for test set

x.train<-train[1:3] 
y.train<-train[,4]
x.test<-test[1:3]
test$label_pred<-knn(x.train, x.test,cl=y.train,k=151) ##classify on test set based on knn

confusionMatrix(test$label,test$label_pred)

data2<-rbind(train,test)
# compare super pixel feature values by organ
###jj ranges from 151 to 170, one can simply pick one
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
# image(flipImg(imgs[,,jj]), col=img.colors)
# image(flipImg(embryo.labels.raw[,,jj]), col=seg.colors)
image(flipImg(embryo.labels[,,jj]), col=seg.colors)
image(flipImg(embryo.seg_new2), col=seg.colors)

