library(dplyr)
library(RColorBrewer)
library(gridExtra)
library(devtools)
library(R.matlab)
library(rgl)
library(ggbiplot)
library(moments)
library(reshape2)
load('embryo_data.Rda')
load('embryo_imgs.Rda')
source('utils.R')

fiji.directory <- '../lab/raw_fiji/'


# set up colors to be plotted
img.colors <- sapply(1:100, function(i) paste0('gray', i)) 
seg.colors <- c(brewer.pal(12, 'Set3'))
fiji.colors <- c(rev(brewer.pal(9, 'Blues')), brewer.pal(9, 'YlOrRd'))

compareToRaw <- function(number){
  # a function that plots embryos with raw labels and with experts labels for comparison
  # Args:
  #  number: which number of image are going to be compared
  par(mfrow = c(1, 2))
  image(flipImg(embryo.labels[, , number]), col = seg.colors)
  image(flipImg(embryo.labels.raw[, , number]), col = seg.colors)
}

compareToRaw(15)

#Comment:
# As we can observe here, we think compared to those plots using raw labels, plots with
# expert labels seems be able to catch more information around the edge. Hence we think 
# use experts labels as our response is a good strategy


# Feature engineering for fiji features
# map first Fiji feature of the first embryo onto super pixels by taking median of raw value within
compareToFiji <- function(number, no.feature, FUN = median){ 
  # a function that plots embryos with expert labels and with fiji feature's values mapped to superpixels
  # with designated function for comparison
  # Args:
  #  number: which number of image are going to be compared
  #  no.feature: which number of fiji feature are going to be compared
  #  FUN: function to map fiji values back to superpixel, median is default
  mat.data <- readMat(sprintf(paste0(fiji.directory, 'emb%s.mat'), number))
  fiji.features <- mat.data$fiji
  fiji.mapped <- mapRawFiji(embryo.seg[,,number], fiji.features[,,no.feature], FUN)
  fiji.seg <- plotOnSP(embryo.seg[,,number], fiji.mapped)
  par(mfrow = c(1, 2))
  image(flipImg(embryo.labels[ , , number]), col = seg.colors)
  image(flipImg(fiji.seg), col = fiji.colors)
}

compareToFiji(5, 13, median)
compareToFiji(5, 13, max)
compareToFiji(5, 13, min)
compareToFiji(5, 13, var)
compareToFiji(12, 2, skewness)
compareToFiji(12, 2, kurtosis)

#Comment:
# We observed that it seems that mapped fiji features sometimes can roughly get edge information like 
# image 15 with 13th fiji feature, but it is actually not consistent across image which makes those fiji
# features look messy and noisy.

# compare super pixel feature values by organ
compareOrgan <- function(number, no.feature, FUN = median){
  # a function that plots density of organs with designated no.feature to 
  # see if the feature can seperate the density well
  # Args:
  #  number: which number of image
  #  no.feature: which number of fiji feature
  #  FUN: function to map fiji values back to superpixel, median is default
  start <- 0
  end <- 0
  if (number == 1) {
    start <- 1
    end <- num.sp[1]
  }else{
    for(i in 1:number){
      end <- end + num.sp[i] 
    }
    start = end - num.sp[i] + 1
  }
  embryo.organ <- as.factor(organ.class[start:end])
  mat.data <- readMat(sprintf(paste0(fiji.directory, 'emb%s.mat'), number))
  fiji.features <- mat.data$fiji
  fiji.mapped <- mapRawFiji(embryo.seg[,,number], fiji.features[,,no.feature], FUN)
  embryo.data <- data.frame(fiji = fiji.mapped, organ = embryo.organ)
  ggplot(embryo.data, aes(x = fiji, fill = organ)) + geom_density(alpha = 0.3)
}
compareOrgan(12, 2, skewness)
compareOrgan(12, 2, kurtosis)
#Comment:
# we observed that those fiji features cannot really seperate those organs well,
# they are all almost overlapped together.



#Try PCA on fiji features
mat.data <- readMat(paste0(fiji.directory, 'emb1.mat'))
fiji.features <- mat.data$fiji
fiji.data <- matrix(0, 180000, 13)

for(i in 1:13){
  fiji.data[, i] <- as.numeric(fiji.features[,,i])
}

fiji.data <- as.data.frame(fiji.data)

fiji.pca <- prcomp(fiji.data,
                 center = TRUE,
                 scale. = TRUE)
plot(fiji.pca, type = "l", main = NULL)
#Take first 3 components based on the plot
fiji.pca.data <- fiji.pca$x[,1:3]
fiji.features.pca <- array(fiji.pca.data, c(300, 600, 3))

embryo.organ <- as.factor(organ.class[1:num.sp[1]])
fiji.pca.mapped1 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,1])
fiji.pca.mapped2 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,2])
fiji.pca.mapped3 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,3])

embryo.pca.data <- data.frame(fiji.pca1 = fiji.pca.mapped1,
                              fiji.pca2 = fiji.pca.mapped2,
                              fiji.pca3 = fiji.pca.mapped3,
                              organ = embryo.organ)
ggplot(embryo.pca.data, aes(x = fiji.pca1, fill = organ)) + geom_density(alpha = 0.3)

plot3d(embryo.pca.data[,1:3], col = embryo.pca.data$organ)





#Try PCA on log fiji features by shifting the original fiji features
mat.data <- readMat(paste0(fiji.directory, 'emb1.mat'))
fiji.features <- mat.data$fiji
fiji.data <- matrix(0, 180000, 13)

for(i in 1:13){
  fiji.data[, i] <- as.numeric(fiji.features[,,i])
}

fiji.data <- as.data.frame(fiji.data)
toadded <- -as.numeric(apply(fiji.data, 2, min)) + 0.1

#Add each column by its maximum to make everthing is positive
fiji.data.shift <- sweep(fiji.data, 2, toadded, FUN = "+")
log.fiji.shift <- log(fiji.data.shift)
ggplot(log.fiji.shift, aes(x = V12)) + geom_density(alpha = 0.3)

fiji.pca <- prcomp(log.fiji.shift,
                   center = TRUE,
                   scale. = TRUE)
plot(fiji.pca, type = "l")
#Take first 3 components
fiji.pca.data <- fiji.pca$x[,1:3]
fiji.features.pca <- array(fiji.pca.data, c(300, 600, 3))

embryo.organ <- as.factor(organ.class[1:num.sp[1]])
fiji.pca.mapped1 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,1])
fiji.pca.mapped2 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,2])
fiji.pca.mapped3 <- mapRawFiji(embryo.seg[,,1], fiji.features.pca[,,3])

embryo.pca.data <- data.frame(fiji.pca1 = fiji.pca.mapped1,
                              fiji.pca2 = fiji.pca.mapped2,
                              fiji.pca3 = fiji.pca.mapped3,
                              organ = embryo.organ)
ggplot(embryo.pca.data, aes(x = fiji.pca1, fill = organ)) + geom_density(alpha = 0.3)

plot3d(embryo.pca.data[,1:3], col = embryo.pca.data$organ)


#Explore spatial features

distToCenter <- function(x){
  # A function to calculate superpixel's distance to center
  # Args:
  #  x: vector to represent coordinate of each superpixel
  return(sqrt((x[1] - 600)^2 + (x[2] - 300)^2))
} 

dist.to.center1 <- apply(as.data.frame(sp.coord[1]), 1, distToCenter)
dist.to.center2 <- apply(as.data.frame(sp.coord[2]), 1, distToCenter)

# compare super pixel feature values by organ
plotClass <- function(n){
  # A function to plot density for each organ for boundary distance and center distance, also
  # scatter plot to see how organs spread between boundary distance and centter distance.
  # Args:
  #  n: which number of image
  start <- 0
  end <- 0
  dist.to.center <- apply(as.data.frame(sp.coord[n]), 1, distToCenter)
  if ( n == 1) {
    start <- 1
    end <- num.sp[1]
  }else{
    for(i in 1:n){
      end <-  end + num.sp[i] 
    }
    start = end - num.sp[i] + 1
  }
  embryo.organ <- as.factor(organ.class[start:end])
  embryo.data <- data.frame(cent.dist = dist.to.center, boundary.dist = boundary.distance[start:end],
                             organ = embryo.organ)
  
  plot1 <- ggplot(embryo.data, aes(x = boundary.dist, fill = organ)) + geom_density(alpha = 0.3)
  plot2 <- ggplot(embryo.data, aes(x = cent.dist, fill = organ)) + geom_density(alpha = 0.3)
  plot3 <- ggplot(embryo.data, aes(x = cent.dist, y = boundary.dist, color = organ)) + geom_point(alpha = 0.5)
  grid.arrange(plot1, plot2, plot3, ncol=3)
}


#We can observe that spatial data classify organs better than fiji features, 
#however, some of area which is hardly classified by spatial data as well. 
# We doubt those areas are located on the edges of those organs.


embryoDataGenerator <- function(n){
  # A function to generate embryo data with spatial features
  # Args:
  #  n: which number of image
  start <- 0
  end <- 0
  dist.to.center <- apply(as.data.frame(sp.coord[n]), 1, distToCenter)
  if ( n == 1) {
    start <- 1
    end <- num.sp[1]
  }else{
    for(i in 1:n){
      end <-  end + num.sp[i] 
    }
    start = end - num.sp[i] + 1
  }
  embryo.organ <- as.factor(organ.class[start:end])
  embryo.data <- data.frame(cent.dist = dist.to.center, boundary.dist = boundary.distance[start:end],
                            organ = embryo.organ)
  return(embryo.data)
}  


#image 2 
embryo.data2 <- embryoDataGenerator(2)
ggplot(embryo.data2, aes(x = cent.dist, y = boundary.dist, color = organ)) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = c(110, 160))

less160 <- which(embryo.data2$boundary.dist < 160)
large110 <- which(embryo.data2$boundary.dist > 110)
suspect <- intersect(less160, large110)
marker <- sort(unique(as.numeric(embryo.seg[,,2])))[-1]
suspect.marker <- marker[suspect]
seg2 <- embryo.labels[,,2]
seg2[which(as.numeric(embryo.seg[,,2]) %in% suspect.marker)] = 4

par(mfrow = c(1, 2))

image(flipImg(embryo.labels[,,2]), col = seg.colors[1:5])
image(flipImg(seg2), col = c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#BEBADA"  ))

  #image 15
embryo.data15 <- embryoDataGenerator(15)
ggplot(embryo.data15, aes(x = cent.dist, y = boundary.dist, color = organ)) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = c(120, 160))

less160 <- which(embryo.data15$boundary.dist < 160)
large120 <- which(embryo.data15$boundary.dist > 120)
suspect <- intersect(less160, large120)
marker <- sort(unique(as.numeric(embryo.seg[,,15])))[-1]
suspect.marker <- marker[suspect]
seg15 <- embryo.labels[,,15]
seg15[which(as.numeric(embryo.seg[,,15]) %in% suspect.marker)] = 4

par(mfrow = c(1, 2))

image(flipImg(embryo.labels[,,15]), col = seg.colors[1:5])
image(flipImg(seg15), col = c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#BEBADA"  ))


#image 117
embryo.data117 <- embryoDataGenerator(117)
ggplot(embryo.data117, aes(x = cent.dist, y = boundary.dist, color = organ)) + geom_point(alpha = 0.5) +
  geom_hline(yintercept = c(130, 160))

less160 <- which(embryo.data117$boundary.dist < 160)
large130 <- which(embryo.data117$boundary.dist > 130)
suspect <- intersect(less160, large130)
marker <- sort(unique(as.numeric(embryo.seg[,,117])))[-1]
suspect.marker <- marker[suspect]
seg117 <- embryo.labels[,,117]
seg117[which(as.numeric(embryo.seg[,,117]) %in% suspect.marker)] = 4

par(mfrow = c(1, 2))

image(flipImg(embryo.labels[,,117]), col = seg.colors[1:5])
image(flipImg(seg117), col = c("#8DD3C7", "#FFFFB3", "#FB8072", "#80B1D3", "#BEBADA"  ))


#Comment:
# After taking a look a few graph, we find that those superpixels that spatial features cannot 
# seperate well almost around the edge of those organs



#Looking into the difference and matching between raw labels and expert labels
diff.mat <- matrix(0, 150, 7)
for(i in 1:150){
  diff.table <- table(abs(as.numeric(embryo.labels.raw[,,i])^2 - as.numeric(embryo.labels[,,i])^2))
  name <- as.numeric(names(diff.table))
  val <- as.vector(diff.table)
  for(j in 1:length(name)){
    if(name[j] < 2){
      diff.mat[i, name[j] + 1] <- val[j]  
    }
    else if(name[j] < 6){
      diff.mat[i, name[j]] <- val[j]
    }else{
      diff.mat[i, name[j] - 2] <- val[j]
    }
  }
}

#To see how many percentage match between super pixels labels and raw labels
match.vec <- diff.mat[,1]/180000
plot(density(match.vec), main = "", xlab = "percentage match between super pixels labels and raw labels")

diff.percent.mat <- sweep(diff.mat[, -1], 1, apply(diff.mat[, -1], 1, sum), FUN = "/")
diff.percent.df <- data.frame(diff.percent.mat)
names(diff.percent.df) <- c("0x1", "1x2","0x2", "2x3", "1x3", "0x3")

#melt the data in order to plot density plots
val.vec <- as.numeric(diff.percent.mat)
class.vec <- rep(c("0-1", "1-2","0-2", "2-3", "1-3", "0-3"), each = 150)
percent.class <- data.frame(percentage = val.vec, class = class.vec )
ggplot(percent.class, aes(x = percentage, fill = class)) + geom_density(alpha = 0.3)


val.vec <- as.numeric(diff.percent.mat[, -3])
class.vec <- rep(c("0-1", "1-2", "2-3", "1-3", "0-3"), each = 150)
percent.class <- data.frame(percentage = val.vec, class = class.vec )
ggplot(percent.class, aes(x = percentage, fill = class)) + geom_density(alpha = 0.3)

val.vec <- as.numeric(diff.percent.mat[, c(-1, -3, -6)])
class.vec <- rep(c("1-2", "2-3", "1-3"), each = 150)
percent.class <- data.frame(percentage = val.vec, class = class.vec )
ggplot(percent.class, aes(x = percentage, fill = class)) + geom_density(alpha = 0.3)
