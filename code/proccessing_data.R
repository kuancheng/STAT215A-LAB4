library(R.matlab)
setwd("/Users/fadikfoury/Downloads/Lab4/lab")
fiji.directory <- './raw_fiji/'

#col: embryoID (1), super_pixel_id(2), x,y(4),  figi(16)
Y <- matrix( nrow = 30600000, ncol = 17)
i = 1
#loop accross every pixel in the the embryo.segments.
for (embryo_id in 1:170)
{
  #Y[1,] = x
  #load the embryo fiji
  mat.data <- readMat(paste0(fiji.directory, paste("emb",embryo_id,".mat", sep = "") ))
  fiji.features <- mat.data$fiji
  
  for (x_i in 1:600)
  {
    
    for (y_i in 1:300)
    {
      
    super_pixel_id = embryo.seg[y_i,x_i,embryo_id]
    
    x <- vector( length = 17)
    x[1] = embryo_id
    x[2] = super_pixel_id
    x[3] = x_i
    x[4] = y_i
    
      for(figi_id in 1:13)
      {
        x[4+figi_id] = fiji.features[y_i,x_i,figi_id]
        
      }
      
    Y[i,] = x
    i = i + 1
    #i is the index of the current pixel
    }
  }
  
  print(paste("Finished" ,embryo_id,"Embryo") )
}
setwd("/Users/fadikfoury/Downloads/Lab4/lab")
fiji.directory <- './raw_fiji/'
mat.data <- readMat(paste0(fiji.directory, paste("emb",1,".mat", sep = "") ))

