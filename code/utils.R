plotOnSP <- function(seg, features) {

  # Function for plotting a vector of features onto an embryo segmented into
  # superpixels
  # Args:
  #   seg: a matrix of integers, entries being superpixel segmentation labels
  #   features: a numeric vector with length equal to the number of superpixels
  #   in seg, entries being the feature values to plot
  # Return:
  #   a numeric matrix of the same size as seg with the values of features
  #   mapped to corresponding super pixels
  sp.values <- sort(unique(seg[seg != 0]))
  if (length(sp.values) != length(features)) stop('feature length mismatch')

  template <- matrix(0, nrow=nrow(seg), ncol=ncol(seg))
  
  for (s in sp.values) {
    sp.idcs <- which(seg == s, arr.ind=TRUE)
    template[sp.idcs] <- features[sp.values == s]
  }

  return(template)
}

mapRawFiji <- function(seg, raw.fiji, FUN=median) {

  # Function for mapping raw fiji data to superpixels
  # Args:
  #  seg: a matrix of integers, entries being superpixel labels
  #  raw.fiji: a numeric matrix of the same dimensions as raw.fiji
  #  FUN: the function used to map raw fiji data within a superpixel to a single
  #  value representing that sp
  # Return:
  #  a numeric vector with one entry per super pixel giving value of the mapped
  #  raw feature
  sp.values <- sort(unique(seg[seg != 0]))
  if (!all.equal(dim(seg), dim(raw.fiji))) stop('matrix dimension mismatch')

  fiji.features <- sapply(sp.values, function(s) {
    sp.idcs <- which(seg == s, arr.ind=TRUE)
    sp.value <- FUN(raw.fiji[sp.idcs])
    return(sp.value)
  })

  return(fiji.features)
}

flipImg <- function(img) {

  # Function to reorient matrix for plotting
  # Args:
  #   img: a numeric matrix, entries being pixel values
  # Return:
  #   a numeric matrix of the same size as img, oriented correctly for R's image
  #   function
  img.t <- t(img)
  img.flipped <- img.t[, seq(ncol(img.t), 1, by=-1)]
  return(img.flipped)
}
