
meansfunc <- function(x,y, c1 = 0.01){
  levels <- levels(factor(c(x,y)))
  x = factor(x, levels = levels)
  y = factor(y, levels = levels)
  tablexy <- table(x,y)
  tablex <- rowSums(tablexy)
  tabley <- colSums(tablexy)

  (2*sum(tablex*tabley) + c1)/(sum(tablex^2) + sum(tabley^2) + c1)
}

gini <- function(x){
 1 - (sum(table(x)^2)/length(x)^2)
}

ginicorr <- function(x){
  k <- length(table(x))
  if(k > 1){
  gini(x)/(1-1/k)
  } else gini(x)
}


sqrtgini <-  function(x){
  1 - sqrt(sum(table(x)^2)/length(x)^2)
}

sqrtginicorr <- function(x){
  k <- length(table(x))
  if(k > 1){
  sqrtgini(x)/(1-1/sqrt(k))
  } else sqrtgini(x)
}

cfunc <- function(x, y, c2 = 0.01){
  varx <- ginicorr(x)
  vary <- ginicorr(y)

  (2*sqrt(varx * vary)+c2)/(varx + vary + c2)

}

sfunc <- function(x, y){
  # should take care of this beforehand!
  numx <- as.numeric(x)
  numy <- as.numeric(y)
  EMCluster::RRand(numx, numy)$adjRand
}


binssim <- function(x, y, alpha = 1, beta = 1, gamma = 1, ...){
  (meansfunc(x, y,...)^alpha)*(cfunc(x, y, ...)^beta)*(sfunc(x, y)^gamma)
}

ssimcomponents <- function(x, y, alpha = 1, beta = 1, gamma = 1, ...){
  c((meansfunc(x, y,...)^alpha),(cfunc(x, y, ...)^beta),(sfunc(x, y)^gamma))
}


randmode <- function(x){
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  #sample(ux[tab == max(tab)],1)
  ux[which.max(tab)]
}

downsample_2d <- function(x){
  dims = dim(x)
  newdims = floor(dims/2)
  if(any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- matrix(nrow = newdims[1], ncol = newdims[2])
  for(i in 1:newdims[1]){
    for(j in 1:newdims[2]){
      xstart <- 2*i -1
      ystart <- 2*j - 1
      newx[i,j] = randmode(c(x[xstart:(xstart+1), ystart:(ystart+1)]))
    }
  }
  newx
}

downsample_3d_slice <- function(x){
  dims = dim(x)
  newdims = floor(c(dims[1:2]/2,dims[3]))
  if(any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- array(dim = newdims)
  for(i in 1:dims[3]){
    newx[,,i] = downsample_2d(x[,,i])
  }
  newx
}

downsample_3d_cube <- function(x){
  dims = dim(x)
  newdims = floor(dims/2)
  if(any(newdims < 1)){
    warning("Cannot downsample, at least one dimension is of length 1")
    return(x)
  }
  newx <- array(dim = newdims)
  for(i in 1:newdims[1]){
    for(j in 1:newdims[2]){
      for(k in 1:newdims[3]){
      xstart <- 2*i -1
      ystart <- 2*j - 1
      zstart <- 2*k - 1
      newx[i,j] = randmode(c(x[xstart:(xstart+1), ystart:(ystart+1), zstart:(zstart+1)]))
      }
    }
  }
  newx
}


catssim_2d <- function(x,y,...){
  if(any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels <- levels(factor(c(x,y)))
  dims = dim(x)
  nrow = dims[1]
  ncol = dims[2]
  resultmatrix = matrix(0, nrow = (nrow-7)*(ncol-7),3)
  for(i in 1:(nrow-7)){
    for(j in 1:(ncol-7)){
      place = j+(ncol-7)*(i-1)
      subx = x[i:(i+7),j:(j+7)]
      suby = y[i:(i+7),j:(j+7)]

      subx = factor(subx, levels = levels)
      suby = factor(suby, levels = levels)
      resultmatrix[place,] = ssimcomponents(subx, suby, ...)
    }
  }
  colMeans(resultmatrix)
  #resultmatrix
}



catmssim_2d <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), ...){
  if(any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels = length(weights)
  if(min(dim(x))<32) levels = min(c(floor(min(dim(x)/8)),levels))
  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)
  results[1,] = catssim_2d(x,y)

  if( levels > 1){
    for(i in 2:levels){
        x = downsample_2d(x)
        y = downsample_2d(y)
        results[i,] = catssim_2d(x,y)
      }
    }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults

}


catssim_3d_slice <- function(x,y,...){
  if(any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  dims = dim(x)
  sliceresults = matrix(0, nrow = dims[3], ncol = 3)
  for(i in 1:dims[3]) sliceresults[i,] = catssim_2d(x[,,i],y[,,i],...)
  colMeans(sliceresults)
}


catmssim_3d_slice <- function(x, y, weights = c(0.0448, 0.2856, 0.3001, 0.2363, 0.1333), ...){
  if(any(dim(x) != dim(y))) stop('x and y have nonconformable dimensions')
  levels = length(weights)
  dims = dim(x)
  if(length(dims) < 3) stop('x and y are not 3-dimensional')
  if(min(dims[1:2])<32) levels = min(c(floor(min(dims[1:2]/8)),levels))
  weights = weights[1:levels]
  results = matrix(0, nrow = levels, ncol = 3)

  results[1,] = catssim_3d_slice(x,y,...)

  if( levels > 1){
    for(j in 2:levels){
      x = downsample_3d_slice(x)
      y = downsample_3d_slice(y)
      results[j,] = catssim_3d_slice(x,y,...)
    }
  }

  # use "luminosity" only from top level, use C and S from all levels
  csresults <- prod(results[,2:3]^(weights))
  #print(results)

  (results[levels,1]^weights[levels])*csresults


}

