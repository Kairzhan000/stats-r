bubleSort <- function(data){
  n <- length(data)
  if(n < 2) return (data)
  for (i in 1:(n-1)){
    for(j in 1:(n-i)){
      if(data[j] > data[j+1]){
        temp <- data[j]
        data[j] <- data[j+1]
        data[j+1] <- temp
      }
    }
  }
  return (data)
}
#_____________


CFmomentOfOrder <- function(data, r, c = 0){
  n <- length(data)
  return (sum((data - c)^r) / n)
}

CFmean <- function(data){
  answer <- CFmomentOfOrder(data, 1)
  return (answer)
}

CFvariance <- function(data,ownMean){
  answer <- CFmomentOfOrder(data,2,ownMean)
  return (answer)
}

CFsd <- function(data){
  return (sqrt(CFvariance(data)))
}

CFquasiVariance <- function(data){      
  n <- length(data)
  return (CFvariance(data) * n / (n - 1))
}
RMS <- function(data){
  n <- length(data)
  answer <- (sum(data^2) / n)^(1/2)
  return (answer)
}
CFweightedMean <- function(data, weights){
  return (sum(data * weights) / sum(weights))
}




findMode <- function(data){
  n <- length(data)
  if(n == 0) return (NULL)
  values <- c()
  counts <- c()
  for (i in 1:n){
    found <- FALSE
    if(length(values) > 0){
      for(j in 1:length(values)){
        if(data[i] == values[j]){
          counts[j] <- counts[j] + 1
          found <- TRUE
          break
        }
      }
    }
    if(!found){
      values <- c(values, data[i])
      counts <- c(counts, 1)
    }
  }
  maxCount <- counts[1]
  for(j in 1:length(counts)){
    if(counts[j] > maxCount) maxCount <- counts[j]
  }
  modes <- c()
  for(j in 1:length(values)){
    if(counts[j] == maxCount) modes <- c(modes, values[j])
  }
  return (modes)
}

CFpercentiles <- function(data,percent){
  n <- length(data)
  if(n == 0) return (NULL)
  if(percent < 0 || percent > 100)return (NULL)
  x <- bubleSort(data)
  c <- percent / 100
  cN <- c * n
  E <- floor(cN)
  D <- cN - E
  if(D != 0){
    return (x[E+1])
  }else{
    if(E == 0) return(x[1])
    if(E == n)return (x[n])
    return ((x[E] + x[E+1])/2)
  }
}

CFmedian <- function(data) CFpercentiles(data, 50)
CFrange  <- function(data) max(data) - min(data)
CFiqr    <- function(data) CFpercentiles(data, 75) - CFpercentiles(data, 25)

CFmeanDeviation <- function(data, p = CFmean(data)){
  return (CFmomentOfOrder(abs(data - p), 1))
}

CFmse <- function(data, p){
  return (CFmomentOfOrder(data, 2, p))
}

#_______________

CFmomentTable <- function(xi, ni, r, c = 0){
  return (sum((xi - c)^r * ni) / sum(ni))
}

CFmeanTable <- function(xi, ni){
  return (CFmomentTable(xi, ni, 1))
}

CFvarianceTable <- function(xi, ni){
  return (CFmomentTable(xi, ni, 2, CFmeanTable(xi, ni)))
}

CFsdTable <- function(xi, ni){
  return (sqrt(CFvarianceTable(xi, ni)))
}

CFquasiVarianceTable <- function(xi, ni){
  N <- sum(ni)
  return (CFvarianceTable(xi, ni) * N / (N - 1))
}

CFcvTable <- function(xi, ni, percent = FALSE){
  answer <- CFsdTable(xi, ni) / abs(CFmeanTable(xi, ni))
  if(percent) return (answer * 100)
  return (answer)
}

CFmodeTable <- function(xi, ni){
  return (xi[ni == max(ni)])
}

CFcumFrequencyTable <- function(ni, relative = FALSE){
  if(relative) return (cumsum(ni) / sum(ni))
  return (cumsum(ni))
}

CFpercentilesTable <- function(xi, ni, percent){
  return (CFpercentiles(rep(xi, ni), percent))
}

CFmedianTable <- function(xi, ni) CFpercentilesTable(xi, ni, 50)

CFmeanDeviationTable <- function(xi, ni, p = CFmeanTable(xi, ni)){
  return (CFmomentTable(abs(xi - p), ni, 1))
}

CFskewnessFisherTable <- function(xi, ni){
  mu3 <- CFmomentTable(xi, ni, 3, CFmeanTable(xi, ni))
  return (mu3 / CFsdTable(xi, ni)^3)
}

CFkurtosisTable <- function(xi, ni){
  mu4 <- CFmomentTable(xi, ni, 4, CFmeanTable(xi, ni))
  return (mu4 / CFvarianceTable(xi, ni)^2 - 3)
}

CFskewnessPearsonTable <- function(xi, ni){
  Mo <- CFmodeTable(xi, ni)[1]
  return ((CFmeanTable(xi, ni) - Mo) / CFsdTable(xi, ni))
}

CFclassMarks <- function(L){
  k <- length(L) - 1
  return ((L[1:k] + L[2:(k+1)]) / 2)
}

CFwidths <- function(L){
  k <- length(L) - 1
  return (L[2:(k+1)] - L[1:k])
}
CFgroupedPercentile <- function(L, ni, percent){
  N <- sum(ni)
  cN <- round(percent / 100 * N, 10)
  Ncum <- cumsum(ni)
  if(cN >= N) return (L[length(L)])
  i <- 1
  while(Ncum[i] <= cN) i <- i + 1
  if(i == 1) prev <- 0 else prev <- Ncum[i-1]
  a <- L[i+1] - L[i]
  return (L[i] + (cN - prev) / ni[i] * a)
}

CFmodalInterval <- function(L, ni){
  h <- ni / CFwidths(L)
  i <- 1
  for(j in 1:length(h)){
    if(h[j] > h[i]) i <- j
  }
  return (c(lower = L[i], upper = L[i+1], height = h[i]))
}
CFgroupedMean     <- function(L, ni) CFmeanTable(CFclassMarks(L), ni)
CFgroupedVariance <- function(L, ni) CFvarianceTable(CFclassMarks(L), ni)
CFgroupedSd       <- function(L, ni) sqrt(CFgroupedVariance(L, ni))