TECAxx <- function(xx,startN,frequencyN){
  uk <- ts(xx,start=c(startN , 1) , frequency = frequencyN)
  return (uk)
}
autoCorelation <- function(data,n){
  xx1 <- head(data,-n)
  xx2 <- tail(data,-n)
  
  sdXx1 <- xx1 - mean(xx1)
  sdXx2 <- xx2 - mean(xx2)
  
  answer <- sum(sdXx1 * sdXx2)/(sqrt(sum(sdXx1^2) * sum(sdXx2^2)))
  return (answer)
}
movingAverage <- function(series, k){
  if(k %% 2 == 1){
    w <- rep(1, k) / k
  } else {
    w <- c(0.5, rep(1, k - 1), 0.5) / k
  }
  return (stats::filter(series, w))
}