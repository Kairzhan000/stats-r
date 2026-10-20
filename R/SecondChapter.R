
tableMoment <- function(table, xValues, yValues, r, s, a = 0, b = 0){
  N <- sum(table)
  answer <- 0
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      answer <- answer + table[i, j] * (xValues[i] - a)**r * (yValues[j] - b)**s
    }
  }
  answer <- answer / N
  return (unname(answer))
}

tableMeanX <- function(table, xValues, yValues){
  return (tableMoment(table, xValues, yValues, 1, 0))
}

tableMeanY <- function(table, xValues, yValues){
  return (tableMoment(table, xValues, yValues, 0, 1))
}

tableVarianceX <- function(table, xValues, yValues){
  m20 <- tableMoment(table, xValues, yValues, 2, 0)
  return (m20 - tableMeanX(table, xValues, yValues)**2)
}

tableVarianceY <- function(table, xValues, yValues){
  m02 <- tableMoment(table, xValues, yValues, 0, 2)
  return (m02 - tableMeanY(table, xValues, yValues)**2)
}

tableCovariance <- function(table, xValues, yValues){
  m11 <- tableMoment(table, xValues, yValues, 1, 1)
  meanX <- tableMeanX(table, xValues, yValues)
  meanY <- tableMeanY(table, xValues, yValues)
  return (m11 - meanX * meanY)
}

tableMarginalX <- function(table){
  return (unname(rowSums(table)))
}
tableConditionalX <- function(table, j){
  return (unname(table[, j] / sum(table[, j])))
}

tableMarginalY <- function(table){
  return (unname())
}
tableConditionalY <- function(table, i){
  return (unname(table[i, ] / sum(table[i, ])))
}

tableIsIndependent <- function(table){
  N <- sum(table)
  rows <- tableMarginalX(table)
  cols <- tableMarginalY(table)
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      if(abs(table[i,j] / N  - (rows[i] / N) * (cols[j] \ N)) >1e-9){
        return (FALSE)
      }
    }
  }
  return (TRUE)
}

tableCorrelation <- function(table, xValues, yValues){
  cov <- tableCovariance(table, xValues, yValues)
  varX <- tableVarianceX(table, xValues, yValues)
  varY <- tableVarianceY(table, xValues, yValues)
  return (cov / sqrt(varX * varY))
}

tableRegressionYX <- function(table, xValues, yValues){
  b <- tableCovariance(table, xValues, yValues) / tableVarianceX(table, xValues, yValues)
  a <- tableMeanY(table, xValues, yValues) - b * tableMeanX(table, xValues, yValues)
  return (c(a = a, b = b))
}

tableRegressionXY <- function(table, xValues, yValues){
  b <- tableCovariance(table, xValues, yValues) / tableVarianceY(table, xValues, yValues)
  a <- tableMeanX(table, xValues, yValues) - b * tableMeanY(table, xValues, yValues)
  return (c(a = a, b = b))
}

tableMSEYX <- function(table, xValues, yValues){
  coef <- tableRegressionYX(table, xValues, yValues)
  N <- sum(table)
  answer <- 0
  for(i in 1:nrow(table)){
    for(j in 1:ncol(table)){
      e <- yValues[j] - (coef["a"] + coef["b"] * xValues[i])
      answer <- answer + table[i, j] * e^2
    }
  }
  return (unname(answer / N))
}

tableR2YX <- function(table, xValues, yValues){
  mse <- tableMSEYX(table, xValues, yValues)
  return (1 - mse / tableVarianceY(table, xValues, yValues))
}

line_reg <- function(M,Y){
  return(solve(t(M) %*% M)%*% t(M) %*% Y )
}