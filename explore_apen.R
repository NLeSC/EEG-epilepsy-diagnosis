rm(list=ls())
graphics.off()

library(pracma)
# apen = function (ts, edim = 2, r = 0.2 * sd(ts), elag = 1) {
#   N <- length(ts)
#   result <- numeric(2)
#   for (j in 1:2) {
#     m <- edim + j - 1
#     phi <- zeros(1, N - m + 1)
#     dataMat <- zeros(m, N - m + 1)
#     for (i in 1:m) dataMat[i, ] <- ts[i:(N - m + i)]
#     for (i in 1:(N - m + 1)) {
#       tempMat <- abs(dataMat - repmat(dataMat[, i, drop = FALSE], 
#                                       1, N - m + 1))
#       boolMat <- apply(tempMat > r, 2, max)
#       phi[i] <- sum(!boolMat)/(N - m + 1)
#     }
#     result[j] <- sum(phi)/(N - m + 1)
#   }
#   apen <- log(result[1]/result[2])
#   return(apen)
# }

apen2 = function (ts, edim = 2, r = 0.2 * sd(ts), elag = 1) {
  N <- length(ts)
  result <- numeric(2)
  for (j in 1:2) {
    m <- edim + j - 1
    phi <- zeros(1, N - m + 1)
    dataMat <- zeros(m, N - m + 1)
      # k0 = Sys.time()  
      for (i in 1:m) dataMat[i, ] <- ts[i:(N - m + i)]
#     k0 = Sys.time()  
#     for (i in 1:(N - m + 1)) {
#       tempMat <- abs(dataMat - repmat(dataMat[, i, drop = FALSE], 
#                                       1, N - m + 1))
#       boolMat <- apply(tempMat > r, 2, max)
#       phi[i] <- sum(!boolMat)/(N - m + 1)
#     }
    k1 = Sys.time()
    
    nn = N - m + 1
#     for (i in 1:nn) {
#       test2 = matrix(dataMat[,i,drop=FALSE],m,nn)
#       tempMat <- abs(dataMat - test2)
#       boolMat <- apply(tempMat > r, 2, max)
#       phi[i] <- sum(!boolMat)/nn
#     }
    myfun = function(x){
      tempMat <- abs(dataMat - matrix(as.matrix(x),m,nn))
      boolMat <- apply(tempMat > r, 2, max)
      phi[i] <- sum(!boolMat)/nn
    }
    phi = apply(dataMat,2,myfun)
    k2 = Sys.time()
    result[j] <- sum(phi)/(N - m + 1)
    # print(k1-k0)
    print(k2-k1)
  }
  apen <- log(result[1]/result[2])
  return(apen)
}


x = randn(n=1000,m=1)
b = sin((1:1000)/20)
x = x + b
# x11()
# plot(x)

t0 = Sys.time()
print(sample_entropy(x))
t1 = Sys.time()
print(apen2(x))
t2 = Sys.time()

cat("approx_entropy:\n")
print(t1-t0)
cat("myversion:\n")
print(t2-t1)