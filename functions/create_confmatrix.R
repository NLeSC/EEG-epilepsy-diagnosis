create_confmatrix = function(predicted,reference) {
  makeconfmatsquare = function(x) { # function to turn confusion 1x2 or 2x1 table() matrix into a square matrix
    ihave = which(colnames(x) == rownames(x))
    if (length(ihave) == 0) {
      print(colnames(x))
      print(rownames(x))
    }
    if (ncol(x) > nrow(x)) {
      if (ihave == 2) {
        x = rbind(c(0,0),x)
      } else {
        x = rbind(x,c(0,0))
      }
    } else if (ncol(x) < nrow(x)) {
      if (ihave == 2) {
        x = cbind(c(0,0),x)
      } else {
        x = cbind(x,c(0,0))
      }
    }
    return(x)
  }
  confmat = makeconfmatsquare(table(predicted,reference))
  return(confmat)
}