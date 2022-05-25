predict.far2 <- function (object, newdata = NULL, label = NULL) 
{
  if ((!is.null(class(object))) && (class((object)) != "far")) 
    stop("object is not of class far")
  if ((!is.null(class(newdata))) && (class((newdata)) != "fdata")) 
    stop("newdata is not of class fdata")
  y <- object$y
  ny <- length(y)
  r <- ny
  n <- ncol(newdata[[object$y[1]]])
  if(is.null(label)) label <- c(colnames(newdata[[y[1]]])[-1])
  data <- list()
  for (i in 1:ny) {
        data[[y[i]]] <- sweep((newdata[[y[i]]]), 1, object$databar[[y[i]]], "-")
      }
  class(data) <- "fdata"
  kn <- object$kn
  listobs <- rep(TRUE, n)
  nbobs <- sum(listobs == TRUE)
  datacent <- matrix(0, ncol = nbobs, nrow = sum(kn))
  kkn <- c(0, kn)
  for (k in (1:r)) {
      datacent[sum(kkn[1:k]) + (1:kkn[k + 1]), ] <- t(object$v[[k]]) %*% 
        ((data[[k]])[, listobs, drop = FALSE])
    }
    pred <- list()
    length(pred) <- ny
    for (i in (1:ny)) {
      pred[[i]] <- object$v[[i]] %*% (object$rho %*% datacent)[sum(kkn[1:i]) + 
                                                                 (1:kkn[i + 1]), , drop = FALSE]
    }
  for (i in (1:ny)) {
    if (!is.null(object$databar)) 
      pred[[i]] <- sweep((pred[[i]]), 1, object$databar[[i]], "+")
    rownames(pred[[i]]) <- rownames(data[[i]])
    colnames(pred[[i]]) <- label[listobs]
  }
  names(pred) <- object$y
  class(pred) <- "fdata"
  return(pred)
}


