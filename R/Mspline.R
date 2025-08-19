#' Generate M-Spline Basis Functions
#'
#' Computes the basis functions for M-splines, which are non-negative splines
#' useful for modeling hazard functions. [cite_start]This function is based on the work of Ramsay (1988)[cite: 95].
#'
#' @param x A numeric vector of values at which to evaluate the spline basis.
#' @param order An integer specifying the order of the spline (e.g., 2 for quadratic).
#' @param knots A numeric vector of knot locations for the spline.
#'
#' @return A matrix where each column is an M-spline basis function evaluated at the points in `x`.

Mspline<-function(x,order,knots){
  # get M spline and I spline matrix with order
  # x is a row vector
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.


  ### get Mspline bases ###
  k1=order
  m=length(knots)
  n1=m-2+k1 # number of parameters
  t1=c(rep(1,k1)*knots[1], knots[2:(m-1)], rep(1,k1)*knots[m]) # newknots

  tem1=array(rep(0,(n1+k1-1)*length(x)),dim=c(n1+k1-1, length(x)))
  for (l in k1:n1){
    tem1[l,]=(x>=t1[l] & x<t1[l+1])/(t1[l+1]-t1[l])
  }

  if (order==1){
    mbases=tem1
  }else{
    mbases=tem1
    for (ii in 1:(order-1)){
      tem=array(rep(0,(n1+k1-1-ii)*length(x)),dim=c(n1+k1-1-ii, length(x)))
      for (i in (k1-ii):n1){
        tem[i,]=(ii+1)*((x-t1[i])*mbases[i,]+(t1[i+ii+1]-x)*mbases[i+1,])/(t1[i+ii+1]-t1[i])/ii
      }
      mbases=tem
    }
  }

  return(mbases)
}
