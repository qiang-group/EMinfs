#' Generate I-Spline Basis Functions
#'
#' Computes the basis functions for I-splines, which are integrated M-splines.
#' I-splines are non-decreasing functions useful for modeling cumulative hazard functions.
#' This function is based on the work of Ramsay (1988).
#'
#' @param x A numeric vector of values at which to evaluate the spline basis.
#' @param order An integer specifying the order of the spline (e.g., 2 for quadratic).
#' @param knots A numeric vector of knot locations for the spline.
#'
#' @return A matrix where each column is an I-spline basis function evaluated at the points in `x`.


Ispline <-
  function(x,order,knots){
    # M Spline function with order k=order+1. or I spline with order
    # x is a row vector
    # k is the order of I spline
    # knots are a sequence of increasing points
    # the number of free parameters in M spline is the length of knots plus 1.

    k=order+1
    m=length(knots)
    n=m-2+k # number of parameters
    t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots

    yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
    for (l in k:n){
      yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
    }

    yytem1=yy1
    for (ii in 1:order){
      yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
      for (i in (k-ii):n){
        yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
      }
      yytem1=yytem2
    }

    index=rep(0,length(x))
    for (i in 1:length(x)){
      index[i]=sum(t<=x[i])
    }

    yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))

    if (order==1){
      for (i in 2:n){
        yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
      }
    }else{
      for (j in 1:length(x)){
        for (i in 2:n){
          if (i<(index[j]-order+1)){
            yy[i-1,j]=1
          }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
            yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
          }else{
            yy[i-1,j]=0
          }
        }
      }
    }
    return(yy)
  }
