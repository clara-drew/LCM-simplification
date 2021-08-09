poly <- function(x,cov){
  p <- c(x^4,x^3,x^2,x,1)
  return(t(cov) %*% p)
}

t2_poly <- function(maxT){
  
  min <- (19/32)
  last <- 1
  y1 <- 0.875
  y2 <- 0.3
  y3 <- 0.125
  
  A = rbind(c(0,0,0,1,0),c(4*min^3,3*min^2,2*min,1,0),c(0,0,0,0,1),
            c(min^4,min^3,min^2,min,1),c(last^4,last^3,last^2,last,1))
  
  b = c(0,0,y1,y2,y3)
  
  
  cov <- solve(A) %*% b
  
  tps <- ((1:maxT)-1)/(maxT-1)
  
  return(unlist(lapply(tps,poly,cov=cov)))
  
}


