summarynz=function(object, EPS=0){
    degrees = object$degrees
    p = length(degrees)
    counter = rep(1:p, degrees)
    counter = diag(p)[counter, ]
    isalpha=abs(object$alpha) > EPS
    isbeta = (t(counter) %*% abs(object$beta))>EPS
    numnlin = apply(isbeta, 2, sum)
    islin= isalpha& !isbeta
    numlin = apply(islin, 2, sum)
    numnz = apply(isalpha | isbeta, 2, sum)
    cbind(Linear = numlin, Nonlinear = numnlin, Nonzero = numnz)
}

