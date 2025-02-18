
K.cos <- function(s,t){
    k1 = abs(s-t)/2
    k2 = (s+t)/2
    ans1 = k1^4-2*k1^3+k1^2-1/30
    ans2 = k2^4-2*k2^3+k2^2-1/30
    
    return(-1/3 *(ans1+ans2))
}



K.gauss <- function(s,t, gamma){
    return(exp(-0.5 * (s-t)^2/ (gamma^2)))
}



################################
##### Sobolev RK function ######
################################

k1 <- function(t){
  return(t-.5)
}
k2 <- function(t){
  return( (k1(t)^2-1/12)/2 )
}
k4 <- function(t){
  return( (k1(t)^4-k1(t)^2/2+7/240)/24 )
}
K.sob <- function(s,t){
  ans <- 1 + k1(s)*k1(t) + k2(s)*k2(t) - k4(abs(s-t))
  return(ans)
}





#############################
#### svec transformation ####
#############################
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

svec <- function(X){
    # symmetric X of size n times n
    ind <- as.vector(row(X)> col(X))
    ind2 <- as.vector(row(X)>= col(X))
    X[ind] <- X[ind] * sqrt(2)
    return(X[ind2])
}

svec.inv <- function(x){
    nn <- length(x)
    n <- (-1 + sqrt(1+4 *2*nn))/2
    X <- matrix(nr=n, nc=n)
    ind <- as.vector(row(X)> col(X))
    ind2 <- as.vector(row(X)>= col(X))
    ind3 <- as.vector(row(X)< col(X))
    X[ind2] <- x
    X[ind] <- X[ind]/sqrt(2)
    X[ind3] <- t(X)[ind3]
    return(X)
}

smat <- function(X){
    # symmetric X of size n^2 times n^2
    n <- sqrt(dim(X)[1])
    if (!is.wholenumber(n)) stop("n is not an integer")
    temp <- matrix(1:(n*n), nr=n, nc=n)
    index <- temp[as.vector(row(temp)> col(temp))]
    index2 <- t(temp)[as.vector(row(temp)> col(temp))]
    ind2 <- as.vector(row(temp)>= col(temp))
    Y <- X
    Y[,index]  <- (Y[,index] + Y[,index2])/sqrt(2)
    Y[index,]  <- (Y[index,] + Y[index2,])/sqrt(2)
    return(Y[ind2, ind2])
}


#####################
#### preparation ####
#####################

# using pivoted cholesky (warning is expected from chol)
getM.chol <- function(tt1, tt2, tol=1e-6){
    K1 <- getK(tt1)
    K2 <- getK(tt2)
    K=K1*K2
    temp <- chol(K, pivot=T, tol=tol)
    r <- attr(temp, "rank")
    oo <- order(attr(temp, "pivot"))
    return(t(temp[1:r, oo]))
}

# sampling a subset of time points
getM.sam <- function(tt1, tt2,nsam){
    K1 <- getK(tt1)
    K2 <- getK(tt2)
    K=K1*K2
    ind <- sample(length(tt), nsam)
    u <- svd(K[,ind])$u
    a <- chol(t(u) %*% K %*% u)
    return(u %*% t(a))
}
# getM<- function(tt1, tt2, tol=1e-6){
#     K1 <- getK(tt1)
#     K2 <- getK(tt2)
#     K=K1*K2
#     Eig=eigen(K)
#     M= Eig$vectors[,1:(q1*q2)]%*%diag(sqrt(Eig$values[1:(q1*q2)]))
#     return(M)
# }

prep2 <- function(time1, time2, x, subject, Mmethod="eig", tol=1e-4, nsam=50, gamma = NULL){
    
    # prepare for fitting
    Xs <- list()
    # tts1 <- list()
    # tts2 <- list()

    if(is.null(gamma)){
        dis2 = get_distance2(time1, time2)
        gamma = median(sqrt(dis2))
    }
    ms <- NULL
    i <- 0
    for (zz in unique(subject)){
        if (sum(subject==zz)>1){
            i <- i+1
            # tts1[[i]] <- time1[subject == zz]
            # tts2[[i]] <- time2[subject == zz]
            Xs[[i]] <- as.double(x[subject == zz])
            ms <- c(ms, length(time1))
        }
    }

    for (zz in unique(subject)){
        if (sum(subject==zz)>1){
            i <- i+1
            # tts1[[i]] <- time1[subject == zz]
            # tts2[[i]] <- time2[subject == zz]
            Xs[[i]] <- as.double(x[subject == zz])
            ms <- c(ms, length(time1))
        }
    }

    print('done!')
    # getting Mi
    n <- length(Xs)
    # note: some element of time may be removed, so use unlist(tts) instead of time
    if (Mmethod=="eig"){
        if(Kernel == 'gauss')
       { M <- getM(unlist(time1),unlist(time2), tol= tol, gamma = gamma)}
       else{
        M <- getM(unlist(time1),unlist(time2), tol= tol)
       }
    } else if (Mmethod=="chol"){
        # M <- getM.chol(unlist(tts1), unlist(tts2),tol=tol)
    } else if (Mmethod=="sam"){
        # M <- getM.sam(unlist(tts1), unlist(tts2),nsam=nsam)
    }
    ii <- c(0, cumsum(ms))
    Ms <- list()
    for (i in (1:n)){
        # Ms[[i]] <- M[(ii[i]+1):ii[i+1],]
        Ms[[i]] <- M
    }
    print(dim(M))
    obj <- Qlossprep_cpp(Xs, Ms, 1:n)
    return(list(Xs=Xs, Ms=Ms, ms=ms, M=M, tts1=time1, tts2=time2, R=obj$R, Qv=obj$Qv, c=obj$c, gamma = gamma))
}

rkhscov.alg.control <- function(Mmethod="eig", preptol=1e-4, nsam=50, L=1,
eta=2, alpha=0.9, maxit=10000, traceit=FALSE,
max.rank=-1, tol=1e-12, variant=2, cond=2, gamma = NULL){
    return(list(Mmethod=Mmethod, preptol=preptol, nsam=nsam, L=L, eta=eta, alpha=alpha, maxit=maxit, traceit=traceit, max.rank=max.rank, tol=tol, variant=variant, cond=cond))
}

rkhscov <- function(time1, time2, x, subject, lam, gam=1, weight=NULL, centered=FALSE, B0v=NULL, pos=TRUE, control=list())
{
    
    control <- do.call(rkhscov.alg.control, control)
    
    pout <- prep2(time1, time2, x, subject, Mmethod=control$Mmethod, tol=control$preptol, nsam=control$nsam)
    r <- ncol(pout$M)
    if (is.null(weight)){
        weight <- rep(1, r)
    } else {
        if (length(weight)<r) weight <- c(rep(1e6*weight[1],r-length(weight)), weight)
        else stop("length of weight does not match!")
        if (any(sort(weight, decreasing=T)!=weight))
        warning("weight should be decreasing to match the order of the eigenvalues!")
    }
    # weight may not work with accelerated gradient descent
    
    # run accelarated proximal gradient method
    if (is.null(B0v)) B0v <- svec(diag(rep(1e-8,r)))
    res <- rkhscov_pg_cpp(RR=pout$R, RQv=pout$Qv, c=pout$c, lam=lam, gam=gam, B0v=B0v,
    weight=weight, L=control$L, eta=control$eta,
    alpha=control$alpha, maxit=control$maxit,
    traceit=as.logical(control$traceit), tol=control$tol,
    max_rank=control$max.rank, pos=as.logical(pos),
    variant=control$variant, cond=control$cond)
    return(list(e=res$e, obj=res$obj, dBs=res$dBs, conv=res$conv, r=r,
    k=length(res$e$values), Ms=pout$Ms, Minv=MASS::ginv(pout$M), t1=pout$tts1, t2=pout$tts2,
    call=match.call(), pos=pos, weight=res$weight, lam=lam, gam=gam))
}

###################
#### k-fold CV ####
###################

gen.groups <- function(n, nfold){
    leave.out <- trunc(n/nfold)
    o <- sample(1:n)
    groups <- vector("list", nfold)
    for (j in (1:(nfold-1))){
        jj <- (1+(j-1)*leave.out)
        groups[[j]]<-(o[jj:(jj+leave.out-1)])
    }
    groups[[nfold]] <- o[(1+(nfold-1)*leave.out):n]
    return(groups=groups)
}

rkhscov.cv.control <- function(lams=NULL, gams=1, lamu=NULL, lam2.rtol=1e-1,
lam.min.ratio=1e-10, nlam=20, ngam=20, fine=TRUE, nfine=20){
    return(list(lams=lams, gams=gams, lamu=lamu, lam2.rtol=lam2.rtol,
    lam.min.ratio=lam.min.ratio, nlam=nlam, ngam=ngam, fine=fine,
    nfine=nfine))
}

rkhscov.cv <- function(time1, time2, x, subject, nfold=5, weight=NULL, centered=FALSE,
ncpu=nfold, pos=TRUE, control.alg=list(),
control.cv=list()){
    control.alg <- do.call(rkhscov.alg.control, control.alg)
    control.cv <- do.call(rkhscov.cv.control, control.cv)
    
    
    pout0 <- prep2(time1, time2, x, subject, Mmethod=control.alg$Mmethod, tol=control.alg$preptol, nsam=control.alg$nsam, gamma = control.alg$gamma)
    r <- ncol(pout0$M)
    n <- length(pout0$Ms)
    
    # initialize weight
    
    if (is.null(weight)){
        weight <- rep(1, r)
    } else {
        if (length(weight)<r) weight <- c(rep(1e6*weight[1],r-length(weight)), weight)
        else stop("length of weight does not match!")
        if (any(sort(weight, decreasing=T)!=weight))
        warning("weight should be decreasing to match the order of the eigenvalues!")
    }
    # weight may not work with accelerated gradient descent
    
    # preparation for CV
    groups <- gen.groups(n, nfold)
    
    # construction of lam sequence
    if (is.null(control.cv$gams)){
        gams <- seq(1e-3, 1, len=control.cv$ngam)
    } else {
        gams <- control.cv$gams
    }
    
    if (is.null(control.cv$lams)){
        if (is.null(control.cv$lamu)){
            if (any(gams<1e-20)){
                if (length(gams)==1)
                lamu <- start_lambda(pout0$R, pout0$Qv, pout0$c, weight, 0, pos,
                control.cv$lam2.rtol) # when we supply gams as c(0)
                else
                stop("gams includes 0, but not of length 1")
            } else {
                lamu <- start_lambda(pout0$R, pout0$Qv, pout0$c, weight, 1, pos,
                control.cv$lam2.rtol) # When gam>0, it gives the same value
            }
        } else {
            lamu <- control.cv$lamu
        }
        laml <- lamu * control.cv$lam.min.ratio
        lams <- exp(seq(log(laml), log(lamu), len=control.cv$nlam)) # on log scale
    } else {
        lams <- control.cv$lams
    }
    
    
    B0v <- svec(diag(rep(1e-8, r)))
    cvs.grid <- matrix(nr=length(lams), nc=length(gams))
    for (k in (1:length(gams))){
        gam <- gams[k]
        rkhscov_pg_cvj=function(j){
            out=rkhscov_pg_cvj_cpp(pout0$Xs, pout0$Ms,  setdiff(1:n, groups[[j]]),
            groups[[j]], lams, gam, B0v, weight, L=control.alg$L,
            eta=control.alg$eta, alpha=control.alg$alpha, maxit=control.alg$maxit,
            traceit=as.logical(FALSE), tol=control.alg$tol,
            max_rank=control.alg$max.rank, pos=as.logical(pos),
            variant=control.alg$variant, cond=control.alg$cond)
            return (out)
        }
        cvjs=lapply(1:nfold,rkhscov_pg_cvj)
        #cvs.se <- apply(sapply(cvjs, function(x){x[,1]}), 1, function(x){sd(x)})/sqrt(nfold)
        cvs <- apply(matrix(sapply(cvjs, function(x){x[,2]}),nc=nfold), 1, sum)/n
        cvs.grid[,k] <- cvs
    }
    cv.best.ind <- arrayInd(which.min(cvs.grid), dim(cvs.grid))
    
    # fine tuning on lambda
    if ((control.cv$fine && (cv.best.ind[1]>1)) && (cv.best.ind[1]<length(lams))){
        B0v <- svec(diag(rep(1e-8, r)))
        lams2 <- seq(lams[cv.best.ind[1]-1], lams[cv.best.ind[1]+1], len=control.cv$nfine)
        gam <- gams[cv.best.ind[2]]
        rkhscov_pg_cvj=function(j){
            out=rkhscov_pg_cvj_cpp(pout0$Xs, pout0$Ms,  setdiff(1:n, groups[[j]]),
            groups[[j]], lams, gam, B0v, weight, L=control.alg$L,
            eta=control.alg$eta, alpha=control.alg$alpha, maxit=control.alg$maxit,
            traceit=as.logical(FALSE), tol=control.alg$tol,
            max_rank=control.alg$max.rank, pos=as.logical(pos),
            variant=control.alg$variant, cond=control.alg$cond)
            return (out)
        }
        cvjs2=lapply(1:nfold,rkhscov_pg_cvj)
        # combine
        lams <- c(lams, lams2); oo <- order(lams); lams <- lams[oo]
        cvs2 <- apply(matrix(sapply(cvjs2, function(x){x[,2]}),nc=nfold), 1, sum)/n
        cvs <- c(cvs.grid[,cv.best.ind[2]], cvs2); cvs <- cvs[oo]
        
        cv.best.ind[1] <- which.min(cvs)
    }
    
    
    B0v <- svec(diag(rep(1e-8, r)))
    res <- rkhscov_pg_cpp(RR=pout0$R, RQv=pout0$Qv, c=pout0$c, lam=lams[cv.best.ind[1]],
    gam=gams[cv.best.ind[2]], B0v=B0v, weight=weight,
    L=control.alg$L, eta=control.alg$eta, alpha=control.alg$alpha,
    maxit=control.alg$maxit*2, traceit=as.logical(control.alg$traceit),
    tol=control.alg$tol, max_rank=control.alg$max.rank,
    pos=as.logical(pos), variant=control.alg$variant,
    cond=control.alg$cond)
    sobj <- list(e=res$e, obj=res$obj, dBs=res$dBs, conv=res$conv,
    r=r, k=length(res$e$values), Ms=pout0$Ms, Minv=MASS::ginv(pout0$M), t1=pout0$tts1,t2=pout0$tts2,
    weight=res$weight, lam=lams[cv.best.ind[1]], gam=gams[cv.best.ind[2]])
    
    # flag information only for lambda
    if (cv.best.ind[1]==1) flag <- 1
    else if (cv.best.ind[1]==length(lams)) flag <- 2
    else flag <- 0
    
    return(list(sobj=sobj, cvs.grid=cvs.grid, cvs=cvs, cv.best.ind=cv.best.ind,
    lams=lams, gams=gams,
    call=match.call(), pos=pos, weight=weight, flag=flag))
}

#####################
#### predictions ####
#####################

# fitted.rkhscov <- function(sobj){
#     out <- list()
#     p <- ncol(sobj$e$vectors)
#     B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
#     for (i in (1:length(sobj$Ms)))
#     out[[i]] <- sobj$Ms[[i]] %*% B %*% t(sobj$Ms[[i]])
#     return(out)
# }

# compute.A <- function(sobj){
#     p <- ncol(sobj$e$vectors)
#     if(p==1){
#         B= as.numeric(sobj$e$values) * sobj$e$vectors %*% t(sobj$e$vectors)
#     }else{
#         B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
#     }
#     A <- t(sobj$Minv) %*% B %*% sobj$Minv
#     return(A)
# }


# predict.rkhscov <- function(newtt1, newtt2,sobj){
#     A <- compute.A(sobj)
#     K1 <- outer(newtt1, unlist(sobj$t1), K.cos)
#     K2 <- outer(newtt2, unlist(sobj$t2), K.cos)
#     K = khatri_rao(K1, K2)
#     return(K %*% A %*% t(K))
# }



###################################
#### functions related to FPCA ####
###################################
library(rTensor)
# using average
# more stable
Qmat3 <- function(ss, ss2,  o=150, Kernel, fix.pos=F, gamma = NULL){
  N <- length(ss)
  tt <- seq(0, 1, len=o)
  if(Kernel == "sob"){
Z1 <- outer(ss, tt, K.sob)
Z2 <- outer(ss2, tt, K.sob)
  }
if(Kernel == "cos"){
Z1 <- outer(ss, tt, K.cos)
Z2 <- outer(ss2, tt, K.cos)
}   
if(Kernel == "gauss"){
Z1 <- outer(ss, tt, K.gauss, gamma)
Z2 <- outer(ss2, tt, K.gauss, gamma)
}   
  Z = t(khatri_rao(t(Z1), t(Z2)))
  out <- Z %*% t(Z) / (o^2)

  if (fix.pos){
    ee <- eigen(out, symmetric=T)
    ee$values[ee$values<0] <- 0
    out <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
  }
  return(out)
}


# using average
# more stable
# Qmat3 <- function(ss, ss2,  o=100, Kernel, fix.pos=F, gamma = NULL){
#   N <- length(ss)
#   tt <- seq(0, 1, len=o)
#   if(Kernel == "sob"){
# Z1 <- outer(ss, sort(ss), K.sob)
# Z2 <- outer(ss2, sort(ss2), K.sob)
#   }
# if(Kernel == "cos"){
# Z1 <- outer(ss, sort(ss), K.cos)
# Z2 <- outer(ss2, sort(ss2), K.cos)
# }   
# if(Kernel == "gauss"){
# Z1 <- outer(ss, sort(ss), K.gauss, gamma)
# Z2 <- outer(ss2, sort(ss2), K.gauss, gamma)
# }   
#   Z = t(khatri_rao(t(Z1), t(Z2)))
# #   Z = Z1 * Z2
#   out <- Z %*% t(Z) / (length(ss) * length(ss2))

#   if (fix.pos){
#     ee <- eigen(out, symmetric=T)
#     ee$values[ee$values<0] <- 0
#     out <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
#   }
#   return(out)
# }



# Qmat3 <- function(ss, ss2,  o=10000, Kernel, fix.pos=F, gamma = NULL){
#   N <- length(ss)
#   tt <- seq(0, 1, len=o)
#   if(Kernel == "sob"){
# Z1 <- outer(ss, ss, K.sob)
# Z2 <- outer(ss2, ss2, K.sob)
#   }
# if(Kernel == "cos"){
# Z1 <- outer(ss, ss, K.cos)
# Z2 <- outer(ss2, ss2, K.cos)
# }   
# if(Kernel == "gauss"){
# Z1 <- outer(ss, ss, K.gauss, gamma)
# Z2 <- outer(ss2, ss2, K.gauss, gamma)
# }   
#   Z = Z1 * Z2
#   out <- Z %*% t(Z) / o

#   if (fix.pos){
#     ee <- eigen(out, symmetric=T)
#     ee$values[ee$values<0] <- 0
#     out <- ee$vectors %*% diag(ee$values) %*% t(ee$vectors)
#   }
#   return(out)
# }


fpca.rkhscov <- function(sobj, Q=NULL, Kernel, fix.pos = F, gamma = gamma){
  if (is.null(Q)){
    Q <- Qmat3(unlist(sobj$t1), unlist(sobj$t2), Kernel = Kernel, fix.pos=fix.pos, gamma = gamma )
  }
  R <- sobj$Minv %*% Q %*% t(sobj$Minv)
  R <- (R + t(R)) / 2
  #Rroot <- matsq_cpp(R, tol=1e-20) # matsq_cpp produce R =Rroot %*% t(Rroot), not Rroot %*% Rroot
  ee <- eigen(R); ee$values[ee$values < 0] <- 0
  Rroot <- ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)
  p <- ncol(sobj$e$vectors)
  B <- sobj$e$vectors %*% diag(sobj$e$values, nr=p, nc=p) %*% t(sobj$e$vectors)
  e <- eigen(Rroot %*% B %*% Rroot, symmetric=T)
  U <- t(sobj$Minv) %*% MASS::ginv(Rroot) %*% e$vectors

  rank <- length(sobj$e$values)
  return(list(t1=sobj$t1,t2=sobj$t2, U=U[,1:rank], values=e$values[1:rank]))
}



compute.fpca <- function(tt1, tt2,  fpca.obj, Kernel, gamma = NULL){
    if(Kernel == "sob"){
  Z1 <- outer(unlist(fpca.obj$t1), tt1, K.sob)
  Z2 <- outer(unlist(fpca.obj$t2), tt2, K.sob)
    }
    if(Kernel == "cos"){
        Z1 <- outer(unlist(fpca.obj$t1), tt1, K.cos)
        Z2 <- outer(unlist(fpca.obj$t2), tt2, K.cos)
    }
    if(Kernel == "gauss"){
        Z1 <- outer(unlist(fpca.obj$t1), tt1, K.gauss, gamma)
        Z2 <- outer(unlist(fpca.obj$t2), tt2, K.gauss, gamma)
    }
Z = Z1 * Z2
  return(t(fpca.obj$U) %*% Z)
}








################## helper functions ##################

# normalize
normalization22 <- function(counts) { 
  varx = apply(counts, 2, var)
  meanx = apply(counts, 2, mean)
  phi = coef(nls(varx ~ meanx + phi * meanx^2, start = list(phi = 1)))
  
  ## regress out log total counts
  norm_counts <- log(counts + 1/(2 * phi))
  total_counts <- apply(counts, 1, sum)
  
  res_norm_counts <- apply(norm_counts, 2, function(x){resid(lm(x ~ log(total_counts)))} )
  
  return(res_norm_counts)
}

filter_count <- function(count, sample_info, min_total = 10,min_percentage = 0.1){
  gene_num <- ncol(count)
  sample_num <- nrow(count)
  if(sum(rowSums(count) < min_total) == 0){
    if (sum(colSums(count == 0) > (1-min_percentage)*sample_num) > 0){
      sample_f <- sample_info
      count_f <- count[,-(which(colSums(count == 0) > (1-min_percentage)*sample_num))]
    }
    else{
      sample_f <- sample_info
      count_f <- count
    }}
  else{
    if (sum(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num) > 0){
      sample_f <- sample_info[-which(rowSums(count)<min_total),]
      count_f <- count[-which(rowSums(count)<min_total),-(which(colSums(count[-which(rowSums(count)<min_total),] == 0) > (1-min_percentage)*sample_num))]
    }
    else{
      sample_f <- sample_info[-which(rowSums(count)<min_total),]
      count_f <- count[-which(rowSums(count)<min_total),]
    }
  }
  return(list(sample_f,count_f))
}
