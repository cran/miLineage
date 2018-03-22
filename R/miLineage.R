###############
# m is the total number of taxa
###############
#library(GUniFrac, lib.loc="/nas02/home/z/t/ztang/Rlibs/")  
# adding functions for restricted permutation: *_restrict
#library("Matrix")  ## use kronecker product from this library to avoid random NaN in the resulting matrix (when running on Unix)
library("geepack")

.F.test <- function(x){
  
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){

  return( min(length(x) * x/rank(x)) )
  
}

# add 06/19/2016
.diag2 <- function(x){
  
  if(length(x)>1){
    return(diag(x))
  }else{
    return(x)
    
  }
  
}

# assume two subjects in each strata
# last updated 05/16/2016
.sampling.strata<- function(var, strata){
  
  n.strata = length(table(strata))
  strata.uniq = unique(strata)
  var.sample = var
  for(i in 1:n.strata){
    index = which(strata==strata.uniq[i])
    if( rbinom(1, 1, 0.5)==1 ){
      tmp = var.sample[index[2]]
      var.sample[index[2]] = var.sample[index[1]]
      var.sample[index[1]] = tmp
    }
    
  }
  return(var.sample)
  
}

########################################
#                                      #
#           One Part Model             #
#                                      #
########################################

## change on 03/28/2016
.Ei.beta <- function(m, p, beta, X.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  for(j in 1:(m-1)){
    
    Ei.out[j] = exp(beta[((j-1)*p+1):(j*p)] %*% X.i)
  
  }
  
 
  Ei.out[m] = 1
    
 
  
  
  return (Ei.out)
}


## change on 03/28/2016
.fun.neg.loglik.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      loglik = loglik + Y[i,] %*% log(P.i)
    }
    
  }
  
  return (-loglik)
  
}

.fun.neg.score.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)

    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      Score.beta = Score.beta + kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      
    }
    
    return (-Score.beta)
  }
  
  
  
}

.fun.score.i.beta <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
#       if(sum.E.i==0){
#         P.i = rep(0,m)
#       }
      
      Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}

#fun.hessian.beta(est.reduce.beta, data.beta, save.list=TRUE)
#beta = est.reduce.beta
#data = data.beta
.fun.hessian.beta <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; X = data$X
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
    
      Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) ) 
     
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }

    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}

 
.Score.test.stat <- function(Y, X, X.par.index){
  
  p = ncol(X)
  
  nY = rowSums(Y)
  
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    X.reduce = X[,-X.par.index, drop=FALSE]
    p.reduce = p - length(X.par.index)    
    par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    
    data.reduce.beta = list(Y=Y, X=X.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta, gr=.fun.neg.score.beta, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y, X=X)
    Score.reduce.beta = .fun.score.i.beta(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    tmp = .fun.hessian.beta(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # change in 04/08/2016(warning! need to change this step in the resampling function too)
    for(i in 1:n){
      B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
    }
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list )   )
  
  
}

# add 05/31/2016: allow for multiple nuisance covariate and multiple covariates of interest
# Score.test.stat.pos.4resampling is the simplier version for two-group comparison
# rename 06/16/2016 Score.test.stat.pos.4Gresampling to Score.test.stat.4Gresampling and use in both one-part model and postive part of the two-part model
.Score.test.stat.4Gresampling <- function(X.perm, X.par.index, S.beta.list, I.beta.list){
  
  n = nrow(X.perm)  
  p = ncol(X.perm)
  
  m.beta = length(S.beta.list[[1]])
  
  #n.par.interest.beta = m.beta
  n.beta = m.beta*p
  #n.beta = m.beta*2
  par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  #par.interest.index.beta = (1:m.beta )*2
  n.par.interest.beta = length(par.interest.index.beta) 

  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  for(i in 1:n){
    B2 = B2 + Score.reduce.beta.perm.reorg[i,] %o% Score.reduce.beta.perm.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.beta.perm = A %*% ginv(B) %*% A
  
  
  return(score.stat.beta.perm)
  
  
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.one <- function(X, X.par.index, score.stat.beta, S.beta.list, I.beta.list, start.nperm, end.nperm, n.one, one.acc){
  
  n = nrow(X)
  
  n.one.new = n.one
  one.acc.new = one.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    
    score.stat.beta.perm = try( .Score.test.stat.4Gresampling(X.perm, X.par.index, S.beta.list, I.beta.list) )
    
    if(class(score.stat.beta.perm) != "try-error"){
      
      n.one.new = n.one.new + 1
      if(score.stat.beta.perm >= score.stat.beta){
        one.acc.new = one.acc.new + 1
        
      }  
    }
  }
  
  if(one.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(one.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
#   else if(one.acc.new<20){
#     next.end.nperm = ( end.nperm + 1) * 5 - 1;
#     flag = 1;
#     
#   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.one.new=n.one.new, one.acc.new=one.acc.new, flag=flag, next.end.nperm=next.end.nperm))
  
}


# Y: nxm count of microbiomes
# X: covariates
# Z: covariates 
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part
# change 06/16/2016, change the resampling part, to accommodate multiple covariate and potential confounders

.Score.test <- function(Y, X, X.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  
 
  p = ncol(X)
  
  nY = rowSums(Y)
  
  ## remove 03/28/2016
#   nY0.index = which(nY==0)
#   if(length(nY0.index)>0){
#     Y = Y[-nY0.index, , drop=FALSE]
#     X = X[-nY0.index, , drop=FALSE]
#   }
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NULL
    score.pvalue.beta = NA
    n.par.interest.beta = 0
    
  }else{
    
    tmp.one = try( .Score.test.stat(Y, X, X.par.index) )
    if(class(tmp.one) == "try-error"){
      
      score.stat.beta = NA
      score.pvalue.beta = NA
      n.par.interest.beta = NA
      
    }else{
      
      n.par.interest.beta = (m-1)*length(X.par.index)
      score.stat.beta = tmp.one$score.stat.beta
      score.pvalue.beta = 1 - pchisq(score.stat.beta,  n.par.interest.beta) 
     
    }
    

  }
  
 
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue.beta, df =n.par.interest.beta)
#  print(score.stat.beta)
    
  if(resample){
   ########################################
   #                                      #
   #        Resampling Score Test         #
   # change 07/02/2016 for adaptive resampling
   #                                      #
   ########################################  
   set.seed(seed)
   if(!is.na(score.stat.beta)){
     
     n.one = 0
     one.acc = 0
     
     start.nperm = 1;
     end.nperm = min(100,n.replicates);
     flag = 1
     while(flag & end.nperm <= n.replicates){
       
       results = .resample.work.one(X, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, end.nperm, n.one, one.acc)
       n.one = results$n.one.new
       one.acc = results$one.acc.new
       flag = results$flag
       next.end.nperm = results$next.end.nperm
       
       if(flag){
         start.nperm = end.nperm + 1;
         end.nperm = next.end.nperm;
         
       }
       
       if(start.nperm < n.replicates & end.nperm > n.replicates){ 
         warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
         results = .resample.work.one(X, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, n.replicates, n.one, one.acc)
         n.one = results$n.one.new
         one.acc = results$one.acc.new
         
       }
       
     }
     

#      print(paste("Final number of resamplings: ", n.one) )
#      if(n.one<end.nperm/2){
#        print("Number of resamplings too small for one-part test")
#      }
     
     tmp = (one.acc+1)/(n.one+1)
     
     #print(n.one)
     #print(one.acc)
     
   }else{
     
     tmp = NA
   }
   
   beta.results = c(beta.results, score.Rpvalue = tmp)
   
   
  }

  return(beta.results)
  
  
}
# Y: nxm count of microbiomes
# case (case/control status: 1 for cases, 0 for controls)
# score.stat.gamma: score statistics for gamma part
# score.stat.beta: score statistics for beta part
# S.gamma.list, S.beta.list: list with length n, the element i is the score statistics (length = m-1) of the subject i 
# I.gamma.list, I.beta.list: list with length n, the element i is the Hessian matrix ( (m-1)x(m-1) ) of the subject i 

.Tstat <- function(Y, case){
  
  m = ncol(Y)
  n = nrow(Y)
  n1 = sum(case==1)
  
  Ym = Y[,-m, drop=FALSE]
  nY = rowSums(Y)
  nY1 = sum(nY[case==1])
  nY2 = sum(nY[case==0])
  p0 = colSums(Ym)/sum(nY)
  
  A = colSums(Ym[case==1, , drop=FALSE]) - p0 * sum(nY[case==1])
  
  B = matrix(0, m-1, m-1)
  
  C = matrix(0, m-1, m-1)
  index.case = which(case==1) 
  n.case = length(index.case)
  for(i in 1:n.case){
    
    tmp = Ym[index.case[i], ] - p0 * nY[index.case[i]]
    C.tmp = tmp %o% tmp
    C = C + C.tmp
    
  }
  
  D = matrix(0, m-1, m-1)
  for(i in 1:n){
    
    tmp = Ym[i, ] - p0 * nY[i]
    D.tmp = tmp %o% tmp
    D = D + D.tmp
    
  }
  
  B = ((nY2-nY1)/sum(nY))*C + ((nY1/sum(nY))^2)*D
  
  #B = ((nY2-nY1)/sum(nY))*C + ((nY1/sum(nY))^2)*D - ((n1*n1)/n)*((A/n1)%o%(A/n1))
  
  #   B1 = cbind(diag(m-1), -(sum(Y[1:25,])/sum(Y))*diag(m-1) )
  # 
  #   tmp = matrix(NA, nrow=50, ncol=m-1)
  #   for(i in 1:50){
  #     tmp[i,] = Ym[i,] - p0 * nY[i]
  #     
  #   }
  #   #print(colSums(tmp))
  #   #S1 = colMeans(tmp[1:25,])
  #   #print("S1"); print(S1)
  # 
  #   Score.reduce.reorg = rbind(cbind(tmp[1:25,], tmp[1:25,]), cbind(matrix(0, nrow=25, ncol=m-1), tmp[26:50,]))
  # 
  #   Score.reduce.mean = colMeans(Score.reduce.reorg)
  #   print("Score.reduce.mean"); print(Score.reduce.mean)
  ## change 04/14/2016
  #   Score.reduce.center = Score.reduce.reorg   
  #     for(i in 1:n){
  #       Score.reduce.center[i,] = Score.reduce.center[i,] - Score.reduce.mean
  #       
  #     }
  #   B2 =  matrix(0, (m-1)*2, (m-1)*2)
  #   for(i in 1:n){
  #     B2 = B2 + Score.reduce.center[i,] %o% Score.reduce.center[i,] 
  #   }
  #print(nrow(Score.reduce.reorg))
  #print(ncol(Score.reduce.reorg))
  #   D = matrix(0, nrow=m-1, ncol=m-1)
  #   for(i in 1:25){
  #     D  = D + tmp[i,] %o% tmp[i,]
  #   }
  # 
  #   tmp.mean = colMeans(tmp)
  #   print("tmp.mean"); print(tmp.mean)
  #   E = matrix(0, nrow=m-1, ncol=m-1)
  #   for(i in 1:50){
  #     E  = E + (tmp[i,]-tmp.mean) %o% (tmp[i,]-tmp.mean)
  #   }
  #   print("E"); print(E)
  # 
  #   S1 = colMeans(tmp[1:25,])
  #   tmp = tmp - (25/50)*S1
  #   C = matrix(0, nrow=m-1, ncol=m-1)
  #   for(i in 1:25){
  #     C  = C + tmp[i,] %o% tmp[i,]
  #   }
  #   C = C + ((25*25*25)/50*50)*(S1 %o% S1)
  # 
  #   B2 = rbind(cbind( C, D ), cbind( D, E))
  
  # change in 04/08/2016(warning! need to change this step in the resampling function too)
  #     for(i in 1:n){
  #       B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
  #     }
  
  #  B = B1 %*% B2 %*% t(B1)
  
  
  score.stat = A %*% ginv(B) %*% A
  
  #print("A"); print(A)
  #print("B1"); print(B1)
  #print("B2"); print(B2)
  #print("B"); print(B)
  #print("simulated stat:")
  #print(score.stat.beta)
  
  return(score.stat)
}

.Score.test.simple <- function(Y, case, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  #print(A)
  #print(B)
  m = ncol(Y)
  
  score.stat.beta = try( .Tstat(Y, case) )
  if(class(score.stat.beta) == "try-error"){
    
    score.stat.beta = NA
    score.pvalue = NA
    df = NA
    
  }else{
    
    score.pvalue = 1 - pchisq(score.stat.beta, m-1) 
    df = m-1
  }
  
  
  #print("observed score statistics:")
  #print(score.stat.beta)
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue, df = df )
  
  if(resample){
    
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      #print("simulated stat:")
      #set.seed(16)
      
      for(k in 1:n.replicates){
        
        case.perm = sample(case)
        score.stat.beta.perm = try( .Tstat(Y, case.perm) )
        
        #       if(k<10){
        #         print("case.perm"); print(case.perm)
        #         #       print("A"); print(A)
        #         #       print("Hess.reduce.beta.perm.reorg"); print(Hess.reduce.beta.perm.reorg)
        #         #       print("B1"); print(B1)
        #         #       print("B2"); print(B2)
        #         #       print("B"); print(B)
        #         print(score.stat.beta.perm)
        #         
        #       }
        if(class(score.stat.beta.perm) != "try-error"){
          
          n.one = n.one + 1
          if(score.stat.beta.perm >= score.stat.beta){
            one.acc = one.acc + 1
            
          }           
          
        }
        
        
      }
      
      if(n.one<n.replicates/2){
        print("#replicate too small for one-part test")
      }
      
      score.Rpvalue = one.acc/n.one
      
      
    }else{
      
      score.Rpvalue = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = score.Rpvalue)
    
  }
  
  
  return(beta.results)
  
}


.Score.test.simple.restrict <- function(Y, case, strata, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  #print(A)
  #print(B)
  m = ncol(Y)
  
  score.stat.beta = try( .Tstat(Y, case) )
  if(class(score.stat.beta) == "try-error"){
    
    score.stat.beta = NA
    score.pvalue = NA
    df = NA
    
  }else{
    
    score.pvalue = 1 - pchisq(score.stat.beta, m-1) 
    df = m-1
  }
  
  
  #print("observed score statistics:")
  #print(score.stat.beta)
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue, df = df )
  
  if(resample){
    
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      #print("simulated stat:")
      #set.seed(16)
      
      for(k in 1:n.replicates){
        
        case.perm = .sampling.strata(case, strata)
        
        score.stat.beta.perm = try( .Tstat(Y, case.perm) )
        
        
        
        #       if(k<10){
        #         print("case.perm"); print(case.perm)
        #         #       print("A"); print(A)
        #         #       print("Hess.reduce.beta.perm.reorg"); print(Hess.reduce.beta.perm.reorg)
        #         #       print("B1"); print(B1)
        #         #       print("B2"); print(B2)
        #         #       print("B"); print(B)
        #         print(score.stat.beta.perm)
        #         
        #       }
        if(class(score.stat.beta.perm) != "try-error"){
          
          n.one = n.one + 1
          if(score.stat.beta.perm >= score.stat.beta){
            one.acc = one.acc + 1
            
          }           
          
        }
        
        
      }
      
      if(n.one<n.replicates/2){
        print("#replicate too small for one-part test")
      }
      
      score.Rpvalue = one.acc/n.one
      
      
    }else{
      
      score.Rpvalue = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = score.Rpvalue)
    
  }
  
  
  return(beta.results)
  
}

########################################
#                                      #
#    Two Part Model: positive Part     #
#    ## add on 04/25/2016              #
#                                      #
########################################


.Ei.beta.pos <- function(m, p, beta, X.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  I.i = as.numeric(Y.i>0)
  
  for(j in 1:(m-1)){
    
    Ei.out[j] = I.i[j]*exp(beta[((j-1)*p+1):(j*p)] %*% X.i)
    
  }
  
  
  Ei.out[m] = I.i[m]
  
  
  return (Ei.out)
}

.fun.neg.loglik.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      index = which(Y[i,]>0)
      loglik = loglik + Y[i,index] %*% log(P.i[index])
#       if(is.na(loglik)){
#         print(i); break;
#       }
    }
    
  }
  
  return (-loglik)
  
}

.fun.neg.score.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      Score.beta = Score.beta + kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      
    }
    
    return (-Score.beta)
  }
  
  
  
}

.fun.score.i.beta.pos <- function(beta, data){
  
  Y = data$Y; X = data$X; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  
  n.beta = (m - 1)*p
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      Score.beta.i[i,] =  kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1)) 
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}

.fun.hessian.beta.pos <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; X = data$X
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(X)
  n.beta = (m-1)*p
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos(m, p, beta, X[i,], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      Hessian.beta = Hessian.beta + kronecker(tmp.beta, ( X[i,] %o% X[i,] ))
   
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


.Score.test.stat.pos <- function(Y, X, X.par.index){

    p = ncol(X)
    
    nY = rowSums(Y)
    
    n = nrow(Y)
    m = ncol(Y)  
    n.beta = (m - 1)*p
    
    if(sum(X.par.index == 1)){
      stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
    }  
    
    if(is.null(X.par.index) || n==0){
      
      score.stat.beta = NA
      
    }else{
      
      X.reduce = X[,-X.par.index, drop=FALSE]
      p.reduce = p - length(X.par.index)    
      par.interest.index.beta =  kronecker(((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
      
      n.par.interest.beta = length(par.interest.index.beta) 
      
      beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
      
      data.reduce.beta = list(Y=Y, X=X.reduce)
      
      est.reduce.beta = rep(NA, n.beta)
      est.reduce.beta[par.interest.index.beta] = 0 
      #est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta)$par
      # change 04/08/2016
      est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta, method="BFGS")$par
      
      
      data.beta = list(Y=Y, X=X)
      Score.reduce.beta = .fun.score.i.beta.pos(est.reduce.beta, data.beta)
      
      # for resampling: S.beta.list, I.beta.list
      S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
      tmp = .fun.hessian.beta.pos(est.reduce.beta, data.beta, save.list=TRUE)
      I.beta.list = tmp$I.beta.list  
      
      Hess.reduce.beta = tmp$Hessian.beta
      #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
      
      # re-organized the score statistics and Hessian matrix
      Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
      Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
      
      
      A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
      
      B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
      
      
      B2 =  matrix(0, n.beta, n.beta)
      
      # change in 04/08/2016(warning! need to change this step in the resampling function too)
      for(i in 1:n){
        B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
      }
      
      B = B1 %*% B2 %*% t(B1)
      score.stat.beta = A %*% ginv(B) %*% A
      
      
    }

  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list )   )
  
  
}

.Score.test.stat.pos.4resampling <- function(case.perm, S.beta.list, I.beta.list){
  
  n = length(case.perm)  
  
  m.beta = length(S.beta.list[[1]])
  n.par.interest.beta = m.beta
  n.beta = m.beta*2
  par.interest.index.beta = (1:m.beta )*2
  beta.acc = 0
  
  #print("simulated stat:")
  #set.seed(16)
  
  
  XZ = cbind(1,case.perm)
  
  
  Score.reduce.beta.perm = matrix(0, n, m.beta*2 )
  Hess.reduce.beta.perm = matrix(0, m.beta*2, m.beta*2 )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(XZ[i,], ncol=1))  
    
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  XZ[i,] %o% XZ[i,] ) )
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  for(i in 1:n){
    B2 = B2 + Score.reduce.beta.perm.reorg[i,] %o% Score.reduce.beta.perm.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.beta.perm = A %*% ginv(B) %*% A
  
  
  return(score.stat.beta.perm)
  
  
  
}


########################################
#                                      #
#          Two Part Model              #
#            (zero part)               #
#                                      #
########################################
.Tstat.zero <- function(Y0, case){
  
  m = ncol(Y0)
  n = nrow(Y0)
  n1 = sum(case==1)  
  n2 = n-n1
   
  A = (colSums(Y0[case==1, ,drop=FALSE])/n1  - colSums(Y0[case==0, ,drop=FALSE])/n2)
  
  BB = matrix(0, m, m)
  for(i in 1:n){
    BB = BB + Y0[i,] %o% Y0[i,]
    
  }
  
  tmp = colSums(Y0)/n
  B = (1/n1 + 1/n2)*(BB/n - tmp%o%tmp)
  
  score.stat = A %*% ginv(B) %*% A
  return(score.stat)
  
}


# add 05/31/2016: function GEE.zero to run general zero-part wald test
#Y0 = Y
#Y0[Y==0] = 1  
#Y0[Y>0] = 0
#remove.index = which(colSums(Y0)==0)
#Y0 = Y0[,-remove.index]
#Z = cbind( 1, c(rep(0, 25), rep(1, 25)), rnorm(50) )
#Z.par.index = 2
#GEE.zero(Y0, Z, Z.par.index, "independence")
.GEE.zero <- function(Y0, Z, Z.par.index, cor.stru){
  
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova = NULL
  cova.reduce = NULL
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p
    
    index.start.reduce = 1
    index.end.reduce = p.reduce
    
    for(j in 1:m){
      tmp = rep(0, m*p)
      tmp[index.start:index.end] = Z[i,]
      cova = rbind(cova, tmp )
      index.start = index.start + p
      index.end = index.end + p
      
      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp )
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  
  gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence", std.err ="jack")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence", std.err ="jack")
  wald.test = anova(gee.full, gee.reduce)

  
  return(wald.test)
  
}


.Pi.alpha<-function(m, p, alpha, X.i){
  
  Pi.out = rep(NA,m)
  
  for(j in 1:m){
    
    tmp = exp(alpha[((j-1)*p+1):(j*p)] %*% X.i)
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }
    
  }
  
  
  return (Pi.out)    
}

#fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
#alpha = est.reduce.alpha
#data = data.alpha 
.fun.score.i.alpha <- function(alpha, data, save.list=FALSE){
  
  Y = data$Y; Z = data$Z; 
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  
  vA.list = list()
  Vinv.list = list()
  VY.list = list()
  
  n.alpha = m*p
  
  if(length(alpha)!=n.alpha){
    
    warning("Dim of initial alpha does not match the dim of covariates")
    
  }else{
    
    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(vA.tmp) 
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure
      
      tmp.V.i = ginv(V.i)
      tmp.VY = tmp.V.i %*% (Y[i,] - Pi.i)
      Score.alpha.i[i,] = t.D.i %*% tmp.VY
      
      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }
    
    
  }
  
  if(save.list){
    
    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )
    
  }else{
    
    return (Score.alpha.i)
  }  
  
}

#fun.hessian.alpha(est.reduce.alpha, data.alpha)
.fun.hessian.alpha <- function(alpha, data){
  
  Y = data$Y; Z = data$Z
  
  n = nrow(Y)
  m = ncol(Y)
  p = ncol(Z)
  n.alpha = m*p
  
  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")
    
  }else{
    
    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)
    
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha(m, p, alpha, Z[i,])
      tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(tmp)
      t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      V.i = A.i # independent cor structure        
      
      
      Hessian.alpha = Hessian.alpha + t.D.i %*% ginv(V.i) %*% t(t.D.i) 
      
      
    }
    
    return (Hessian.alpha)
    
    
  }
  
  
}

# add 06/01/2016: function to run general GEE score test for the zero-part
#Y0 = Y
#Y0[Y==0] = 1  
#Y0[Y>0] = 0
#remove.index = which(colSums(Y0)==0)
#Y0 = Y0[,-remove.index]

#Z.par.index = 2
#GEE.zero(Y0, Z, Z.par.index, "independence")

.Score.test.stat.zero <- function(Y0, Z, Z.par.index, cor.stru){
  
  Z.reduce = Z[,-Z.par.index,drop=FALSE]
  n = nrow(Y0)
  m = ncol(Y0)
  p = ncol(Z)
  p.reduce = ncol(Z.reduce)
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    index.start = 1
    index.end = p
    
    index.start.reduce = 1
    index.end.reduce = p.reduce
    
    for(j in 1:m){
      
      tmp = rep(0, m*p.reduce)
      tmp[index.start.reduce:index.end.reduce] = Z.reduce[i,]
      cova.reduce = rbind(cova.reduce, tmp )
      index.start.reduce = index.start.reduce + p.reduce
      index.end.reduce = index.end.reduce + p.reduce
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  #data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  #gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)
  
  
  ########### perform score test
  n.alpha = m *p
  par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha) 
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce
  
  data.alpha = list(Y=Y0, Z=Z)

  tmp = .fun.score.i.alpha(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list
    
  Hess.reduce.alpha =  .fun.hessian.alpha(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha) 
  
   
  
  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list )   )
  
}

.Score.test.stat.zero.4Gresampling <- function(Z.perm, Z.par.index, vA.list, Vinv.list, VY.list){
  
  n = nrow(Z.perm)
  p = ncol(Z.perm)
  m.alpha = length(vA.list[[1]])
  n.alpha = m.alpha*p
    
  par.interest.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  n.par.interest.alpha = length(par.interest.index.alpha) 
  
  Score.reduce.alpha.perm = matrix(0, n, n.alpha )
  Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         alpha part: resampling Score test        #
    #                                                 #
    ################################################### 
    tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))
    
    Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + tD.tmp %*% VY.list[[i]]
    
    Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)

    
  }
  
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha.perm[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  for(i in 1:n){
    B2 = B2 + Score.reduce.reorg[i,] %o% Score.reduce.reorg[i,] 
  }
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  
  
  return(score.stat.alpha)
  
}

########################################
#                                      #
#          Two Part Model              #
#                                      #
#                                      #
########################################

# add 05/31/2016: function Score.test2.zerowald to run two-part test in the general setting (zero part: use wald GEE test)
# Y: nxm count of microbiomes
# X: covariates for positive part: first column is always intercept
# Z: covariates for zero part: first column is always intercept
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part
.Score.test2.zerowald <- function(Y, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  n = nrow(X)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      index.cova = 1 + which(apply(X1[,-1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      X1.par.index.ava = as.numeric( !is.na(match(X.par.index, index.cova))) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( length(index.cova)<1 |  sum(X1.par.index.ava)==0 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        X1 = X1[,c(1, index.cova), drop=FALSE] 
        tmp = match(X.par.index, index.cova) + 1
        X1.par.index = tmp[!is.na(tmp)]
        d1.pos = length(X1.par.index)
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, X1.par.index) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos*(m-1) ) 
          df.pos = d1.pos * (m-1)
        }
        
      }    
      
    }
    

    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .GEE.zero(Y0, Z, Z.par.index, "independence") )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = score.stat.zero[[3]]
        df.zero = score.stat.zero[[1]]
        score.stat.zero = score.stat.zero[[2]]
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      df.comb = df.zero + df.pos
      score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = df.pos
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = df.zero
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      #set.seed(16)
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          perm.index = sample(1:n)
          X.perm = X
          X.perm[,X.par.index,drop=FALSE] = X.perm[perm.index,X.par.index,drop=FALSE]
          X1.perm = X.perm[index.subj.pos, , drop=FALSE]
          X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
          
          Z.perm = Z
          Z.perm[,Z.par.index,drop=FALSE] = Z.perm[perm.index,Z.par.index,drop=FALSE]
          
          score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .GEE.zero(Y0, Z.perm, Z.par.index, "independence")[[2]] )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          #case.perm = sample(case)
          #case1.perm = case.perm[index.subj.pos]
          perm.index = sample(1:n)
          X.perm = X
          X.perm[,X.par.index,drop=FALSE] = X.perm[perm.index,X.par.index,drop=FALSE]
          X1.perm = X.perm[index.subj.pos, , drop=FALSE]
          X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
          
          score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          #case.perm = sample(case)
          perm.index = sample(1:n)
          Z.perm = Z
          Z.perm[,Z.par.index,drop=FALSE] = Z.perm[perm.index,Z.par.index,drop=FALSE]
          
          score.stat.zero.perm = try( .GEE.zero(Y0, Z.perm, Z.par.index, "independence")[[2]] )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}



# add 07/02/2016 for adaptive resampling
.resample.work.two <- function(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, pos.S.beta.list, pos.I.beta.list, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc){
  
  n = nrow(X)
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  n.comb.new = n.comb
  comb.acc.new = comb.acc
  
  score.stat.comb = score.stat.pos + score.stat.zero
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    X1.perm = X.perm[index.subj.pos, , drop=FALSE]
    X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
    
    Z.perm = Z
    Z.perm[,Z.par.index] = Z.perm[perm.index,Z.par.index]
    
    score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, pos.S.beta.list,pos.I.beta.list) )
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling(Z.perm, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    if(class(score.stat.pos.perm) != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    if(class(score.stat.zero.perm) != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
      
      score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
      n.comb.new = n.comb.new + 1
      if(score.stat.comb.perm >= score.stat.comb){
        comb.acc.new = comb.acc.new + 1
        
      }  
    }
    
  }
  
  if(comb.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(comb.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              n.comb.new=n.comb.new, comb.acc.new=comb.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 07/02/2016 for adaptive resampling
.resample.work.pos <- function(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, pos.S.beta.list, pos.I.beta.list, start.nperm, end.nperm, n.pos, pos.acc){
  
  n = nrow(X)
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    X.perm = X
    X.perm[,X.par.index] = X.perm[perm.index,X.par.index]
    X1.perm = X.perm[index.subj.pos, , drop=FALSE]
    X1.perm = X1.perm[,c(1, index.cova), drop=FALSE] 
    
    
    score.stat.pos.perm = try( .Score.test.stat.4Gresampling(X1.perm, X1.par.index, pos.S.beta.list,pos.I.beta.list) )
    
    
    if(class(score.stat.pos.perm) != "try-error"){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    
    
  }
  
  if(pos.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(pos.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 07/02/2016 for adaptive resampling
.resample.work.zero <- function(Z, Z.par.index, score.stat.zero, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.zero, zero.acc){
  
  n = nrow(Z)
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  for(k in start.nperm:end.nperm){
    
    perm.index = sample(1:n)
    
    Z.perm = Z
    Z.perm[,Z.par.index] = Z.perm[perm.index,Z.par.index]
    
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling(Z.perm, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    
    if(class(score.stat.zero.perm) != "try-error"){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    
  }
  
  if(zero.acc.new < 1){
    next.end.nperm = (end.nperm + 1) * 100 - 1;
    flag = 1;
    
  }else if(zero.acc.new<10){
    next.end.nperm = ( end.nperm + 1) * 10 - 1;
    flag = 1;
    
  }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = ( end.nperm + 1) - 1;
    flag = 0;  
  }
  
  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}

# add 06/01/2016: function Score.test2 to run two-part test in the general setting(zero part: use score GEE test)
# Y: nxm count of microbiomes
# X: covariates for positive part: first column is always intercept
# Z: covariates for zero part: first column is always intercept
# X.par.index: index for the parameter of interest for the X part
# Z.par.index: index for the parameter of interest for the Z part

# Score.test2(Y.rff, X, 2, X, 2, seed=11, resample=TRUE, n.replicates=1000)   
# Y = Y.rff
# X.par.index = 2
# Z = X
# Z.par.index = 2
.Score.test2 <- function(Y, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL){
  

  n = nrow(X)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      index.cova = 1 + which(apply(X1[,-1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      X1.par.index.ava = as.numeric( !is.na(match(X.par.index, index.cova))) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( length(index.cova)<1 |  sum(X1.par.index.ava)==0 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        X1 = X1[,c(1, index.cova), drop=FALSE] 
        tmp = match(X.par.index, index.cova) + 1
        X1.par.index = tmp[!is.na(tmp)]
        d1.pos = length(X1.par.index)
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, X1.par.index) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos*(m-1) ) 
          df.pos = d1.pos * (m-1)
        }
        
      }    
      
    }
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      
      tmp.zero = try( .Score.test.stat.zero(Y0, Z, Z.par.index, "independence") )
      if(class(tmp.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        df.zero = tmp.zero[[1]]
        score.stat.zero = tmp.zero[[2]]
        score.pvalue.zero = tmp.zero[[3]]
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      df.comb = df.zero + df.pos
      score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = df.pos
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = df.zero
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
          
          results = .resample.work.two(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc)
            
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          n.comb = results$n.comb.new
          comb.acc = results$comb.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            results = .resample.work.two(X, X.par.index, X1.par.index, Z, Z.par.index, index.subj.pos, index.cova, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            n.comb = results$n.comb.new
            comb.acc = results$comb.acc.new
            
          }
          
        }
        
#         if(n.pos<n.replicates/2){
#           print("#replicate too small for pos test")
#         }
#         if(n.zero<n.replicates/2){
#           print("#replicate too small for zero test")
#         }
#         if(n.comb<n.replicates/2){
#           print("#replicate too small for comb test")
#         }
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        score.Rpvalue.comb = (comb.acc+1)/(n.comb+1)
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
           
          results = .resample.work.pos(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, end.nperm, n.pos, pos.acc)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings")) 
            results = .resample.work.pos(X, X.par.index, X1.par.index, index.subj.pos, index.cova, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, n.replicates, n.pos, pos.acc)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            
          }
          
        }
        
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        while(flag & end.nperm <= n.replicates){
          
            
          results = .resample.work.zero(Z, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.zero, zero.acc)
          
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          if(flag){
            start.nperm = end.nperm + 1;
            end.nperm = next.end.nperm;
            
          }
          
          if(start.nperm < n.replicates & end.nperm > n.replicates){ 
            warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            
            results = .resample.work.zero(Z, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.zero, zero.acc)
            
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            
          }
          
        }
                
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

.Score.test.simple2 <- function(Y, case, seed=11, resample=FALSE, n.replicates=NULL){
  
  set.seed(seed)
  
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    X = cbind(1, case)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects are all case/control
      if( length(table(X1[,2]))<=1 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, 2) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, m-1) 
          df.pos = m-1
        }
        
      }
      
      
      
    }
    
    
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .Tstat.zero(Y0, case) )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = 1 - pchisq(score.stat.zero, m0) 
        df.zero = m0
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      score.pvalue.comb = 1 - pchisq(score.stat.comb, (m0 + m -1) ) 
      df.comb = m0 + m -1
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = m - 1
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = m0
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      #set.seed(16)
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = sample(case)
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
    

  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

.Score.test.simple2.restrict <- function(Y, case, strata, seed=11, resample=FALSE, n.replicates=NULL){
  
  
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      comb.results = c(comb.results, score.Rpvalue = NA)
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE] 
    }
    
    m = ncol(Y)
    X = cbind(1, case)
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      X1 = X[index.subj.pos, , drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects are all case/control
      if( length(table(X1[,2]))<=1 | m<=1){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        tmp.pos = try( .Score.test.stat.pos(Y1, X1, 2) )
        if(class(tmp.pos) == "try-error"){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          
          score.stat.pos = tmp.pos$score.stat.beta
          score.pvalue.pos = 1 - pchisq(score.stat.pos, m-1) 
          df.pos = m-1
        }
        
      }
      
      
      
    }
    
    
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
      }
      m0 = ncol(Y0)
      
      score.stat.zero = try( .Tstat.zero(Y0, case) )
      if(class(score.stat.zero) == "try-error"){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        score.pvalue.zero = 1 - pchisq(score.stat.zero, m0) 
        df.zero = m0
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos + score.stat.zero
      score.pvalue.comb = 1 - pchisq(score.stat.comb, (m0 + m -1) ) 
      df.comb = m0 + m -1
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.pos
      score.pvalue.comb = score.pvalue.pos
      df.comb = m - 1
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      score.stat.comb = score.stat.zero
      score.pvalue.comb = score.pvalue.zero  
      df.comb = m0
      
    }else{
      
      score.stat.comb = NA
      score.pvalue.comb =  NA
      df.comb = NA
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos, df = df.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero, df = df.zero)
    comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb, df = df.comb )
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            } 
          }
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }  
          }
          
          if(class(score.stat.pos.perm) != "try-error" & class(score.stat.zero.perm) != "try-error"){
            
            score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
            n.comb = n.comb + 1
            if(score.stat.comb.perm >= score.stat.comb){
              comb.acc = comb.acc + 1
              
            }  
          }
          
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        if(n.comb<n.replicates/2){
          print("#replicate too small for comb test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        score.Rpvalue.zero = zero.acc/n.zero
        score.Rpvalue.comb = comb.acc/n.comb
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          case1.perm = case.perm[index.subj.pos]
          
          score.stat.pos.perm = try( .Score.test.stat.pos.4resampling(case1.perm, tmp.pos$S.beta.list, tmp.pos$I.beta.list) )
          
          if(class(score.stat.pos.perm) != "try-error"){
            
            n.pos = n.pos + 1
            if(score.stat.pos.perm >= score.stat.pos){
              pos.acc = pos.acc + 1
              
            }  
          }
        }
        
        if(n.pos<n.replicates/2){
          print("#replicate too small for pos test")
        }
        
        score.Rpvalue.pos = pos.acc/n.pos
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        for(k in 1:n.replicates){
          
          case.perm = .sampling.strata(case, strata)
          score.stat.zero.perm = try( .Tstat.zero(Y0, case.perm) )
          
          if(class(score.stat.zero.perm) != "try-error"){
            
            n.zero = n.zero + 1
            if(score.stat.zero.perm >= score.stat.zero){
              zero.acc = zero.acc + 1
              
            }          
          }
          
          
          
        }
        
        if(n.zero<n.replicates/2){
          print("#replicate too small for zero test")
        }
        
        score.Rpvalue.zero = zero.acc/n.zero
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}

# one-part test
# no missing value is allowed in any of the inputs
#
# OTU: a matrix contains counts with each row corresponds to a sample and each column corresponds to an OTU or a taxa. Column name is mandatory
#
# Tax: a matrix define the taxonomy ranks with each row corresponds to an OTU or a taxa and each column corresponds to a rank (start from the higher taxonomic level). Row name is mandatory and should be consistent with the column name of the OTU table,  Column name should be formated as "Rank1", "Rank2 ()"... etc 
#       If provided, tests will be performed for lineages based on the taxonomic rank. The output contains P-values for all lineages; a list of significant lineages controlling the false discovery rate (based on resampling p-value if resampling test was performed); p-values of the global tests (Fisher- and Simes-combined the p-values for testing lineages).
#       If not provided, one test will be performed with all the OTUs and one p-value will be output
#
# X: a matrix contains covariates for the positive-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders)
#
# X.index: a vector indicate the columns in X for the covariate(s) of interest
#
# Z: a matrix contains covariates for the zero-part test with each column pertains to one variable (pertains to the covariate of interest or the confounders)
#
# Z.index: a vector indicate the columns in X for the covariate(s) of interest
#
# min.depth: keep samples with depths >= min.depth
#
# n.resample: perform asymptotic test is n.resample is null, other perform resampling tests using the specified number of resamplings.
#
# fdr.alpha: false discovery rate for multiple tests on the lineages.
 
QCAT <- function(OTU, X, X.index, Tax=NULL, min.depth=0, n.resample=NULL, fdr.alpha=0.05){
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep, drop=FALSE]
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      pval = as.matrix( .Score.test(count, X, X.index, resample=FALSE, n.replicates=NULL)$score.pvalue )
      colnames(pval) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test(count, X, X.index, resample=TRUE, n.replicates=n.resample)
      pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    
    tax = Tax[keep, ,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              pval = cbind(pval, .Score.test(Y, X, X.index, resample=FALSE, n.replicates=NULL)$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test(Y, X, X.index, resample=TRUE, n.replicates=n.resample)
              pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval) = "Asymptotic"
      score.tmp = pval[1,]
      
    }else{
      
      rownames(pval) = c("Asymptotic", "Resampling")
      score.tmp = pval[2,]
    }
    
    #print(pval)
    
    # identify significant lineages
    subtree.tmp = subtree
    index.na = which(is.na(score.tmp))
    if(length(index.na)>0){
      score.tmp = score.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }
    
    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.tmp)
    
    # Benjamini-Hochberg FDR control
    index.p = order(score.tmp)
    p.sort = sort(score.tmp)
    #fdr.alpha = 0.05
    
    # change 04/17/2016
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1      
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    
    # perform global test
    global.score.fisher = .F.test(score.tmp)
    global.score.min =  .simes.test(score.tmp)
    global.pval = c(global.score.fisher, global.score.min)
    names(global.pval) = c("Fisher", "Simes")
    
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
        
    
  }
  
}

#results.two = QCAT_GEE(count.rff, X, 1, X, 1, tax, n.resample=1000, fdr.alpha=0.05)
QCAT_GEE <- function(OTU, X, X.index, Z, Z.index, Tax=NULL, min.depth=0, n.resample=NULL, fdr.alpha=0.05){
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  Z = cbind(1, Z) # add the intercept term
  Z.index = Z.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      tmp = .Score.test2(count, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL)
      pval.comb = as.matrix( tmp$comb.results$score.pvalue )
      pval.zero = as.matrix( tmp$zero.results$score.pvalue )
      pval.pos = as.matrix( tmp$pos.results$score.pvalue )
      colnames(pval.comb) = "Asymptotic"
      colnames(pval.zero) = "Asymptotic"
      colnames(pval.pos) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test2(count, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample)
      pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
      pval.zero = c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue)
      pval.pos = c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue)
      names(pval.comb) = c("Asymptotic", "Resampling")
      names(pval.zero) = c("Asymptotic", "Resampling")
      names(pval.pos) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval.comb, pval.zero=pval.zero, pval.pos=pval.pos) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    tax = Tax[keep,,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    
    n.level = n.rank-1
    
    subtree = NULL
    pval.comb = NULL
    pval.zero = NULL
    pval.pos = NULL
    
    for(k in 1:n.level){
      
      #print(k)
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              tmp = .Score.test2(Y, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL)   
              
              pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
              pval.zero = cbind(pval.zero, tmp$zero.results$score.pvalue)
              pval.pos = cbind(pval.pos, tmp$pos.results$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test2(Y, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample)   
              
              
              pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
              pval.zero = cbind(pval.zero, c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue) )
              pval.pos = cbind(pval.pos, c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue) )
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval.comb) = subtree
    colnames(pval.zero) = subtree
    colnames(pval.pos) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval.comb) = "Asymptotic"
      score.comb.tmp = pval.comb[1,]
      
      rownames(pval.zero) = "Asymptotic"
      score.zero.tmp = pval.zero[1,]
      
      rownames(pval.pos) = "Asymptotic"
      score.pos.tmp = pval.pos[1,]
      
    }else{
      
      rownames(pval.comb) = c("Asymptotic", "Resampling")
      score.comb.tmp = pval.comb[2,]
      
      rownames(pval.zero) = c("Asymptotic", "Resampling")
      score.zero.tmp = pval.zero[2,]
      
      rownames(pval.pos) = c("Asymptotic", "Resampling")
      score.pos.tmp = pval.pos[2,]
    }
    
    #print(pval.comb)
    
    # identify significant lineages
    
    global.pval = NULL
    
    sig.lineage <- vector("list",3)
    names(sig.lineage) <- c("Two-Part", "Zero-Part", "Positive-Part")
    
    for(i in 1:3){
      
      if(i==1){score.tmp = score.comb.tmp}
      if(i==2){score.tmp = score.zero.tmp}
      if(i==3){score.tmp = score.pos.tmp}
      
      subtree.tmp = subtree
      index.na = which(is.na(score.tmp))
      if(length(index.na)>0){
        score.tmp = score.tmp[-index.na]
        subtree.tmp = subtree.tmp[-index.na]
      }
      
      #score.tmp[score.tmp==0] = 1e-4
      m.test = length(score.tmp)
      
      # Benjamini-Hochberg FDR control
      index.p = order(score.tmp)
      p.sort = sort(score.tmp)
      #fdr.alpha = 0.05
      
      # change 04/17/2016
      reject = rep(0, m.test)
      tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
      if(length(tmp)>0){
        index.reject = index.p[1:max(tmp)]
        reject[index.reject] = 1      
      }
      
      sig.lineage[[i]] = subtree.tmp[reject==1]
      
      # perform global test
      global.score.min =  .simes.test(score.tmp)
      global.pval = c(global.pval, global.score.min)
      
      
    }
    
    
    
    names(global.pval) = c("Simes_Two-Part", "Simes_Zero-Part", "Simes_Positive-Part")
    
    pval = list(pval.comb, pval.zero, pval.pos)
    names(pval) = c("Two-Part", "Zero-Part", "Positive-Part")
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}




########################################
#                                      #
#             ZIGDM test               #
#             add in v2                #
#             10/17/2017               #
#                                      #
########################################

ZIGDM <- function(OTU, X4freq, X4mean, X4disp, test.type="Mean", X.index, ZI.LB=10, Tax=NULL, min.depth=0, n.resample=NULL, fdr.alpha=0.05){
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.null(X4freq) & !is.matrix(X4freq)){
    warning("Covariate table X4freq is not a matrix")
    X4freq = as.matrix(X4freq)  
  }  

  if(!is.null(X4mean) & !is.matrix(X4mean)){
    warning("Covariate table X4mean is not a matrix")
    X4mean = as.matrix(X4mean)  
  }  
  
  if(!is.null(X4disp) & !is.matrix(X4disp)){
    warning("Covariate table X4disp is not a matrix")
    X4disp = as.matrix(X4disp)  
  }  
  
  flag = 0
  if( !is.null(X4freq) ){
    if(nrow(OTU)!=nrow(X4freq)){
      flag = 1
    }
  }
  if( !is.null(X4mean) ){
    if(nrow(OTU)!=nrow(X4mean)){
      flag = 1
    }
  }
  if( !is.null(X4disp) ){
    if(nrow(OTU)!=nrow(X4disp)){
      flag = 1
    }
  }
  
  if(flag){
    stop("Samples in the OTU table and the covariate tables should be the same.")  
  }
  
  if( test.type == "Freq" & is.null(ZI.LB)){
    stop("Test for differential presence-absence frequency cannot be performed under GDM model (ZI.LB is NULL).")  
  }
  
  if( (test.type == "Mean" & is.null(X4mean)) | (test.type == "Disp" & is.null(X4disp)) | (test.type == "Freq" & is.null(X4freq)) | (test.type == "Omni" & (is.null(X4mean) | is.null(X4disp) | is.null(X4freq)) ) ){
    stop("The specified test cannot be performed because the required covariate matrix is not provided.")
  }
  
  if( test.type == "Omni" & (!identical(X4mean, X4disp) | !identical(X4mean, X4freq) ) ){
    stop("X4mean, X4disp, and X4freq should be the same in order to perform Omnibus test.")
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X4freq = X4freq[-remove.subject,  ,drop=FALSE]
    X4mean = X4mean[-remove.subject,  ,drop=FALSE]
    X4disp = X4disp[-remove.subject,  ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  
  n = nrow(count)
  
  if( is.null(X4freq) ){
    X4freq = matrix(1, nrow=n, 1)
  }else{
    X4freq = cbind(1, X4freq) # add the intercept term
  }
  
  if( is.null(X4mean) ){
    X4mean = matrix(1, nrow=n, 1)
  }else{
    X4mean = cbind(1, X4mean) # add the intercept term
  }
  
  if( is.null(X4disp) ){
    X4disp = matrix(1, nrow=n, 1)
  }else{
    X4disp = cbind(1, X4disp) # add the intercept term
  }
  
  
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    

    aa = order( colSums(count), decreasing = TRUE )
    Y = count[,c(aa[-1],aa[1])]
    
    if(is.null(ZI.LB)){
      zindex = NULL
    }else{
      zindex = which(colSums(count[,1:(ncol(count)-1),drop=FALSE]==0)>ZI.LB)
      if(length(zindex)==0){
        zindex = NULL
      } 
    }
    
    if(is.null(n.resample)){ # asymptotic test only
      
      if(test.type == "Freq"){
        tmp = .ZIGDM.freq.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.freq )
      }
      if(test.type == "Mean"){
        tmp = .ZIGDM.abun.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.abun )
      }      
      if(test.type == "Disp"){
        tmp = .ZIGDM.disp.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.disp )
      } 
      if(test.type == "Omni"){
        tmp = .ZIGDM.omni.test.PAR2(Y, X4mean, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.omni )
      }
      
      colnames(pval) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      if(test.type == "Freq"){
        tmp = .ZIGDM.freq.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.freq, tmp$score.pval.perm.freq )
      }
      if(test.type == "Mean"){
        tmp = .ZIGDM.abun.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.abun, tmp$score.pval.perm.abun )
      }      
      if(test.type == "Disp"){
        tmp = .ZIGDM.disp.test.PAR2(Y, X4freq, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.disp, tmp$score.pval.perm.disp )
      } 
      if(test.type == "Omni"){
        tmp = .ZIGDM.omni.test.PAR2(Y, X4mean, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.omni, tmp$score.pval.perm.omni )
      }
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    
    tax = Tax[keep, ,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            aa = order( colSums(Y), decreasing = TRUE )
            Y2 = Y[,c(aa[-1],aa[1])]
            
            if(is.null(ZI.LB)){
              zindex = NULL
            }else{
              zindex = which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>ZI.LB)
              if(length(zindex)==0){
                zindex = NULL
              } 
            }
            
            # remove subject with depth ==0
            remove.sub.index = which(rowSums(Y2)==0)
            
            if(length(remove.sub.index)==0){
              X2.freq = X4freq
              X2.mean = X4mean
              X2.disp = X4disp
              
            }else{
              
              Y2 = Y2[-remove.sub.index, ,drop=FALSE]
              X2.freq = X4freq[-remove.sub.index, ,drop=FALSE]
              X2.mean = X4mean[-remove.sub.index, ,drop=FALSE]
              X2.disp = X4disp[-remove.sub.index, ,drop=FALSE]
              
            }
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              if(test.type == "Freq"){
                tmp = .ZIGDM.freq.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.freq )
              }
              if(test.type == "Mean"){
                tmp = .ZIGDM.abun.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.abun )
              }      
              if(test.type == "Disp"){
                tmp = .ZIGDM.disp.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.disp )
              } 
              if(test.type == "Omni"){
                tmp = .ZIGDM.omni.test.PAR2(Y2, X2.mean, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.omni )
              }
              
              
              
            }else{ # resampling test + asymptotic test 
              
              if(test.type == "Freq"){
                tmp = .ZIGDM.freq.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.freq, tmp$score.pval.perm.freq) )
              }
              if(test.type == "Mean"){
                tmp = .ZIGDM.abun.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.abun, tmp$score.pval.perm.abun) )
              }      
              if(test.type == "Disp"){
                tmp = .ZIGDM.disp.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.disp, tmp$score.pval.perm.disp) )
              } 
              if(test.type == "Omni"){
                tmp = .ZIGDM.omni.test.PAR2(Y2, X2.mean, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.omni, tmp$score.pval.perm.omni) )
              }
             
              
            }
            
          
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval) = "Asymptotic"
      score.tmp = pval[1,]
      
    }else{
      
      rownames(pval) = c("Asymptotic", "Resampling")
      score.tmp = pval[2,]
    }
    
    #print(pval)
    
    # identify significant lineages
    subtree.tmp = subtree
    index.na = which(is.na(score.tmp))
    if(length(index.na)>0){
      score.tmp = score.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }
    
    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.tmp)
    
    # Benjamini-Hochberg FDR control
    index.p = order(score.tmp)
    p.sort = sort(score.tmp)
    #fdr.alpha = 0.05
    
    # change 04/17/2016
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1      
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    
    # perform global test
    global.score.fisher = .F.test(score.tmp)
    global.score.min =  .simes.test(score.tmp)
    global.pval = c(global.score.fisher, global.score.min)
    names(global.pval) = c("Fisher", "Simes")
    
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}



.lbeta2 <- function(a,b){
  a[a<0] = 0
  b[b<0] = 0
  return(lbeta(a,b))
}

# parameters of posterior probability P_1...P_K | Y_1...Y_K Y_{K+1}
.rP.Y.par <- function(pv, av, bv, Y){
  
  K = length(Y)-1
  N = sum(Y)
  
  pv.post = rep(0, K)
  
  av.prim = av + Y[1:K]  
  bv.prim = bv
  for(j in 1:K){
    bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))
    
    if(Y[j]==0){
      
      if(beta(av[j], bv[j])==0){
        
        pv.post[j] = pv[j]/(pv[j] + (1-pv[j]))
        
      }else{
        tmp = pv[j] + (1-pv[j])*beta(av.prim[j], bv.prim[j])/beta(av[j], bv[j])
        if(tmp==0){
          pv.post[j] = 1
        }else{
          pv.post[j] = pv[j]/tmp
        }
        
      }
      
      
    }
    
  }
  
  return(list(pv.post = pv.post, av.post = av.prim, bv.post = bv.prim))
}

# maximize the objective function for logistic regression
# par = Logit.par.ini
# data = Logit.data
.Logit.neg.loglik <- function(par, data){
  
  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  tmp = data$Del * log(p) + (1-data$Del) * log(1-p)
  index = which(p==0 | p==1)
  if(length(index)>0){
    tmp[index] = 0
  }
  #tmp[is.na(tmp)] = 0   ## remove in 12/28/2016 will make optim function output a large estimate  
  return( -sum( tmp) )
  
}

.Logit.neg.score <- function(par, data){
  
  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  return( -colSums( (data$Del - p) * data$X ) )
  
}


#.Logit.optim(Del.R[,j], X.null, c.last[j, ])
#X=X.null
#Del = Del.R[,j]
#c.ini = c.last[j,]
.Logit.optim <- function(Del, X, c.ini){
  
  Logit.par.ini = c.ini
  Logit.data = list(Del=Del, X=X)
  
  #.Logit.neg.loglik(Logit.par.ini, Logit.data)
  #.Logit.neg.score(Logit.par.ini, Logit.data)
  return( optim(par=Logit.par.ini, fn=.Logit.neg.loglik, gr=.Logit.neg.score, data = Logit.data, method="BFGS")$par)
  #return( optim(par=Logit.par.ini, fn=.Logit.neg.loglik, gr=.Logit.neg.score, data = Logit.data)$par)
  
}


# maximize the objective function for Beta regression
# par = Beta.par.ini
# data = Beta.data
.Beta.neg.loglik.PAR2 <- function(par, data){
  
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  
  tmp = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp/(1+tmp) )
  
  tmp = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp/(1+tmp) )
  
  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp) 
  a[a<0] = 0
  b[b<0] = 0
  
  return( -sum( (1-data$Del) * ( -.lbeta2(a,b) + data$A * (a-1) + data$B *(b-1) )  ) )
  
}


.Beta.neg.score.PAR2 <- function(par, data){
  
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  tmp.a = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
  
  tmp.b = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp.b/(1+tmp.b) )
  
  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp) 
  a[a<0] = 0
  b[b<0] = 0
  
  a.a = (1/tmp.b) * mu.tmp * (1/(1+tmp.a))
  a.b = -(1/tmp.b) * mu.tmp
  b.a = - a.a
  b.b = -(1/tmp.b) *  (1/(1+tmp.a))
  
  one = (1-data$Del)*(digamma(a+b)-digamma(a) + data$A)
  two = (1-data$Del)*(digamma(a+b)-digamma(b) + data$B)
  
  return( -c( colSums( (one*a.a + two*b.a) * data$Xa ), colSums( (one*a.b + two*b.b) * data$Xb ) ) )
}


.Beta.optim.PAR2 <- function(Del, A, B, Xa, Xb, alpha.ini, beta.ini){
  
  Beta.par.ini = c(alpha.ini, beta.ini) 
  Beta.data = list(Del=Del, A=A, B=B, Xa=Xa, Xb=Xb)
  
  #Beta.neg.loglik(Beta.par.ini, Beta.data)
  #Beta.neg.score(Beta.par.ini, Beta.data)
  return( optim(par=Beta.par.ini, fn=.Beta.neg.loglik.PAR2, gr=.Beta.neg.score.PAR2, data = Beta.data, method="BFGS")$par)
  
}


################################# Reparameterize: Estimation and Testing #################################

#make a simulated data
# EM algorithm
# Y: taxa count matrix (length: K+1 )
# X: design matrix with the first component as intercept (dim: n x d)
# c0: initial value for c (dim: K x d)
# alpha0: initial value for alpha (dim: K x d)
# beta0: initial value for beta (dim: K x d)
# tol: convergence criterion (scalar 0-1)
# max.iter: maximum number of iteration (scalar integer)
#ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)

.ZIGDM.EM.PAR2 <- function(Y, Xc, Xa, Xb, c0, alpha0, beta0, tol, max.iter){
  
  CONV = 0
  CONV.iter = max.iter
  n = nrow(Y)
  K = ncol(Y)-1
  dc = ncol(Xc)
  da = ncol(Xa)
  db = ncol(Xb)
  
  #c.all = list(c0); alpha.all = list(alpha0); beta.all = list(beta0); # save estimate from all the iterations (can be comment out in the final code)
  c.last = c0; alpha.last = alpha0; beta.last = beta0;
  c.now = c0; alpha.now = alpha0; beta.now = beta0;
  
  Del.R = matrix(NA, n, K)
  A.R = matrix(NA, n, K)
  B.R = matrix(NA, n, K)
  
  for(l in 1:max.iter){
    
    # E-step
    
    # print(paste("====== ", l, "th ======", sep=""))
    for(i in 1:n){
      
      tmp = exp(c.last %*% Xc[i,])
      tmp[is.na(tmp)] = 0    
      pv = as.numeric( tmp/(1+tmp) )
      pv[is.infinite(tmp) & tmp>0] =1  # positive inf
      
      tmp = exp(alpha.last %*% Xa[i,])
      mv = as.numeric( tmp/(1+tmp) )
      tmp = exp(beta.last %*% Xb[i,])
      sv = as.numeric( tmp/(1+tmp) )
      
      av = (1/sv - 1) * mv 
      bv = (1/sv - 1) * (1-mv)  
      
      par.post = .rP.Y.par(pv, av, bv, Y[i,])
      
      Del.R[i, ] = par.post$pv.post
      tmp = par.post$av.post + par.post$bv.post
      A.R[i, ] = digamma(par.post$av.post) - digamma(tmp)
      B.R[i, ] = digamma(par.post$bv.post) - digamma(tmp)
      
    }
    
    
    
    # M-step
    for(j in 1:K){
      
      if( !is.infinite(c.last[j,1]) ){
        c.now[j, ] = .Logit.optim(Del.R[,j], Xc, c.last[j, ])
      }
      
      tmp = .Beta.optim.PAR2(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[j, ], beta.last[j, ])
      alpha.now[j, ] = tmp[1:da]; beta.now[j, ] = tmp[-(1:da)]
      
    }
    
    
    diff = 0
    if( sum(!is.infinite(c.now))>0 ){
      
      diff = diff + sum(abs(c.now[!is.infinite(c.now)]-c.last[!is.infinite(c.now)])) 
      
    }
    
    diff = diff + sum(abs(alpha.now-alpha.last))
    diff = diff + sum(abs(beta.now-beta.last))
    if( diff < tol ){  
      
      CONV = 1; CONV.iter = l;
      break;
      
    }else{
      
      c.last = c.now; alpha.last = alpha.now; beta.last = beta.now;
    }
    
  }
  
  return(list(c.est = c.now, alpha.est = alpha.now, beta.est = beta.now, CONV=CONV, CONV.iter = CONV.iter))
  
}


.ZIGDM.stat.omni.perm <- function(X.perm, par.interest.index, S.list, I.list){
  
  n = nrow(X.perm)
  d = ncol(X.perm)
  K3 = length(S.list[[1]])
  Kd3 = K3 * d
  
  ES = rep(0,Kd3)
  I = matrix(0, Kd3, Kd3)
  
  for(i in 1:n){
    
    ES = ES + kronecker(S.list[[i]], X.perm[i,])
    X2 = X.perm[i,] %o% X.perm[i,]
    I = I + kronecker( I.list[[i]] , X2)     
    
  }
  
  S1 = ES[par.interest.index]
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  
  return(score.stat)
  
}

#Y = Y.sim
#X = X.sim
#X.index = 2
#nperm = 1000

# default: only doing asymptotic test
# default: zero.index is the taxa with at least 1 zero observation

# ZIGDM.omni.test(Y2, X2, 2, zero.index=which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>0), nperm=100) 
# Y = Y2
# X = X2
# X.index = 2
# zero.index=which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>0)
# zero.index = NULL
#ZIGDM.omni.test(Y2, X2, 2, zero.index=NULL, nperm=100)
.ZIGDM.omni.test.PAR2 <- function(Y, X, X.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  X <- as.matrix(X[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  d = ncol(X)
  
  d.interest = length(X.index)
  X.null = X[,-X.index, drop=FALSE]
  d.null = ncol(X.null)
  
  # estimate parameters under the null model (omnibus test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, d.null) 
    
  }else{
    c0 = matrix(0, K, d.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, d.null) 
  beta0 = matrix(0, K, d.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, X.null, X.null, X.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  Kd3 = K * d * 3
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% X.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% X.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% X.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2 # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    X2 = X[i,] %o% X[i,]
    # I_1 = I_1 + kronecker(EI, X2)
    # I_2 = I_2 + kronecker(ESS, X2)
    # I_3 = I_3 + kronecker(ES.2, X2)
    I.list[[i]] = (-EI - ESS + ES.2)
    I = I + kronecker( I.list[[i]] , X2)  
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  tmp = c(X.index, X.index+d, X.index+2*d)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+3*d
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  S1 = ES[par.interest.index]
  #I = - I_1 - I_2 + I_3
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat.omni = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  score.pval.omni = 1-pchisq(score.stat.omni, length(par.interest.index))
  
  
  if(!is.null(nperm)){
    
    set.seed(11)
    count.omni = 0
    
    X.perm = X
    for(l in 1:nperm){
      if(is.null(strata)){
        X.perm[,X.index] = X[sample(1:n),X.index]
      }else{
        X.perm[,X.index] = X[.sampling.strata(1:n, strata),X.index]
        
      }
      tmp.omni = .ZIGDM.stat.omni.perm(X.perm, par.interest.index, S.list, I.list)
      
      if(tmp.omni >= score.stat.omni){
        count.omni = count.omni + 1
      }      
      
    }
    
    if(count.omni<=2){
      warnings("permutation p-value for omnibus test may not be accurate.")
    }
    
    score.pval.perm.omni = (count.omni+1)/(nperm+1)
    
  }else{
    
    score.pval.perm.omni = NA
    
  }
  
  
  return(list(score.stat.omni=score.stat.omni, score.pval.omni = score.pval.omni, score.pval.perm.omni = score.pval.perm.omni, df=df))
  
}

## FOR permutation test of each component
# if target == "freq", then only test for frequency
# if target == "abun", then only test for abundance
# if target == "disp", then only test for dispersion
# ZIGDM.stat.target.perm(Xc, Xa, Xb.perm, par.interest.index, S.list, I.list)
.ZIGDM.stat.target.perm <- function(Xc, Xa, Xb, par.interest.index, S.list, I.list){
  
  n = nrow(Xc)
  dc = ncol(Xc)
  da = ncol(Xa)
  db = ncol(Xb)
  dd = da + db + dc
  K3 = length(S.list[[1]])
  K = K3/3
  Kd3 = K * dd
  
  ES = rep(0,Kd3)
  I = matrix(0, Kd3, Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  
  for(i in 1:n){
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])   
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,])   
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,])
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,])
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])
      
      
    }
    
    
  } 
  
  
  
  S1 = ES[par.interest.index]
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  
  return(score.stat)
  
}


.ZIGDM.stat.freq.perm <- function(Y, Xc, Xa, Xb, Xc.index, zero.index=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  dc.interest = length(Xc.index)
  Xc.null = Xc[,-Xc.index, drop=FALSE]
  dc.null = ncol(Xc.null)
  
  # estimate parameters under the null model (absent frequency test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc.null) 
    
  }else{
    c0 = matrix(0, K, dc.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc.null, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xc.index)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.freq.test.PAR2 <- function(Y, Xc, Xa, Xb, Xc.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xa)>1 | ncol(Xb)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  dc.interest = length(Xc.index)
  Xc.null = Xc[,-Xc.index, drop=FALSE]
  dc.null = ncol(Xc.null)
  
  # estimate parameters under the null model (absent frequency test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc.null) 
    
  }else{
    c0 = matrix(0, K, dc.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc.null, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xc.index)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    
    if(!is.null(nperm) & nuis == 0){
      
      set.seed(11)
      count = 0
      
      Xc.perm = Xc
      for(l in 1:nperm){
        
        
        if(is.null(strata)){
          Xc.perm[,Xc.index] = Xc[sample(1:n),Xc.index]
        }else{
          Xc.perm[,Xc.index] = Xc[.sampling.strata(1:n, strata),Xc.index]
          
        }
        tmp = .ZIGDM.stat.target.perm(Xc.perm, Xa, Xb, par.interest.index, S.list, I.list)
        
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
      
    }else if(!is.null(nperm) & nuis == 1){
      
      set.seed(11)
      count = 0
      
      Xc.perm = Xc
      for(l in 1:nperm){
        
        
        if(is.null(strata)){
          perm.idx = sample(1:n)
        }else{
          perm.idx = .sampling.strata(1:n, strata)
        }
        Xc.perm = Xc[perm.idx,]
        Xa.perm = Xa[perm.idx,]
        Xb.perm = Xb[perm.idx,]
        
        tmp = .ZIGDM.stat.freq.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xc.index, zero.index=zero.index)
        
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
      
    }else{
      
      score.pval.perm = NA
      
    }
    
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  
  return(list(score.stat.freq=score.stat, score.pval.freq = score.pval, score.pval.perm.freq = score.pval.perm, df=df))
  
}


.ZIGDM.stat.abun.perm <- function(Y, Xc, Xa, Xb, Xa.index, zero.index=NULL){
  
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  da.interest = length(Xa.index)
  Xa.null = Xa[,-Xa.index, drop=FALSE]
  da.null = ncol(Xa.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da.null) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa.null, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xa.index+dc)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.abun.test.PAR2 <- function(Y, Xc, Xa, Xb, Xa.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xc)>1 | ncol(Xb)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  da.interest = length(Xa.index)
  Xa.null = Xa[,-Xa.index, drop=FALSE]
  da.null = ncol(Xa.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da.null) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa.null, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xa.index+dc)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    if(!is.null(nperm) & nuis == 0){
      
      set.seed(11)
      count = 0
      
      Xa.perm = Xa
      for(l in 1:nperm){
        
        if(is.null(strata)){
          Xa.perm[,Xa.index] = Xa[sample(1:n),Xa.index]
        }else{
          Xa.perm[,Xa.index] = Xa[.sampling.strata(1:n, strata),Xa.index]
          
        }
        
        tmp = .ZIGDM.stat.target.perm(Xc, Xa.perm, Xb, par.interest.index, S.list, I.list)
        
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
    }else if(!is.null(nperm) & nuis == 1){
      
      set.seed(11)
      count = 0
      
      Xa.perm = Xa
      for(l in 1:nperm){
        
        if(is.null(strata)){
          perm.idx = sample(1:n)
        }else{
          perm.idx = .sampling.strata(1:n, strata)
        }
        
        Xa.perm = Xa[perm.idx,]
        Xb.perm = Xb[perm.idx,]
        Xc.perm = Xc[perm.idx,]
        
        tmp = .ZIGDM.stat.abun.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xa.index, zero.index=zero.index)
        
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
    }else{
      
      score.pval.perm = NA
      
    }
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  
  return(list(score.stat.abun=score.stat, score.pval.abun = score.pval, score.pval.perm.abun = score.pval.perm, df=df))
  
}


.ZIGDM.stat.disp.perm <- function(Y, Xc, Xa, Xb, Xb.index, zero.index=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  db.interest = length(Xb.index)
  Xb.null = Xb[,-Xb.index, drop=FALSE]
  db.null = ncol(Xb.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa, Xb.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2) 
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xb.index+dc+da)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.disp.test.PAR2 <- function(Y, Xc, Xa, Xb, Xb.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xc)>1 | ncol(Xa)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  db.interest = length(Xb.index)
  Xb.null = Xb[,-Xb.index, drop=FALSE]
  db.null = ncol(Xb.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa, Xb.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2) 
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xb.index+dc+da)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    
    if(!is.null(nperm) & nuis==0){
      
      set.seed(11)
      count = 0
      
      Xb.perm = Xb
      #score.stat.perm = rep(0, nperm)
      for(l in 1:nperm){
        
        if(is.null(strata)){
          Xb.perm[,Xb.index] = Xb[sample(1:n),Xb.index]
        }else{
          Xb.perm[,Xb.index] = Xb[.sampling.strata(1:n, strata),Xb.index]
          
        }
        tmp = .ZIGDM.stat.target.perm(Xc, Xa, Xb.perm, par.interest.index, S.list, I.list)
        #score.stat.perm[l] = tmp
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
    }else if(!is.null(nperm) & nuis==1){
      
      set.seed(11)
      count = 0
      
      Xb.perm = Xb
      #score.stat.perm = rep(0, nperm)
      for(l in 1:nperm){
        
        if(is.null(strata)){
          perm.idx = sample(1:n)
          
        }else{
          perm.idx = .sampling.strata(1:n, strata)
          
          
        }
        Xc.perm = Xc[perm.idx,]
        Xa.perm = Xa[perm.idx,]
        Xb.perm = Xb[perm.idx,]
        
        tmp = .ZIGDM.stat.disp.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xb.index, zero.index = zero.index)
        #score.stat.perm[l] = tmp
        if(tmp >= score.stat){
          count = count + 1
        }      
        
      }
      
      if(count<=2){
        warnings("permutation p-value for frequency test may not be accurate.")
      }
      
      score.pval.perm = (count+1)/(nperm+1)
      
    }else{
      
      score.pval.perm = NA
      
    }
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  return(list(score.stat.disp=score.stat, score.pval.disp = score.pval, score.pval.perm.disp = score.pval.perm, df=df))
  
}

