library(MASS)
library(TMB)
library(data.table)
dyn.load(dynlib(paste0(system.file("./rates/src/TMB",package="rates"),"/tau2f")))
#' Fit the RATES model.
#' @param betajk MxK matrix of effect size estimates, where row j corresponds to j-th SNP, and column k corresponds to k-th study .  Missing values can be represented with NA.
#' @param sjk2 MxK matrix of effect size estimate variances, where row j corresponds to SNP j, and column k corresponds to study k.  Missing values can be represented with NA.
#' @param PC Kx4 matrix of PCs (PC0,PC1,PC2,PC3) of genome-wide AF for K studies.
#' @param p initial value for EM algorithm, for the proportion of non-zero SNPs
#' @param lambda initial value for EM algorithm, for the proportion (P(O_jk=1|Rj=0)) of non-replicable SNPs which are outliers with inflated variance.
#' @param tau2 initial value for EM algorithm, a (L+1)-vector for the variance of coefficients of PCs (gamma_jl,l= 0,1,2,3,if L=3) for replicable non-zero effect SNPs.
#' @param alpha initial value for EM algorithm, for the variance inflation of outlier summary statistics at non-replicable SNPs.
#' @param rel.eps threshold for when to end EM algorithm. rel.eps = (ll[i] - ll[i-1])/ll[i], where ll[i] is the log-likelihood at iteration i.
#' @param maxIter maximum # of EM iterations.
#' @param minIter minimum # of EM iterations.
#' @return
#' A list.
#' @export
rates.fit<-function(
  betajk,sjk2,PC,
  p=0.003,
  lambda=0.04,
  tau2=NULL,
  alpha=3,
  SNP=NULL,
  rel.eps=1e-8,
  verbose=1,
  maxIter=10^4L,
  minIter=50){
  betajk<-as.matrix(betajk)
  sjk2<-as.matrix(sjk2)
  PC<-as.matrix(PC)
  Lp1<-ncol(PC)
  if (missing(tau2)) {
    tau2<-rep(0.0002,Lp1)
  }
  if (Lp1 != length(tau2)) {stop("unmatched dimension for PC matrix and tau2")}
  ###Handle missing value###
  ##replace sd2=0 with NA
  zeroL<-lapply(1:nrow(sjk2), function(j){
    which(sjk2[j,]==0)
  })
  chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2[j,zeroL[[j]]]<-NA
      betajk[j,zeroL[[j]]]<-NA
    }
  }
  ## Replace any summary stats with standard error = inf or beta=inf with missing
  infL<-lapply(1:nrow(sjk2), function(j){
    which(is.infinite(sjk2[j,]))
  })
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2[j,infL[[j]]]<-NA
      betajk[j,infL[[j]]]<-NA
    }
  }
  infL<-lapply(1:nrow(betajk), function(j){
    which(is.infinite(betajk[j,]))
  })
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2[j,infL[[j]]]<-NA
      betajk[j,infL[[j]]]<-NA
    }
  }
  ### Identify indices of studies with non-missing summary stats at each SNP
  mis.inds<-lapply(1:nrow(sjk2), function(j){
    sdjk2.non.miss<-which(!is.na(sjk2[j,]))
    beta.non.miss<-which(!is.na(betajk[j,]))
    return(intersect(sdjk2.non.miss,beta.non.miss))
  })
  ######it's possible that sd2 (can be approximated) is not missing but beta is missing##
  ll<-NA
  # # of SNPs 
  MM<-nrow(betajk)
  ##prepare input for tau2f object function
  z2invsj <-lapply(1:MM,function(j){
    t(PC[mis.inds[[j]],])%*%diag(1/sjk2[j,mis.inds[[j]]])%*%PC[mis.inds[[j]],]
  })
  dj <-lapply(1:MM,function(j){
    t(PC[mis.inds[[j]],])%*%diag(1/sjk2[j,mis.inds[[j]]]) %*% betajk[j,mis.inds[[j]]]
  })
  
  # # of studies with non-missing summary stat for each variant
  k_j<-sapply(mis.inds, length)
  strt<-Sys.time()
  ## Start of EM algorithm
  for(iter in 1:(maxIter)){
    #print(dim(betajk))
    ## E-step start
    # calculate Ojk.p matrix
    Ojk.p<-lapply(1:MM,function(j){
      rates::deltis(betajk_j = betajk[j,],sdjk2_j = sjk2[j,],lambda = lambda,alpha = alpha)
    })
    Ojk.p<-matrix(unlist(Ojk.p),nrow = MM,byrow = TRUE)
    missingdelta<-unlist(lapply(1:MM, function(j){ 
      mean(is.na(Ojk.p[j,mis.inds[[j]]]))
    }))
    #print(missingdelta)
    if(sum(missingdelta) > 0) {
      print("na's in deltajk")
      break
    }
    #Calculate Rj.p vector
    if(iter==1){
      llbR0<-unlist(lapply(1:nrow(betajk), function(j){
        rates::bjR0(betajk_j = betajk[j,mis.inds[[j]]],sdjk2_j = sjk2[j,mis.inds[[j]]],alpha = alpha,lambda = lambda)
      }))
      
      llbR1<-unlist(lapply(1:nrow(betajk),function(j){
        rates::bjR1(betajk_j = betajk[j,mis.inds[[j]]],sdjk2_j = sjk2[j,mis.inds[[j]]],zkl =as.matrix(PC[mis.inds[[j]],]),tau2 = tau2)
      }))
    }
    #print(nrow(betajk))
    Rj.p<- unlist(lapply(1:nrow(betajk),function(j){
      #print(p)
      Rj_j<-p*exp(llbR1[j])/exp(rates::logsumexp(c(log(p)+llbR1[j],log(1-p)+llbR0[j])))
      if(is.na(Rj_j)){
        m<-which.max(c(log(1-p)+llbR0[j],log(p)+llbR1[j]))
        Rj_j<-(m-1)
      }
      #print(paste0("Rj_j length:",length(Rj_j)))
      return(Rj_j)
    }))
    #print(length(Rj.p))
    ### E-step finished
    if (verbose>0) print("e-step finished")
    #Deal with posterior prob as NA
    if(sum(is.na(Rj.p))>0) break
    if(sum(1-Rj.p)==0) break
    #print(dim(Ojk.p))
    ### M-step start
    # Update p
    p<-sum(Rj.p)/MM
    # Update lambda
    fn<-sum(unlist(lapply(1:MM,function(j){
      (1-Rj.p[j])*sum(Ojk.p[j,mis.inds[[j]]])
    })))
    fd<-sum((1-Rj.p)*k_j)
    lambda<-fn/fd
    #update alpha
    alpha.fn<-sum(unlist(lapply(1:MM,function(j){
      (1-Rj.p[j])*sum(Ojk.p[j,mis.inds[[j]]]*(betajk[j,mis.inds[[j]]])^2/sjk2[j,mis.inds[[j]]])
    })))
    alpha.fd<-as.vector(crossprod(1-Rj.p,rowSums(Ojk.p,na.rm = TRUE)))
    alpha<-alpha.fn/alpha.fd
    #print(betajk)
    ###########
    #update tau2
    nllk<-MakeADFun(data = list(z2invsj=z2invsj,dj=dj,gammaj=Rj.p),
                    parameters = list(logTau2=log(tau2)),
                    DLL="tau2f",silent = TRUE)
    fit <- nlminb ( start = nllk $par , objective = nllk $fn , gradient = nllk $gr ,
                    lower = -40 ,upper = Inf,control = list(trace=0))
    
    #fit<-nlminb(start = log(tau2),objective = nTau2f,
    #            betajk=betajk,sjk2=sjk2,PC=PC,Rj.p=Rj.p)
    tau2<-exp(fit$par)
    #tau2<-unname(sapply(tau2,max,1e-17))
    if(sum(is.na(tau2))>0) break
    if (verbose>0) print("m-step finished")
    ######finished review above######
    get.ll<-function(betajk,sjk2,PC,p,lambda,alpha,tau2,mis.inds){
      llbR0<-unlist(lapply(1:nrow(betajk), function(j){
        rates::bjR0(betajk_j =  betajk[j,mis.inds[[j]]],sdjk2_j = sjk2[j,mis.inds[[j]]],alpha = alpha,lambda=lambda)
      }))
      
      llbR1<-unlist(lapply(1:nrow(betajk),function(j){
        rates::bjR1(betajk_j =  betajk[j,mis.inds[[j]]],sdjk2_j = sjk2[j,mis.inds[[j]]],zkl = as.matrix(PC[mis.inds[[j]],]),tau2 = tau2)
      }))
      ll<-sum(unlist(lapply(1:nrow(betajk),function(j){
        rates::logsumexp(c(log(p)+llbR1[j],log(1-p)+llbR0[j]))
      })))
      return(list(ll=ll,llbR0=llbR0,llbR1=llbR1))
    }
    ll.list<-get.ll(betajk = betajk,sjk2 = sjk2,PC = PC,p = p,lambda = lambda,alpha = alpha,tau2 = tau2,mis.inds = mis.inds)
    ll[iter]<-ll.list$ll
    llbR0<-ll.list$llbR0
    llbR1<-ll.list$llbR1
    if(iter > 1){
      llMat<-rbindlist(list(llMat,
                            data.table(p=p,lambda=lambda,tau2=tau2,alpha=alpha,ll=ll[iter])))
    } else{
      llMat<-data.table(p=p,lambda=lambda,tau2=paste(tau2,collapse = ":"),alpha=alpha,ll=ll[iter])
    } 
    
    #Convergence criteria
    if(iter >=2){
      if(!is.na(ll[iter] - ll[iter-1])){
        if(((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && ll[iter] - ll[iter-1] < 0.01) || 
           ((ll[iter]-ll[iter-1]) < 0.001 && iter > minIter) || 
           ((ll[iter] - ll[iter-1])/ll[iter] < rel.eps && iter > minIter)){break}
      }
    }
  }
  end<-Sys.time()
  
  #####Posterior mean and variance for gamma and wald test
  post_var_list<-lapply(1:MM,
                        function(j){
                          raw.V_gamma_j<-diag(1/diag(t(PC[mis.inds[[j]],])%*%diag(1/sjk2[j,mis.inds[[j]]])%*%PC[mis.inds[[j]],]
                          +diag(1/tau2)))}
                        )
  post_mean_list<-lapply(1:MM,function(j){
    post_var_list[[j]]%*%dj[[j]]
  })
  ##Wald statistics
  wald_stat<-unlist(mapply(function(E_gamma_j,V_gamma_j)
    {t(E_gamma_j)%*%ginv(V_gamma_j)%*%E_gamma_j},
        post_mean_list,post_var_list))

  P_WALD<-pchisq(wald_stat,df = Lp1,lower.tail = F)
  ###Use wald test in R
  
  
  
  
  
  return(list(ll=ll,
              ppr=Rj.p,
              WALD_STAT=wald_stat,
              P_WALD=P_WALD,
              p=p,
              alpha=alpha,
              lambda=lambda,
              tau2=unname(tau2),
              E_gamma_j=post_mean_list,
              V_gamma_j=post_var_list,
              PC = PC,
              time=difftime(end,strt,units = "mins"),
              SNP=SNP,
              mis.inds=mis.inds
  ))
}





#'RATES generative model
#' @param mod a fitted rates model
#' @param sjk2 a M x K matrix of effect size variance used for fitting the original rates model/for rates generative model 
#' @param Mnull Approximate number of null variants to generate
#' @return 
#' a list of generated data
#' @export


rates.generative<-function(mod,sjk2,seed=NULL,Mnull=NULL){
  if(missing(seed)){
    seed<-sample.int(10^6,1)
  }
  sjk2<-as.matrix(sjk2)
  p.hat<-mod[["p"]]
  lambda.hat<-mod[["lambda"]]
  alpha.hat<-mod[["alpha"]]
  tau2.hat<-mod[["tau2"]]
  PC<-mod[["PC"]]
  Lp1<-ncol(PC)
  ########
  if(missing(Mnull)){
    Mnull<-ceiling(nrow(sjk2)*(1-p.hat))
  }
  ####
  M<-ceiling(Mnull/(1-p.hat))
  Rj<-rbinom(M,1,p.hat)
  gamma_jl.sampled<-t(sapply(Rj,mvrnorm,n=1,mu=rep(0,Lp1),Sigma=diag(tau2.hat)))
  gamma_jl<-diag(Rj)%*%gamma_jl.sampled
  mu_jk<-gamma_jl%*%t(PC)
  ###sample sjk2
  sample_inds<-sample(nrow(sjk2),M,replace = TRUE)
  sjk2_sample<-sjk2[sample_inds,]
  ##identify missing studies from sampled s2 matrix
  zeroL<-lapply(1:nrow(sjk2_sample), function(j){
    which(sjk2_sample[j,]==0)
  })
  chk<- which(sapply(zeroL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2_sample[j,zeroL[[j]]]<-NA
    }
  }
  ## Replace any summary stats with standard error = inf or beta=inf with missing
  infL<-lapply(1:nrow(sjk2_sample), function(j){
    #print(is.infinite(sjk2_sample[j,]))
    which(is.infinite(sjk2_sample[j,]))
  })
  chk<- which(sapply(infL, length) > 0)
  if(length(chk) > 0){
    for(j in chk){
      sjk2_sample[j,infL[[j]]]<-NA
    }
  }
  ### Identify indices of studies with non-missing summary stats at each SNP
  mis.inds_sample<-lapply(1:nrow(sjk2_sample), function(j){
    sdjk2.non.miss<-which(!is.na(sjk2_sample[j,]))
  })
  
  ###simulate Ojk
  #k_j<-sapply(1:nrow(s2_sample), function(j){sum(!is.na(s2_sample[j,]))})
  ##make the sampled s2 matrix a list by row
  sj2<-lapply(1:M, function(j){sjk2_sample[j,]})
  Ojk<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sjk2_sample)[2])
    vec[mis.inds_sample[[j]]]<-rbinom(length(mis.inds_sample[[j]]), size=1, prob=lambda.hat)
    vec
  })
  ###
  betaj<-lapply(1:M, function(j){
    vec<-rep(NA, dim(sjk2_sample)[2])
    vec[mis.inds_sample[[j]]]<-
      Rj[j] * rnorm(length(mis.inds_sample[[j]]), mu_jk[j,mis.inds_sample[[j]]], sd=sqrt(sj2[[j]][mis.inds_sample[[j]]])) + 
      (1-Rj[j]) * ((Ojk[[j]][mis.inds_sample[[j]]])*rnorm(length(mis.inds_sample[[j]]), 0, sqrt(alpha.hat)*sqrt(sj2[[j]][mis.inds_sample[[j]]]))+ 
                     (1-Ojk[[j]][mis.inds_sample[[j]]])*rnorm(length(mis.inds_sample[[j]]), 0, sqrt(sj2[[j]][mis.inds_sample[[j]]])))
    vec
  })
  betajk<-matrix(unlist(betaj), ncol=dim(sjk2_sample)[2], byrow=TRUE)
  return(list(betajk=betajk,sjk2_sample=sjk2_sample,
              sample_inds=sample_inds,Rj=Rj,Ojk=Ojk,label.seed=seed))
}

#test<-rates.generative(res,sjk2 = sd2.selected)
#test$sjk2_sample



#' parametric boostrap using RATES generative model
#' @param mod a fitted rates model
#' @param sjk2 a M x K matrix of effect size variance used for fitting the original rates model 
#' @param total_null_score total number of ppr of simulated null SNPs
#' @param seed seed used in bootstrap sampling
#' @param nullSNPsPerModel apprixamet number of null variants generated in each fold of parametric bootstrap
#' @return 
#' a list of null scores and used seed.
#' @export

rates.getnullscore<-function(mod,sjk2,total_null_scores,out.prefix,
                             seed=NULL,nullSNPsPerModel=NULL,clean.existing.score=F){
  if (missing(seed)) {
    seed<-sample.int(10^6,size = 1)
  }
  if(missing(nullSNPsPerModel)){
    nullSNPsPerModel<-ceiling(nrow(sjk2)*(1-mod$p))
  }
  if (clean.existing.score) {
    existing.file<-paste0(out.prefix,".nullscores.txt")
    log.file<-paste0(out.prefix,".nullscores.log")
    system(paste0("rm ", existing.file))
    system(paste0("rm ", log.file))
  }
  set.seed(seed)
  system(paste0("touch ", out.prefix, ".nullscores.txt"))
  total_score_count<-as.numeric(system(paste0("wc -l ", out.prefix, ".nullscores.txt | awk '{print $1}'"),
                                       intern=T))
  print(paste0(total_score_count, " scores already generated in ", out.prefix, ".nullscores.txt."))
  stopifnot(total_score_count < total_null_scores)
  nModels<- ceiling((total_null_scores-total_score_count) / nullSNPsPerModel)
  message(paste0("Fitting ", nModels, 
               " models to generate approximately ", total_null_scores ," non-replicable SNPs using seed:\n",seed))
  stop<-FALSE
  nullMod<-list()
  for (i in 1:nModels){
    gData<-rates.generative(mod =mod,sjk2 = sjk2,Mnull = nullSNPsPerModel)
    nullMod[[i]]<-rates::rates.fit(betajk = gData$betajk,sjk2 = gData$sjk2_sample,PC = mod$PC,
                                   p =mod$p,lambda =mod$lambda ,alpha =mod$alpha ,tau2 = mod$tau2,verbose=0)
    ppr_i<-nullMod[[i]]$ppr
    nullscore_i<-ppr_i[gData$Rj==0]
    fwrite(data.table(nullscore_i),
           file=paste0(out.prefix, ".nullscores.",gData$label.seed,".txt"),
           col.names = FALSE)
    system(paste0("cat ",
                  paste0(out.prefix, ".nullscores.",gData$label.seed,".txt")," >> ",
                  paste0(out.prefix, ".nullscores.txt")))
    system(paste0(" echo $(wc -l ", 
                  paste0(out.prefix, ".nullscores.txt) label seed ",gData$label.seed, " >> ",
                         paste0(out.prefix, ".nullscores.log"))))
    system(paste0("rm ",paste0(out.prefix, ".nullscores.",gData$label.seed,".txt")))
    
    message(paste0("Model ",i, " of ", nModels, " complete.\n"))
    
    total_score_count<-as.numeric(system(paste0("wc -l ", out.prefix, ".nullscores.txt| awk '{print $1}'"),intern=T))
    if (total_score_count>total_null_scores) stop<-TRUE
    if(stop) break
  }
  nullscores<-fread(paste0(out.prefix, ".nullscores.txt"),header=F)[,"V1"]
  return(list(nullscores=nullscores,seed=seed))
}


