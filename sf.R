#Glen Pridham 2023, Rutengroup Dalhousie University Physics Department
#script for fitting the SF model

# to do:
  #add bootstrap
    #needs complex validatation
  #add preprocessing function to subtract baseline mean and divide by baseline sd
    #needs validation
  #add a function to convert all parameters into natural variable parameters
    #needs validation

library(survival) #needed for tmerge and s
library(MASS) #needed for Sim.SF

################################# high-level functions ######################################
################################# high-level functions ######################################
################################# high-level functions ######################################
FitSF = function(y,
                 x=NULL, #covariates e.g. sex and age
                 dt = matrix(1,nrow=nrow(y),ncol=dim(y)[3]-1),
                 iter=5, #number of times to iterate, converges quickly
                 FitFun=FitMoments.SF,
                 fitOptions = list(), #note: diagonal option set in preamble
                 ErrorFun=AymErrorMoments.SF,
                 errorFunOptions = list(),
                 estimateErrors=T,
                 ImputeFirstFun = Impute, #for initial imputation
                 imputeFirstFunOptions = list(),
                 imputeFirst=T, #impute intrawave-conditional mean as first step? i.e. before computing anything else                 
                 ImputeEMFun = ImputeMean.SF, #for expectation-maximization (iterative) imputation
                 imputeEMFunOptions = list(),
                 imputeMean=T, #impute conditional mean as we iterate through? #has this been validated?
                 imputeAge=(imputeMean | imputeFirst), #impute age columns of x using dt?
                 imputeAgeCols = c("age"),
                 imputeDead=T, #controls whether or not you impute dead people
                 imputeCensored=T, #controls whether or not you impute censored people
                 age = NULL, #needed for informative censoring
                 s = NULL, #survival data, include if you want to use informative censoring
                 where = is.na(y), #where to impute (if at all)
                 Npc = NULL, #number of PCs to use for fitting; NULL to avoid PCA entirely
                 imputeInPCSpace = FALSE, #do you want to perform imputation in latent space? May imputed too many values #may give warning messages if F because can't invert reduced rank (happens occasionally)
                 pcaIndex=1, #time index to use for PCA
                 center=!is.null(Npc), #pre-processing step, center by mean
                 scale=F, #pre-processing step, scale by sd
                 addIntercept=T, #add column of 1 to x, first index
                 fixedQ = diag(1,ncol(y),ncol(y)), #set to NA if you want to estimate
                 diagonalW = F, #force W to be diagonal? works reasonably well when paired with PCA
                 verbose=T #prints updates for user (just a few)
)
{
  #generalization of FitSF.moments
  #structure:
  #1) initial imputation (optional) using ImputeFirstFun
  #2) transform using PCA (optional)
  #3) estimate parameters using FitFun
  #4) impute expectation-maximization estimates using parameters from 3)
  #5) repeat 3-4 iter times
  #6) if PCA used, transform back into observed space 
  #7) compute basic fit metrics
  
  #notes:
  #x and y must be 3d arrays
  #imputing dead/censored individuals can reduce bias for regularly sampled data, this requires a survival object, s
  #it is generally inadvisable to impute in PC-space
  
  #to do:
  #update command line parameters
  #update return values (shuld include command line parameters)
  #check preamble to see if it needs editing
  
  
  
  R2 = function(obs,est,na.rm=F)
  {
    #est: estimated y
    #obs: observed y

    r2 =  1 - mean((obs-est)^2,na.rm=na.rm)/mean((obs-mean(obs,na.rm=na.rm))^2,na.rm=na.rm) 

    return(r2)
  }
  
  if(sum(is.na(dt))>0) stop("NAs not permitted in dt")
  
  #diagonal option
  fitOptions[["wOptions"]][["diagonal"]] = diagonalW
  errorFunOptions[["wOptions"]][["diagonal"]] = diagonalW
  
  
  doNotImpute = array(F,dim=dim(y))
  if(!imputeDead)
  {
    if(verbose) print("Not imputing dead...")
    if(is.null(s)) stop("s must be supplied if you don't want to impute the dead")
    for (k in 1:dim(doNotImpute)[3]) for (j in 1:ncol(doNotImpute)) doNotImpute[,j,k][age[,k] > s[,"stop"] + 1e-10 & s[,"status"] == 1] = T
    #drop NAs too
    for (j in 1:ncol(doNotImpute)) doNotImpute[,j,][is.na(age)] = T
  }
  if(!imputeCensored)
  {
    if(verbose) print("Not imputing censored...")
    if(is.null(s)) stop("s must be supplied if you don't want to impute the censored")
    for (k in 1:dim(doNotImpute)[3]) for (j in 1:ncol(doNotImpute)) doNotImpute[,j,k][age[,k] > s[,"stop"] + 1e-10 & s[,"status"] == 0] = T
    #drop NAs too
    for (j in 1:ncol(doNotImpute)) doNotImpute[,j,][is.na(age)] = T
  }
  where[doNotImpute] = F
  
  x00 = x
  
  if(ncol(dt) >= dim(x)[3])
  {
    warning("dt has too many columns, check to make sure you've entered it correctly")
    dt = dt[,1:(dim(x)[3]-1)]
  }
  
  wherex = NULL
  if(imputeAge & !is.null(x))
  {
    wherex = is.na(x)
    for (j in 1:length(imputeAgeCols))
    {
      if (imputeAgeCols[j] %in% colnames(x))
      {
        if(verbose) print(sprintf("Imputing age in column %s",imputeAgeCols[j]))
        
        for (k in 2:dim(x)[3])
        {
          logi = is.na(x[,imputeAgeCols[j],k])
          x[logi,imputeAgeCols[j],k] = x[logi,imputeAgeCols[j],k-1]+dt[logi,k-1]
        }
        
        if(verbose)
        {
          print("Age imputation correct?")
          print(isTRUE(all.equal(dt,x[,imputeAgeCols[j],-1]-x[,imputeAgeCols[j],-dim(x)[3]],check.attributes=F)))
          if(!isTRUE(all.equal(dt,x[,imputeAgeCols[j],-1]-x[,imputeAgeCols[j],-dim(x)[3]],check.attributes=F)))
          {
            print(str(dt))
            print(str(x[,imputeAgeCols[j],-1]-x[,imputeAgeCols[j],-dim(x)[3]]))
            print(mean(dt==(x[,imputeAgeCols[j],-1]-x[,imputeAgeCols[j],-dim(x)[3]])))
          } 
        }

      }
    }
    
  }
  
  
  if(nrow(dt) != nrow(y)) stop("dt and y not same length (number of individuals)!")  
  
  y0 = y
  x0 = x
  
  #add one to x for estimating baseline mean (intercept)
  if(addIntercept)
  {
    if(!is.null(x0))
    {
      x = array(1,dim=c(dim(x0)[1],dim(x0)[2]+1,dim(x0)[3]))
      x[,-1,] = x0
      if(is.null(colnames(x0))) x0cols = sprintf("x%d",1:ncol(x0))
      else x0cols = colnames(x0)
      colnames(x) = c("mu0",x0cols)
    }
    else 
    {
      x = array(1,dim=c(dim(y)[1],1,dim(y)[3]-1))
      colnames(x) = "mu0"
    }
  }
  else if (is.null(x0)) stop("null x + no intercept is not yet supported (to do)")
  
  #print("xpp:")
  #print(dim(x))
  
  if(verbose) print(sprintf("Iterations selected: %.1f",iter))
  
  calcNoise = T
  fixedQ0 = fixedQ
  if(!all(is.na(fixedQ))) #won't calculate noise, just uses whatever fixedQ is
  {
    calcNoise = F
    Q = fixedQ
    if(verbose) print("Noise precision matrix has been supplied...")
  }
  else 
  {
    Q = diag(1,ncol(y),ncol(y))
  }
  
  if(verbose) print(sprintf("NAs to replace: %d, total NAs: %d (diff: %d)",sum(where),sum(is.na(y)),sum(is.na(y))-sum(where)))
  if(imputeFirst) #impute unconditioned mean as first step
  {
    y = ImputeFirstFun(y=y,x=x,dt=dt,where=where,options=imputeFirstFunOptions)[["imp"]]
    if(verbose) 
    {
      print(sprintf("NAs to replace: %d, NAs replaced by initial imputation: %d (imputed: %.0f%%)",sum(where),sum(!is.na(y[where])),100*mean(!is.na(y[where]))))
      print("initial imputation options:")
      print(imputeFirstFunOptions)
    }
  }
  
  mupp = rep(0,ncol(y))
  sdpp = rep(1,ncol(y))
  if(center)
  {
    mupp = apply(y[,,1,drop=F],2,mean,na.rm=T)
    yc = y
    for (j in 1:ncol(yc)) yc[,j,] = y[,j,] - mupp[j]
    y = yc
    rm(yc)
  }
  if(scale)
  {
    stop("scale not fully impemented yet, need to figure out how to scale parameters")
    sdpp = apply(y[,,1,drop=F],2,sd,na.rm=T)
    ys = y
    for (j in 1:ncol(ys)) ys[,j,] = (y[,j,])/sdpp[j]
    y = ys
    rm(ys)
  }
  
  
  doPCA = F
  #trans takes you to latent space
  #multiply vector on left or matrix on right
  trans = diag(1,ncol(y))
  transinv = t(trans)
  lpc = NULL
  pc = NULL
  if (!is.null(Npc))
  {
    if(!imputeFirst)
    {
      warning("You should consider an initial imputation since you're doing PCA (you can also run a full model without PCA, then use the imputed values from it!)")
    }
    if(!center) warning("You should center if you're doing PCA")
    doPCA = T
    #update: # allow multiple pcaIndexes #Dec 15, 2022
    yc = y[,,pcaIndex]
    flat = list()
    for (i in 1:length(pcaIndex)) flat[[i]] = yc
    flat=do.call(rbind,flat)
    
    pc = prcomp2(flat,center=F,scale.=F)
    
    #trans takes you to latent space
    #multiply vector on left or matrix on right
    trans = t(pc$rotation)
    transinv = pc$rotation
    
    #confirm transinv is correct
    tt = transinv%*%trans
    test = diag(1,ncol=ncol(tt),nrow=nrow(tt)) - tt
    if(sum(abs(test)) > 1e-6) 
    {
      print("this inverse doesn't look right. is this the identity matrix?")
      print(round(tt,3))
      stop("bad inverse")
    }
    
    trans = trans[1:Npc,,drop=F]
    transinv = transinv[,1:Npc,drop=F]
    
    for (k in 1:dim(y)[3]) 
    {
      y[,1:Npc,k] = y[,,k]%*%t(trans)
    }
    y = y[,1:Npc,,drop=F]
    colnames(y) = rownames(trans)[1:Npc]
    
    Q = diag(1/apply(y,2,sd,na.rm=T)^2,Npc,Npc)
    
    if(imputeInPCSpace)
    {
      where = is.na(y0)
      for (k in 1:dim(y0)[3]) 
      {
        where[,1:Npc,k] = is.na(y0[,,k]%*%t(trans))
      }
      where[doNotImpute] = F #censored/dead excludes
      where = where[,1:Npc,,drop=F]
      
      if(verbose) 
      {
        print("NAs in PC space")
        print(sum(is.na(y)))
      }
    }
  }
  
  #initial fit
  fit = FitFun(y=y,x=x,dt=dt,W=NULL,lambda=NULL,Q=Q,calcNoise=calcNoise,options=fitOptions)

  W = fit[["W"]]
  lambda = fit[["lambda"]]
  Q = fit[["Q"]]
  sigma = solve(fit[["Q"]])
  

  for (i in 1:iter)
  {
    if(verbose) cat(".")
    if (i < 1) break #because R is a bastard
    
    #(a) impute
    if(imputeMean) 
    {
      if(doPCA & !imputeInPCSpace) #must be done separately in case we pick Npc < initial size
      {
        ysave = y
        #undo transformations
        #general caess
        Wp = transinv[,1:Npc,drop=F]%*%W%*%trans[1:Npc,,drop=F]
        
        yp = y0
        for (k in 1:dim(yp)[3]) yp[,,k] = y[,,k]%*%t(transinv)[1:Npc,,drop=F]
        
        lambdap = transinv[,1:Npc,drop=F]%*%lambda
        
        Qp = transinv[,1:Npc,drop=F]%*%Q%*%trans[1:Npc,,drop=F]
        sigmap = transinv[,1:Npc,drop=F]%*%sigma%*%trans[1:Npc,,drop=F] #added June 2024
        #problem: this is probably not going to be invertible but ImputeMean.SF needs to invert a subset of sigmap...

        yp = ImputeEMFun(y=yp,x=x,dt=dt,where=where,W=Wp,lambda=lambdap,Q=Qp,sigma0=sigmap,options=imputeEMFunOptions)[["imp"]]
        
        #impute skipped values:
        yp[is.na(yp)] = ysave[is.na(yp)]
        rm(ysave)
        #now update latent space:
        for (k in 1:dim(yp)[3]) 
        {
          y[,1:Npc,k] = yp[,,k]%*%t(trans)
        }
      }
      else
      {
        ysave = y
        y = ImputeEMFun(y=y,x=x,where=where,dt=dt,Q=Q,W=W,lambda=lambda,options=imputeEMFunOptions)[["imp"]]
        y[is.na(y)] = ysave[is.na(y)]
        rm(ysave)
      }
    }    
    
    #(b) estimate
    fit = FitFun(y=y,x=x,dt=dt,W=W,lambda=lambda,Q=Q,calcNoise=calcNoise,options=fitOptions)
    W = fit[["W"]]
    lambda = fit[["lambda"]]
    Q = fit[["Q"]]
    sigma = solve(fit[["Q"]])
  }
  
  
  if(verbose)
  {
    print("NAs after iterating:")
    print(sum(is.na(y)))
  }

  #estimate fitting errors:
  if(estimateErrors)
  {
    efit = ErrorFun(y=y,x=x,dt=dt,W=W,lambda=lambda,Q=Q,options=errorFunOptions)
    dW = efit[["dW"]]
    dlambda = efit[["dlambda"]]
    dQ = efit[["dQ"]]
    dsigma = efit[["dsigma"]]
  }
  else #update, was missing (will crash without) #added June 2024
  {
    dW = NA*W
    dlambda = NA*lambda
    dQ = NA*Q
    dsigma = NA*sigma
  }
  
  
  #compute determinant of sigma in latent space if reduced rank
  #do this instead of using loglikz because loglikz has imputed values
  if(is.null(Npc)) logdetQ = NULL
  else if(Npc < ncol(y0))
  {
    warning("Q is reduced rank, taking determinant in latent space (this is fine so long as the transformation is unitary e.g. PCA)...")
    logdetQ = determinant(Q,T)
    logdetQ = as.numeric(logdetQ$sign*logdetQ$modulus)    
  }
  else logdetQ = NULL
  
  
  loglikz = tryCatch(loglik.SF(y=y,x=x,W=W,lambda=lambda,Q=Q,dt=dt,logdetQ=logdetQ),error=function(e) {return(NA)})
  
  
  #undo transformations
  zimp = y
  Qz = Q
  sigmaz = sigma
  lambdaz = lambda
  dlambdaz = dlambda
  Wz = W
  dWz = dW
  sf = list(W=Wz,lambda=lambdaz,trans=trans,transinv=transinv)
  class(sf)="SF"
  prz = predict.SF(sf,y=y,x=x,dt=dt)
  resz = y[,,-1,drop=F]-prz
  if(doPCA)
  {
    W = transinv[,1:Npc,drop=F]%*%W%*%trans[1:Npc,,drop=F]

    if(!is.null(dW))  dW = sqrt(transinv[,1:Npc,drop=F]^2%*%dW^2%*%trans[1:Npc,,drop=F]^2)
    
    yt = y0
    for (k in 1:dim(yt)[3]) yt[,,k] = y[,,k]%*%t(transinv)[1:Npc,,drop=F]
    y = yt
    
    lambda = transinv[,1:Npc,drop=F]%*%lambda
    if(!is.null(dlambda)) dlambda = sqrt(transinv[,1:Npc,drop=F]^2%*%dlambda^2)
    
    Q = transinv[,1:Npc,drop=F]%*%Q%*%trans[1:Npc,,drop=F]
    sigma = transinv[,1:Npc,drop=F]%*%sigma%*%trans[1:Npc,,drop=F]
  }
  
  
  if(scale)
  {
    for (j in 1:ncol(y)) y[,j,] = y[,j,]*sdpp[j]
  }  
  if(center)
  {
    for (j in 1:ncol(y)) y[,j,] = y[,j,] + mupp[j]
    lambda[,"mu0"] = lambda[,"mu0"] + mupp #math works out like this
  }
  
  #rownames get dropped for some reason
  rownames(lambda) = colnames(y)
  
  
  #post amble
  
  #estimate training error
  sf = list(W=W,lambda=lambda,trans=trans,transinv=transinv)
  class(sf)="SF"
  pr = predict.SF(sf,y=y,x=x,dt=dt)
  res = y0[,,-1,drop=F]-pr
  
  #compute accuracy metrics for fit error
  #NULL model should be y_n+1 - y_n = const
  yc = y0
  for (j in 1:ncol(y0)) yc[,j,] = yc[,j,] - mean(yc[,j,],na.rm=T)
  r2 = 1-mean((y0[,,-1,drop=F]-pr)^2,na.rm=T)/mean(yc^2,na.rm=T)
  
  r2_vec = numeric(ncol(y0))
  for (j in 1:ncol(y0)) r2_vec[j] = R2(c(y0[,j,-1]),c(pr[,j,]),na.rm=T)
  names(r2_vec) = colnames(y0)
  
  mse_vec = numeric(ncol(y0))
  for (j in 1:ncol(y0)) mse_vec[j] = mean(c((y0[,j,-1]-pr[,j,])^2),na.rm=T)
  names(mse_vec) = colnames(y0)
  
  mse = mean(c((y0[,,-1]-pr[,,])^2),na.rm=T)
  rmse = sqrt(mse)
  mae = mean(c(abs(y0[,,-1]-pr[,,])),na.rm=T)
  
  #carry forward
  mse_carry = mean(c((y0[,,-1]-y0[,,-dim(y0)[3]])^2),na.rm=T)
  rmse_carry = sqrt(mse_carry)
  mae_carry = mean(c(abs(y0[,,-1]-y0[,,-dim(y0)[3]])),na.rm=T)
  
  #imputed carry forward
  mse_carry_imp = mean(c((y0[,,-1]-y[,,-dim(y0)[3]])^2),na.rm=T)
  rmse_carry_imp = sqrt(mse_carry_imp)
  mae_carry_imp = mean(c(abs(y0[,,-1]-y[,,-dim(y0)[3]])),na.rm=T)
  
  e = tryCatch(eigen(W),error=function(e){return(NA)})
  
  #will fail if Q is reduced rank
  loglik = tryCatch(loglik.SF(y=y0,x=x,W=W,lambda=lambda,Q=Q,dt=dt,logdetQ=logdetQ),error=function(e) {return(NA)})
  
  #fix names
  colnames(W) = colnames(y)
  rownames(W) = colnames(y)
  colnames(lambda) = colnames(x)
  rownames(lambda) = colnames(y)
  
  oneoverWz = 1/Wz #inverse if Wz is diagonal
  oneoverdWz = dWz/Wz^2

  
  #notes for return values:
  #y: what you provided
  #yimp: what you provided with new, imputed values
  #x0: what you provided
  #x: ????
  #xpp: what you provided with preprocessing, including imputing the age variable and adding an intercept (optionally)
  results = list(lambda=lambda,W=W,Q=Q,sigma=sigma,dlambda=dlambda,dW=dW,dQ=dQ,dsigma=dsigma,eigenvalues=e[[1]],
                 y=y0,x=x,dt=dt,yimp=y,x0=x00,x=x0,xpp=x,where=where,age = age, s = s,iter=iter,
                 mupp=mupp,sdpp=sdpp,center=center,scale=scale,
                 loglik=loglik,r2=r2,mse=mse,rmse=rmse,mae=mae,r2_vec=r2_vec,mse_vec=mse_vec,
                 mse_carry=mse_carry,rmse_carry=rmse_carry,mae_carry=mae_carry,
                 mse_carry_imp=mse_carry_imp,rmse_carry_imp=rmse_carry_imp,mae_carry_imp=mae_carry_imp,
                 Npc=Npc,trans=trans,transinv=transinv,
                 addIntercept=addIntercept,
                 FitFun=FitFun,
                 fitOptions = fitOptions,
                 ErrorFun=ErrorFun,
                 errorFunOptions = errorFunOptions,
                 estimateErrors=T,
                 ImputeFirstFun = ImputeFirstFun, #for initial imputation
                 imputeFirstFunOptions = imputeFirstFunOptions,
                 imputeFirst=imputeFirst, #impute intrawave-conditional mean as first step? i.e. before computing anything else                 
                 ImputeEMFun = ImputeEMFun, #for expectation-maximization (iterative) imputation
                 imputeEMFunOptions = imputeEMFunOptions,
                 imputeMean=imputeMean, #impute conditional mean as we iterate through? #has this been validated?
                 imputeAge=imputeAge, #impute age columns of x using dt?
                 imputeAgeCols = imputeAgeCols,
                 imputeDead=imputeDead, #controls whether or not you impute dead people
                 imputeCensored=imputeCensored, #controls whether or not you impute censored people
                 age = age, #needed for informative censoring
                 s = s, #survival data, include if you want to use informative censoring
                 where = where, #where to impute (if at all)
                 Npc = Npc, #number of PCs to use for fitting; NULL to avoid PCA entirely
                 imputeInPCSpace = imputeInPCSpace, #do you want to perform imputation in latent space? May imputed too many values
                 pcaIndex=pcaIndex, #time index to use for PCA
                 center=center, #pre-processing step, center by mean
                 scale=scale, #pre-processing step, scale by sd
                 addIntercept=addIntercept, #add column of 1 to x, first index
                 fixedQ = fixedQ, #set to NA if you want to estiamte
                 doNotImpute=doNotImpute
  )
  class(results)="sf"
  return(results)
  
}


################################# bootstrap function  ######################################
################################# bootstrap function  ######################################
################################# bootstrap function  ######################################

BootSF = function(y, #data... 3d array individuals x var x times
                  x, #covariate data
                  dt, #time data
                  s = NULL, #survival data, for censorship options
                  age = NULL, #needed for censorship options
                  ygt=y, #ygt is used for accuracy measures, generally should be different for imputed data
                  xgt=x, #should match ygt
                  dtgt=dt, #should match ygt
                  where = is.na(y), #where you want to impute
                  nboot=100,
                  mc.cores=1, #number of cores to use (for parallelization)
                  maxRefit=2, #maximum times you'll try to fit again if a fold fails
                  plot=F,
                  options=list(), #sent to FitSF
                  weights=NULL, #use to pick certain individuals more often (seems to work for rare outcomes)
                  quantiles = NULL, #this are very slow
                  saveBoot=F, #keep bootstrapped samples?
                  storePar=c("W","lambda","sigma_noise_scale","eig_re","diagW","diagW_sort"), #parameters to save after bootstrapping
                  keepImp=T, #save imputed values?
                  ygtimp=NULL, #you can supply pre-imputed ygt (will not affect GT, just prediction)
                  imputePredict=T #impute ygt before making prediction (will not affect GT) #don't change unless you're sure
                  )
{

  if(mc.cores > 1) library(parallel) #needed for parallel computing (surprise, surprise)

  if(any(grepl("imputeDead",names(options))) | any(grepl("imputeCensored",names(options))))
  {
    if(is.null(s))
    {
      stop("You must provide s if you want to apply censorship options for imputation")
    }
    if(is.null(age))
    {
      stop("You must provide age if you want to apply censorship options for imputation")
    }
  }
  
  R2bs = function(y,pr)
  {
    residual = y-pr
    yc = y
    for (i in 1:dim(y)[2]) for (j in 1:dim(y)[3]) yc[,i,j] = y[,i,j]-mean(y[,i,j],na.rm=T)
    
    r2  = 1 - mean(residual^2,na.rm=T)/mean(yc^2,na.rm=T)
    return(r2)
  }
  
  wrapper = function(ind)
  {
    #library(Matrix) # I don't think I need this
    l = list()
    l$y  = y[ind,,,drop=F]
    l$x  = x[ind,,,drop=F]
    l$dt = dt[ind,,drop=F]
    l$s = s[ind,,drop=F]
    l$where = where[ind,,,drop=F]
    l$age = age[ind,,drop=F]
    l = c(l,options)

    fit = tryCatch(FitSF.wrapper(l),error=function(e){return(NULL)})
    
    
    #fit = tryCatch(FitSF.wrapper(l),error=function(e){return(e)}) #for debugging
    #print(fit)
    if(is.null(fit)) 
    {
      warning("Fit failed, setting to NA (check column names)")
      return(NA) #fit failed
    }
    res = fit[c("mu0","lambda","W","sigma","Q","loglik","r2","mse","rmse","mae",
                "r2_vec","mse_vec","dW","dlambda"
                )]
    
    #imputed values
    if(keepImp) 
    {
        yimpsum = array(NA,dim=dim(y)) #sum #changed to NA may 2023
        dimnames(yimpsum) = dimnames(y)
        y2impsum = yimpsum #^2
        #Nimp = rep(0,nrow(y)) 
        Nimp = yimpsum
        inds = unique(ind)
        for (ii in 1:length(inds)) 
        {
          logi = ind == inds[ii]
          if(any(is.na(logi)))
          {
            print(logi)
            stop("NAs in logi")
          }
          yimpsum[inds[ii],,] = apply(fit$yimp[logi,,,drop=F],2:3,sum,na.rm=T)
          y2impsum[inds[ii],,] = apply(fit$yimp[logi,,,drop=F]^2,2:3,sum,na.rm=T) #update, May 2023 - for variance
          #sum(logi) is number of times you appeared, we need to consider how often you were IMPUTED
          #Nimp[inds[ii]] = sum(logi)
          Nimp[inds[ii],,] = apply(!is.na(fit$yimp[logi,,,drop=F]),2:3,sum)
          
        }
        res[["yimpsum"]] = yimpsum
        res[["y2impsum"]] = y2impsum
        res[["Nimp"]] = Nimp
        rm(yimpsum)
        rm(y2impsum)
    }
    

    W = res$W
    res[["diagW"]] = diag(res$W)
    res[["diagW_sort"]] = sort(diag(res$W),decreasing=T)
    res[["diagW_sort_space"]] = diff(sort(diag(res$W),decreasing=T))
    eig = eigen(res$W)
    #sort eigenvalues 
    l = sort.list(eig[[1]],decreasing=T)
    eig[[1]] = eig[[1]][l]
    eig[[2]] = eig[[2]][,l]
    
    res[["eig_re"]] = Re(eig[[1]])
    res[["eig_im"]] = Im(eig[[1]])
    res[["diag_min_eig"]] = Re(sort(diag(res$W),decreasing=T) - Re(eig[[1]])) #imaginary part is just -Im(eig[[1]])
    res[["eig_space_re"]] = Re(diff(eig[[1]]))
    res[["eig_space_im"]] = Im(diff(eig[[1]]))
    res[["timescale"]] = -1/Re(eig[[1]])
    
    res[["drift"]] = -res[["W"]]%*%res[["lambda"]]
    
    #out-of-sample prediction
    ytrain = ygt[ind,,,drop=F]
    ytest = ygt[-ind,,,drop=F]
    yraw = y[-ind,,,drop=F]
    if(is.null(xgt))
    {
      if("mu0" %in% colnames(fit$xpp))
      {
        xtrain = array(1,dim=c(nrow(ytrain),1,dim(ytrain)[3]))
        colnames(xtrain) = c("mu0") 
        
        xtest = array(1,dim=c(nrow(ytest),1,dim(ytest)[3]))
        colnames(xtest) = c("mu0") 
        
        xraw = array(1,dim=c(nrow(yraw),1,dim(yraw)[3]))
        colnames(xraw) = c("mu0") 
      }
      else 
      {
        xtrain = NULL
        xtest = NULL
        xraw = NULL
      }
    }
    else if("mu0" %in% colnames(fit$xpp))
    {
      #add 1s if missing
      if(!any(grepl("mu0",colnames(xgt))))
      {
        xtrain = xgt[ind,c(1,1:ncol(xgt)),,drop=F]
        xtrain[,1,] = 1
        colnames(xtrain) = c("mu0",colnames(xgt)) 
        
        xtest = xgt[-ind,c(1,1:ncol(xgt)),,drop=F]
        xtest[,1,] = 1
        colnames(xtest) = c("mu0",colnames(xgt)) 
        
        xraw = x[-ind,c(1,1:ncol(x)),,drop=F]
        xraw[,1,] = 1
        colnames(xraw) = c("mu0",colnames(x)) 
      }
      else #mu0 already in x
      {
        xtrain = xgt[ind,,,drop=F] 
        xtest = xgt[-ind,,,drop=F] 
        xraw = x[-ind,,,drop=F] 
      }
    }
    else    
    {
      xtrain = xgt[ind,,,drop=F]
      xtest = xgt[-ind,,,drop=F]
    }
    dttrain = dtgt[ind,,drop=F]
    dttest = dtgt[-ind,,drop=F]
    dtraw = dt[-ind,,drop=F]
    
    
    if(is.null(ygtimp))
    {
      ytestimp = ytest
      ytrainimp = ytrain
    }
    else
    {
      ytrainimp = ygtimp[ind,,,drop=F]
      ytestimp = ygtimp[-ind,,,drop=F]
    }
    #do you want to impute the test/train data before prediction?
      #typically you want to since otherwise NAs will get propagated and you won't have many GT valuse
      #you will never compare against imputed values, this is purely to get the prediction
    if(imputePredict) #seems like we should want to do this...
    {
      #train data
      wheretrain = is.na(ytrain)
      ytrainimp = fit$ImputeFirstFun(y=ytrain,x=xtrain,dt=dttrain,where=wheretrain,options=fit$imputeFirstFunOptions)[["imp"]]
      if(fit$imputeMean)
      {
        ysave = ytrainimp
        ytrainimp = fit$ImputeEMFun(y=ytrainimp,x=xtrain,dt=dttrain,where=wheretrain,W=fit$W,lambda=fit$lambda,
                                    Q=fit$Q,options=fit$imputeEMFunOptions)[["imp"]]
        ytrainimp[is.na(ytrainimp)] = ysave[is.na(ytrainimp)]
        rm(ysave)
      }

      #test data
      wheretest = is.na(ytest)
      ytestimp = fit$ImputeFirstFun(y=ytest,x=xtest,dt=dttest,where=wheretest,options=fit$imputeFirstFunOptions)[["imp"]]
      if(fit$imputeMean)
      {
        ysave = ytestimp
        ytestimp = fit$ImputeEMFun(y=ytestimp,x=xtest,dt=dttest,where=wheretest,W=fit$W,lambda=fit$lambda,
                                   Q=fit$Q,options=fit$imputeEMFunOptions)[["imp"]]
        ytestimp[is.na(ytestimp)] = ysave[is.na(ytestimp)]
        rm(ysave)
      }
      
    }
    
    prtest = predict.SF(fit,y=ytestimp,x=xtest,dt=dttest)
    prtrain = predict.SF(fit,y=ytrainimp,x=xtrain,dt=dttrain)
    
    #'raw' accuracy measures from y rather than ygt
    prraw = predict.SF(fit,y=yraw,x=xraw,dt=dtraw)
    residual_raw =  yraw[,,-1,drop=F] - prraw
    mse_raw = mean(residual_raw^2,na.rm=T)
    rmse_raw = sqrt(mse_raw)
    mae_raw = mean(abs(residual_raw),na.rm=T)
    r2_raw = R2bs(yraw[,,-1,drop=F],prraw)
    res[["mse_raw_test"]] = mse_raw
    res[["rmse_raw_test"]] = rmse_raw
    res[["mae_raw_test"]] = mae_raw
    res[["r2_raw_test"]] = r2_raw
    
    
    #train accuracy measures
    residual_train = ytrain[,,-1,drop=F] - prtrain
    mse = mean(residual_train^2,na.rm=T)
    rmse = sqrt(mse)
    mae = mean(abs(residual_train),na.rm=T)
    r2 = R2bs(ytrain[,,-1,drop=F],prtrain)
    res[["mse_train"]] = mse
    res[["rmse_train"]] = rmse
    res[["mae_train"]] = mae
    res[["r2_train"]] = r2
    
    #out-of-sample (test) accuracy measures
    residual_test = ytest[,,-1,drop=F] - prtest
    mse = mean(residual_test^2,na.rm=T)
    rmse = sqrt(mse)
    mae = mean(abs(residual_test),na.rm=T)
    r2 = R2bs(ytest[,,-1,drop=F],prtest)
    res[["mse_test"]] = mse
    res[["rmse_test"]] = rmse
    res[["mae_test"]] = mae
    res[["r2_test"]] = r2
    
    #per-variable metrics
    mse_vec_train = numeric(ncol(ytrain))
    names(mse_vec_train) = colnames(y)
    rmse_vec_train = mse_vec_train
    mae_vec_train = mse_vec_train
    r2_vec_train = mse_vec_train
    mse_vec_test = numeric(ncol(ytest))
    names(mse_vec_test) = colnames(y)
    rmse_vec_test = mse_vec_test
    mae_vec_test = mse_vec_test
    r2_vec_test = mse_vec_test
    for (j in 1:ncol(ytest))
    {
      mse_vec_train[j] = mean(residual_train[,j,]^2,na.rm=T)
      rmse_vec_train[j] = sqrt(mse_vec_train[j])
      mae_vec_train[j] = mean(abs(residual_train[,j,]),na.rm=T)
      r2_vec_train[j] = R2bs(ytrain[,j,-1,drop=F],prtrain[,j,,drop=F])
      
      mse_vec_test[j] = mean(residual_test[,j,]^2,na.rm=T)
      rmse_vec_test[j] = sqrt(mse_vec_test[j])
      mae_vec_test[j] = mean(abs(residual_test[,j,]),na.rm=T)
      r2_vec_test[j] = R2bs(ytest[,j,-1,drop=F],prtest[,j,,drop=F])
    }
    res[["mse_vec_train"]] = mse_vec_train
    res[["rmse_vec_train"]] = rmse_vec_train
    res[["mae_vec_train"]] = mae_vec_train
    res[["r2_vec_train"]] = r2_vec_train
    res[["mse_vec_test"]] = mse_vec_test
    res[["rmse_vec_test"]] = rmse_vec_test
    res[["mae_vec_test"]] = mae_vec_test
    res[["r2_vec_test"]] = r2_vec_test
    #if model was just y_n+1 = y_n + mu + eps
    res[["null_sd_train"]] = apply(ytrain[,,-1,drop=F]-ytrain[,,-dim(ytrain)[3],drop=F],2,sd,na.rm=T)
    res[["null_sd_test"]] = apply(ytest[,,-1,drop=F]-ytest[,,-dim(ytest)[3],drop=F],2,sd,na.rm=T)
    
    
    #out-of-sample log-likelihood
    Q = tryCatch(solve(fit$sigma), error=function(e){return(NA)})
    if(any(is.na(Q))) 
    {
      warning("could not invert sigma, skipping loglik")
      lltrain = NA
      lltest = NA
    }
    else
    {
      logdetQ = NULL
      #can provide prediction so we don't have to recalculate
      lltrain = tryCatch(loglik.SF(y=ytrain,x=xtrain,W=fit$W,lambda=fit$lambda,mu0=fit$mu0,Q=Q,dt=dttrain,logdetQ=logdetQ,pr0=prtrain),
                           error=function(e){return(NA)})
      if(is.na(lltrain)) warning("log-likelihood failed (train)")
      lltest = tryCatch(loglik.SF(y=ytest,x=xtest,W=fit$W,lambda=fit$lambda,mu0=fit$mu0,Q=Q,dt=dttest,logdetQ=logdetQ,pr0=prtest),
                          error=function(e){return(NA)})
      if(is.na(lltest)) warning("log-likelihood failed (test)")
    }
    res[["loglik_train"]] = lltrain
    res[["loglik_test"]] = lltest
    
    return(res)
  }
  
  l = list()
  for (i in 1:nboot) l[[i]] = sample(1:nrow(y),nrow(y),replace=T,prob=weights)
  
  if(mc.cores > 1)
  {
    cl = makePSOCKcluster(mc.cores)
    
    clusterExport(cl,c( "wrapper","y","x","dt","s","ygt","xgt","dtgt","options","where","age",
                        "FitSF","ImputeMean.SF","Impute","ImputeGaussian","predict.SF","FitSF.wrapper",
                        "FitOP.SF","FitLM.SF",
                        "AymErrorMoments.SF","FitMoments.SF",
                        "FitW.SF","FitLambda.SF","FitLambdaEq.SF","FitQ.SF",
                        "loglik.SF","imputePredict",
                        "prcomp2","cov2","R2bs",
                        "keepImp"
    ), envir=environment())
    
    results = parLapply(cl,l,wrapper)
    
    stopCluster(cl)
  }
  else
  {
    results = list()
    for (i in 1:length(l))
    {
      results[[i]] = wrapper(l[[i]])
    }
  }
  
  #check for failed fits
  drop=numeric()
  for (i in 1:length(results))
  {
    
    if(all(is.na(results[[i]])))
    {
      print(sprintf("boostrap %d failed",i))
      if(maxRefit > 0)
      {
        for (ii in 1:maxRefit)
        {
          print("trying again...")
          results[[i]] = wrapper(sample(1:nrow(y),nrow(y),replace=T))
          if(!all(is.na(results[[i]]))) break
          else print(sprintf("boostrap %d failed again!",i))
        }
        if(all(is.na(results[[i]]))) #just drop
        {
          print(sprintf("dropping sample %d",i))
          drop = c(drop,i)
        }
      }
      else drop = c(drop,i)
    }
  }
  failures = drop
  
  if(length(drop)>0) 
  {
    results = results[-drop]
    print(sprintf("total bootstraps after dropping: %d",length(results)))
  }
  if(length(results) < 1)
  {
    warning("all fits failed!")
    return(NA)
  }
  
  #post-process: rearrange list
  res = results[[1]]
  results2 = list()
  for (j in 1:length(res)) results2[[j]] = list()
  
  for (i in 1:length(results))
  {
    #print(i)
    for (j in 1:length(res))
    {
      #print(j)
      results2[[j]][[i]] = results[[i]][[j]]
    }
  }
  names(results2) = names(res)
  
  
  
  musd = list()
  par = list()
  print("musd")
  for (j in 1:length(res)) 
  {
    #print(names(res)[j])
    #for liar?
    #if(names(results2)[j]%in%c("lambda")) musd[[j]] = results2[[j]]
    #else musd[[j]] = ListMeanSD(results2[[j]],na.rm=T)
    musd[[j]] = ListMeanSD(results2[[j]],na.rm=T)
    if(names(results2)[j]%in%storePar)
    {
      par[[names(results2)[j]]] = do.call(rbind,lapply(results2[[j]],c))
      if(is.null(dim(results2[[j]][[1]]))) nms = names(results2[[j]][[1]])
      else nms = c(outer(rownames(results2[[j]][[1]]),colnames(results2[[j]][[1]]),paste,sep="_x_"))
      colnames(par[[names(results2)[j]]]) = nms
    }
  }
  names(musd) = names(results2)
  
  print("quants")
  quants = list()
  if(!is.null(quantiles))
  {
    for (k in 1:length(res))
    {
      quants[[k]] = list()
      for (j in 1:length(quantiles)) 
      {
        quants[[k]][[j]] = ListDo(results2[[k]],fun=quantile,probs=quantiles[j],na.rm=T)
      }
      names(quants[[k]]) = sprintf("q%.0f",100*quantiles)
    }
    names(quants) = names(results2)
  }
  
  
  #yimp
  #we can always change a sum to an average by diving top and bottom by N
  #so we have < sum (y) >/ < N > = <y> (appropriately weighted)
  print("imputed values...")
  if(keepImp)
  {
      yimp = musd$yimpsum[[1]]
      #yimpsd = musd$yimpsum[[2]] #not right (see above note)
      yimpsd = array(NA,dim=dim(yimp))
      names(yimpsd) = names(yimp)
      y2imp = musd$y2impsum[[1]]
      for (i in 1:nrow(yimp)) #normalize (average) sums by (average) number in each
      {
        yimp[i,,] = yimp[i,,]/musd$Nimp[[1]][i,,]
        yimp[i,,][musd$Nimp[[1]][i,,] == 0] = NA #technically Nimp is number of obs + imps
        
        #sd = sqrt(<x^2> + <x>^2)
        #biased estimator because I don't know correct N (this is ratio of means so don't subtract 1)
        #add small number to prevent NaN
        yimpsd[i,,] = sqrt(y2imp[i,,]/(musd$Nimp[[1]][i,,]) - yimp[i,,]^2 + 1e-10)
      }
  }
  else
  {
    yimp = NULL
    yimpsd = NULL
  }
  
  #musd: mean and sd
  #quants: specified quantiles
  final_results = list(meansd=musd,quants=quants,nboot0=nboot,nboot=length(results),
                       options=options,par=par,s=s,age=age,
                       yimp=yimp,yimpsd=yimpsd,inds=l,failed=failures)
  if(saveBoot) final_results$boot = results
  
  #options:
  final_results$ygt = ygt
  final_results$xgt = xgt
  final_results$dtgt = dtgt
  final_results$where  = where
  final_results$mc.cores = mc.cores
  final_results$maxRefit = maxRefit
  final_results$weights = weights
  final_results$quantiles = quantiles
  final_results$saveBoot = saveBoot
  
  return(final_results)
}

################################# low-level functions ######################################
################################# low-level functions ######################################
################################# low-level functions ######################################

Impute = function(y,
                  x, #covariates e.g. sex and age, for predicting mean #should be same size as y if interpx=T or dim(x)[3] = dim(y)[3]-1 
                  dt,
                  where=is.na(y),
                  options=list()
                  )
{
  #methods:
    #ind
      #imputes individual mean for each variable (if possible, otherwise NA)
    #carryback
      #carrys forward previous values for each variable and individual, then carrys backward future values into any still NA (if possible, otherwise NA)
    #forwardcarry ('carry')
      #carrys forward previous values for each variable and individual
    #both
      #first applies 'ind' then imputes the multivariate Gaussian mean
    #multimean
      #imputes conditional Gaussian mean assuming time-invariant stats
    #popmean
      #imputes population mean
  #default options
  defaults = list(
              method="carryback" #default: carry forward previous values then carry back future values
            )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  firstImputeMethod = options[["method"]]
  firstImputeMethod = tolower(firstImputeMethod)
  
  #unconditioned mean: #looks like shit, but  makes sense for equilibrium dynamics
  if(grepl("ind",firstImputeMethod))
  {
    #used to be flatmean, but I think that's wrong so moving that to flatmean
    #imputed individual mean
    #in equilibrium mean is a stationary statistic
    for (j in 1:ncol(y))
    {
      impmu = apply(y[,j,],1,mean,na.rm=T)
      for (k in 1:dim(y)[3])
      {
        naLogi = where[,j,k]
        y[naLogi,j,k] = impmu[naLogi]
      }
    }
  }
  else if(grepl("carry",firstImputeMethod) & grepl("back",firstImputeMethod))
  {
    #carry forward + back fill

    #carry forward values
    for (j in 1:ncol(y))
    {
      for (k in 2:dim(y)[3])
      {
        naLogi = where[,j,k]
        y[naLogi,j,k] = y[naLogi,j,k-1]
      }
    }
    
    #backfill
    for (j in 1:ncol(y))
    {
      for (k in dim(y)[3]:2)
      {
        naLogi = where[,j,k-1] & is.na(y[,j,k-1])
        y[naLogi,j,k-1] = y[naLogi,j,k]
      }
    }
  }
  else if(grepl("carry",firstImputeMethod))
  {
    #carry forward values
    for (j in 1:ncol(y))
    {
      for (k in 2:dim(y)[3])
      {
        naLogi = where[,j,k]
        y[naLogi,j,k] = y[naLogi,j,k-1]
      }
    }
  }
  else if (firstImputeMethod=="bothv2") #carryback + conditional mean ("cond")
  {
    #carry forward + back fill
    #carry forward values
    for (j in 1:ncol(y))
    {
      for (k in 2:dim(y)[3])
      {
        naLogi = where[,j,k]
        y[naLogi,j,k] = y[naLogi,j,k-1]
      }
    }
    
    #backfill
    for (j in 1:ncol(y))
    {
      for (k in dim(y)[3]:2)
      {
        naLogi = where[,j,k-1] & is.na(y[,j,k-1])
        y[naLogi,j,k-1] = y[naLogi,j,k]
      }
    }
    
    #now fill anything still missing using Gaussian stats
    for (k in 1:dim(y)[3])
    {
      yMat = matrix(y[,,k],nrow=nrow(y),ncol=ncol(y))
      whereMat = matrix(is.na(y[,,k]) & where[,,k],nrow=nrow(y),ncol=ncol(y))
      y[,,k] = ImputeGaussian(yMat,where=whereMat)[["imp"]] 
    }
  }
  else if (firstImputeMethod=="both") #individual mean (ind/mean) + conditional mean ("cond")
  {

    for (j in 1:ncol(y))
    {
      impmu = apply(y[,j,],1,mean,na.rm=T)
      for (k in 1:dim(y)[3])
      {
        naLogi = where[,j,k]
        y[naLogi,j,k] = impmu[naLogi]
      }
    }
    
    for (k in 1:dim(y)[3])
    {
      yMat = matrix(y[,,k],nrow=nrow(y),ncol=ncol(y))
      whereMat = matrix(is.na(y[,,k]) & where[,,k],nrow=nrow(y),ncol=ncol(y))
      y[,,k] = ImputeGaussian(yMat,where=whereMat)[["imp"]] 
    }
    #print("NAs:")
    #print(sum(is.na(c(y))))
  }
  else if (firstImputeMethod=="multimean") #Gaussian statistics (flat/time-invariant)
  {
    for (k in 1:dim(y)[3])
    {
      yMat = matrix(y[,,k],nrow=nrow(y),ncol=ncol(y))
      whereMat = matrix(where[,,k],nrow=nrow(y),ncol=ncol(y))
      y[,,k] = ImputeGaussian(yMat,where=whereMat)[["imp"]] 
    }
  }
  else if (firstImputeMethod=="popmean") #population mean
  {
    #assume all variable j at timepoint k are same
    for(k in 1:dim(y)[3]) for(j in 1:ncol(y))
    {
      naLogi = where[,j,k]
      y[naLogi,j,k] = mean(y[,j,k],na.rm=T)
    }
  }
  else    
  {
    stop("Imputation method not found")
  }
  l = list(imp=y,where=where,m=1,options=options)
  return(l)
}

ImputeGaussian = function(x, #data matrix or dataframe
                      where = is.na(x),
                      mu= NULL, #mean of process that generated x
                      sigma = NULL #covariance 
)
{
  #general Gaussian imputation of conditional mean
  if(is.null(mu)) mu = apply(x,2,mean,na.rm=T)
  if(is.null(sigma)) sigma = cov(x,use="pairwise.complete")
  
  xi = x
  for (i in 1:nrow(xi))
  {
    miss = where[i,]
    nmiss = sum(miss)
    if(nmiss > 0)
    {
      obs = !miss
      nobs = sum(obs)
      if(nobs < 1) #no observations
      {
        muimp = mu[miss]
      }
      else 
      {
        
        sigma_mo = sigma[miss,obs,drop=F]
        sigma_oo = sigma[obs,obs,drop=F]
        Q_oo = tryCatch(solve(sigma_oo),error=function(e){return(NA)})
        if(any(is.na(Q_oo)))
        {
          warning("Could not invert sigma_oo, setting to Q=NA...")     
          Q_oo = matrix(NA,nrow=nrow(sigma_oo),ncol=ncol(sigma_oo))
        }
        
        muimp = mu[miss] + sigma_mo%*%Q_oo%*%(x[i,obs]-mu[obs])
      }
      xi[i,miss] = muimp
    }
    
  }
  return(list(imp=xi,where=where,mu=mu,sigma=sigma,m=1))
}


ImputeMean.SF = function(y,x,dt,
                         where=is.na(y),
                         W,
                         lambda,
                         Q=diag(1,ncol(y)),
                         sigma0=NULL, 
                         options=list()
                         )
{
  #default options
  defaults = list(
    imputeFirst=T, #impute the first timepoint? can be wonky since it inverts W (which can be close to 0)
    firstImputeMethod="carry" #how do you want to impute the first time index? options: carry back next value or backpropagate dynamics in reverse
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  if(is.null(sigma0) & is.null(Q)) stop("need sigma0 or Q")
  else if (is.null(sigma0)) 
  {
    sigma0 = solve(Q)
  }
  
  
  W0 = W #will replace W during iterations
  Id = diag(1,nrow=nrow(W0),ncol=ncol(W0))
  
  ymu= apply(y,2:3,mean,na.rm=T) #i'll use this when I don't know next timepoint

  
  #indices that need to be imputed
  toImp = list()
  for (k in 1:dim(where)[3]) toImp[[k]] = which(apply(where[,,k,drop=F],1,any))
  #useful for using matrix product to match rows
  whereint = 2*where - 1 #T: 1, #F: -1
  
  firstImputeMethod = options[["firstImputeMethod"]]
  imputeFirst = options[["imputeFirst"]]
  
  #update: you have to loop over k THEN nrows now because you're imputing by missingness type, not by row
  #otherwise you'll have missing values in prev timepoint when you try to impute
  yc = y
  k = 1
  for(i in 1:nrow(y))
  {

    #impute first timepoint #validated
    if(imputeFirst) 
    {
      miss = where[i,,k]
      nmiss = sum(miss)
      obsp1 = !is.na(y[i,,k+1]) 
      nobsp1 = sum(obsp1)
      
      if(grepl("carry",tolower(firstImputeMethod)))
      {
        ynp1 = ymu[,k+1]
        ynp1[obsp1] = yc[i,obsp1,k+1]
        muimp =  ynp1
        yc[i,miss,k] = muimp[miss] #was missing
      }
      else if(nmiss > 0 & nobsp1 > 0)
      {
        xvec = x[i,,k]
        W = diag(1,nrow(W0)) + dt[i,k]*W0
        mu = lambda%*%xvec
        delta = -dt[i,k]*as.numeric(W0%*%mu)
        
        obs = !miss
        nobs = sum(obs)
        Winv_full = solve(W)
        Winv = Winv_full[miss,obsp1,drop=F]

        if(nobs > 0) #ad hoc
        {
          sigma = Winv_full%*%sigma0%*%t(Winv_full)*abs(dt[i,1]) #include dt anyways because why not
          sigma_mo = sigma[miss,obs,drop=F]
          sigma_oo = sigma[obs,obs,drop=F]
          Q_oo = tryCatch(solve(sigma_oo),error=function(e){return(NA)})
          if(any(is.na(Q_oo)))
          {
            warning("Could not invert sigma_oo, setting to Q=NA...")     
            Q_oo = matrix(NA,nrow=nrow(sigma_oo),ncol=ncol(sigma_oo))
          }
          
          ynp1 = ymu[,k+1]
          ynp1[obsp1] = yc[i,obsp1,k+1]
          muimp = Winv_full%*%(ynp1-delta) 
          muimp = muimp[miss] + sigma_mo%*%Q_oo%*%(yc[i,obs,k]-muimp[obs]) 
        }
        else
        {
          ynp1 = ymu[,k]
          ynp1[obsp1] = yc[i,obsp1,k+1]
          muimp = Winv_full[miss,]%*%(ynp1 - delta)
        }
        muimp = as.numeric(muimp)
        if(length(muimp) != sum(miss)) stop("different number of imputed values than needed!")
        yc[i,miss,k] = muimp
        rm(muimp)
      }
      else if(nmiss > 0)  # ad hoc solution / guess #we know nothing about dynamics, so impute mean
      {
        muimp = ymu[miss,k] # NA
        yc[i,miss,k] = muimp
        rm(muimp)
      }
    }
    else break
  }
  
  #impute subsequent timepoints
  kstart = 2
  if(grepl("carry",tolower(firstImputeMethod))) kstart = 3 #we carried k=2 back
  for (k in kstart:dim(where)[3])
  {
    while(length(toImp[[k]])>0)
    {
      i = toImp[[k]][1]
      
      miss = where[i,,k]
      missint = 2*miss-1
      nmiss = sum(miss)
      
      missintm1 = 2*where[i,,k-1]-1
      
      #find everybody with same missingess pattern
      #must consider this timepoint and next
      wheremat = matrix(whereint[toImp[[k]],,k],nrow=length(toImp[[k]]),ncol=ncol(whereint))
      wherematm1 = matrix(whereint[toImp[[k]],,k-1],nrow=length(toImp[[k]]),ncol=ncol(whereint))
      impMeInds = which(((wheremat%*%missint)[,1] == ncol(whereint)) & ((wherematm1%*%missintm1)[,1] == ncol(whereint)) )
      impMe = toImp[[k]][impMeInds]
      
      dtind = outer(dt[impMe,k-1],rep(1,nmiss))  
      if(nmiss>0) # & nobsm1 > 0)
      {
        xmat = x[impMe,,k-1]
        mu = xmat%*%t(lambda)
        delta = -(mu%*%t(W0))
        for (j in 1:ncol(delta)) delta[,j] = delta[,j]*dt[impMe,k-1]
        
        obs = !miss
        nobs = sum(obs)

        if(nobs > 0) #not sure about this rule, ad hoc
        {
          dtindobs = outer(dt[impMe,k-1],rep(1,nobs)) 
          sigma_uo = sigma0[miss,obs,drop=F]
          sigma_oo = sigma0[obs,obs,drop=F]
          sigma_uu = sigma0[miss,miss,drop=F]
          Q_oo = tryCatch(solve(sigma_oo),error=function(e){return(NA)})
          if(any(is.na(Q_oo)))
          {
            warning("Could not invert sigma_oo, setting to Q=NA...")     
            Q_oo = matrix(NA,nrow=nrow(sigma_oo),ncol=ncol(sigma_oo))
          }
          transfer = sigma_uo%*%Q_oo  #dts cancel
          
          ycsubm1 = matrix(yc[impMe,,k-1],nrow=length(impMe),ncol(yc))
          ycsub = matrix(yc[impMe,obs,k,drop=F],nrow=length(impMe),nobs)
          muimp = ycsubm1%*%t(Id[miss,,drop=F])+(ycsubm1%*%t(W0[miss,,drop=F]))*dtind + delta[,miss,drop=F] + 
            ((ycsub-ycsubm1%*%t(Id[obs,,drop=F])-(ycsubm1%*%t(W0[obs,,drop=F]))*dtindobs - delta[,obs,drop=F])%*%t(transfer)) #dts cancel in transfer
          #sig = sigma_uu - sigma_uo%*%Q_oo%*%t(sigma_uo) #for multiple imputations, not used
          
        }
        else
        {
          ycsubm1 = matrix(yc[impMe,,k-1],nrow=length(impMe),ncol(yc))
          muimp = ycsubm1%*%t(Id[miss,,drop=F]) + ( ycsubm1%*%t(W0[miss,,drop=F]) )*dtind + delta[,miss,drop=F]
          #sig = sigma[miss,miss,drop=F] #for multiple imputations, not used
        }
        
        yc[impMe,miss,k] = muimp
        rm(muimp)
      }
      
      #drop imputed indices
      toImp[[k]] = toImp[[k]][-impMeInds]
    }
  }
  
  yi = yc
  
  l =list(imp=yi,where=where,m=1,options=options)
  return(l)
}


FitLM.SF = function(y, #3d array individuals x variables x times
                    x, #dim 3 be same size as y or smaller by one 
                    dt,
                    W = NULL,#stub
                    lambda = NULL,#stub
                    Q = NULL, #stub
                    calcNoise=F,
                    options = list()
                    )
{
  #noise must be diagonal
  #dt may be non-uniform so we need to used weighted regression
  #model is: y_n+1 - y_n = lambda*dt*y_n + beta*(dt*x_n) + eps
  #weights are 1/abs(dt)
  #beta_j = lambda*Lambda_j
  
  if(!is.null(options[["wOptions"]][["diagonal"]])) #another way of providing optional diagonal W constraint
  {
    options[["diagonalW"]] = options[["wOptions"]][["diagonal"]]
  }
  
  #default options
  defaults = list(
    diagonalW=T,
    FitFun = lm,
    error=F,
    noiseMethod="y",
    qOptions = list()
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]


  if(is.null(colnames(y))) colnames(y) = sprintf("y%d",1:ncol(y))
  if(is.null(colnames(x))) colnames(x) = sprintf("x%d",1:ncol(x))
  
  FitFun = options[["FitFun"]]
  if(options[["diagonalW"]])
  {
    W = diag(NA,nrow=ncol(y),ncol=ncol(y))
    rownames(W)=colnames(y)
    colnames(W) = rownames(W)
    dW = W
    lambda = matrix(NA,ncol=ncol(x),nrow=ncol(y))
    colnames(lambda) = colnames(x)
    rownames(lambda) = colnames(W)
    dlambda = lambda

    qterm = W
    dqterm = W
    weights = 1/c(abs(dt[,1:(dim(y)[3]-1)])) 

    for (j in 1:ncol(y))
    {
      #note: all predictors must be scaled by dt!
      #also need weights
      data = data.frame(yn=c(y[,j,-dim(y)[3]]*dt))
      data[,"dyn"] = c(y[,j,-1]-y[,j,-dim(y)[3]])
      for (jj in 1:ncol(x)) data[,colnames(x)[jj]] = c(x[,jj,1:(dim(y)[3]-1)]*dt)

      m = FitFun(dyn~.-1,data=data,weights=weights)
      
      #drop any asshole "`" which for some reason lm wants to put on
      C = coef(m)
      for (jj in 1:length(C)) names(C)[jj] = gsub("`","",names(C)[jj])

      W[j,j] = C["yn"]
      lambda[j,] = -C[colnames(lambda)]/W[j,j]
      pr = predict(m,newdata=data)

      if(options[["error"]])
      {
          #validation: lambda error bars are wonky
          #note: variables with NA coefficients will be droppped automatically
          #e.g. if you have a redundant dummy variable or two highly collinear predictors
          #also note: if you have weird characters in your name, e.g. "/", it will add "''" which will cause this to crash
          parSE = summary(m)$coefficients[,"Std. Error"]
          for (jj in 1:length(parSE)) names(parSE)[jj] = gsub("`","",names(parSE)[jj])
          
          if(!all(colnames(lambda)%in%names(parSE)))
          {
            print(rownames(summary(m)$coefficients))
            print(colnames(lambda))
            print(setdiff(colnames(lambda),rownames(summary(m)$coefficients)))
            print(summary(m)$coefficients)
            print(parSE)
            stop("Missing covariate in fit (too many predictors or naming error, check for '')")
          }
          dW[j,j] = parSE["yn"]
          d1 = parSE[colnames(lambda)]/W[j,j]
          d2 = C[colnames(lambda)]/W[j,j]^2*parSE["yn"]
          dlambda[j,] = sqrt(d1^2+d2^2)
        }
        else
        {
          dW = NULL
          dlambda = NULL
        }
    }
  }
  else
  {
    #same as diagonal you just let all previous y be predictors then you gotta invert stuff
    W = matrix(NA,nrow=ncol(y),ncol=ncol(y))
    rownames(W)=colnames(y)
    colnames(W) = rownames(W)
    dW=W
    lambda = matrix(NA,ncol=ncol(x),nrow=ncol(y))
    colnames(lambda) = colnames(x)
    rownames(lambda) = colnames(W)
    dlambda = lambda

    weights = 1/c(abs(dt[,1:(dim(y)[3]-1)])) #check this
    
    data = data.frame(dyn=c(y[,1,-1]-y[,1,-dim(y)[3]]))#we'll update this one as we go
    for (j in 1:ncol(y))
    {
      data[,colnames(y)[j]] = c(y[,j,-dim(y)[3]]*dt)
    }
    for (j in 1:ncol(x))
    {
        data[,colnames(x)[j]] = c(x[,j,1:(dim(y)[3]-1)]*dt)
    }
    
    res = as.data.frame(matrix(NA,nrow=nrow(data),ncol=ncol(y)))
    colnames(res)=colnames(y)
    for (j in 1:ncol(y))
    {
      #note: all predictors must be scaled by dt!
      #also need weights
      
      #update predictor, the rest stay the same
      data[,"dyn"] = c(y[,j,-1]-y[,j,-dim(y)[3]])
      
      m = FitFun(dyn~.-1,data=data,weights=weights)
      
      #drop any asshole "`" which for some reason lm wants to put on
      C = coef(m)
      for (jj in 1:length(C)) names(C)[jj] = gsub("`","",names(C)[jj])
      

      W[j,] = C[colnames(y)]
      lambda[j,] = -C[colnames(lambda)]
      pr = predict(m,newdata=data)
      #sigma[j,j] = mean((data$dyn-pr)^2,na.rm=T) #just in case... I've heard things about residuals() or lm$res
      res[,colnames(y)[j]] = data$dyn-pr
      
      if(options[["error"]])
      {
        print("errors: to do")
        warning("error: to do")
        #not sure how to do dlambda...
        dW[j,] = NA
        dlambda[j,] = NA
        dqterm = NA
      }
      else
      {
        dW = NULL
        dlambda = NULL
      }
    }
    lambda = solve(W)%*%lambda
  }

  #now estimate noise
  if(calcNoise)
  {
    sf = list(W=W,lambda=lambda)
    class(sf)="SF"
    
    if(options[["noiseMethod"]]=="y0") #use unimputed values
    {
      y0 = options[["y0"]]
      if(is.null(y0)) stop("You must provide y0 to options if you want to use noiseMethod='y0'")
      Q = FitQ.SF(y=y0,x=x,dt=dt,sf=sf,options=options[["qOptions"]])
      sigma = Q[["sigma"]]
      Q = Q[["Q"]]
    }
    else if (any(grepl("scale",tolower(options[["noiseMethod"]])))) #scale by # of imputed values (can be biased)
    {
      if(is.null(y0)) stop("You must provide y0 to options if you want to use noiseMethod='scale'")
      Q = FitQ.SF(y=y,x=x,dt=dt,sf=sf,options=options[["qOptions"]])
      sigma = Q[["sigma"]]*sum(!is.na(y))/(sum(!is.na(y0))-1)
      Q = Q[["Q"]]/sum(!is.na(y))*(sum(!is.na(y0))-1)
    }
    else if (any(grepl("y",tolower(options[["noiseMethod"]])))) #use imputed values (can be biased)
    {
      Q = FitQ.SF(y=y,x=x,dt=dt,sf=sf,options=options[["qOptions"]])    
      sigma = Q[["sigma"]]
      Q = Q[["Q"]]
    }
  }
  else
  {
    sigma = solve(Q)
  }
  
  l = list(W=W,lambda=lambda,sigma=sigma,Q=Q,
           dW=dW,dlambda=dlambda)
  l[["options"]] = options
  return(l)
}

AymErrorMoments.SF = function(y,x,dt,
                              W=NULL,lambda=NULL,Q=diag(1,ncol(y),ncol(y)), #initial parameter values
                              calcNoise=T,
                              options=list()
)
{
  #default options
  defaults = list(
    noiseMethod="y", #options: y, y0, scaled #you must provide y0 if you pick y0 (add to options)
    l2W=0,
    initLambdaZero=F, #specifically initialize lambda as 0s (i.e. exponential model)
    lambdaOptions = list(),
    wOptions = list(),
    qOptions = list()
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  options[["lambdaOptions"]][["error"]] = "asymptotic"
  options[["wOptions"]][["error"]] = "asymptotic"
  
  return(FitMoments.SF(y=y, x=x, dt=dt, W=W, lambda=lambda, Q=Q, calcNoise=calcNoise, options=options))
}


FitMoments.SF = function(y,x,dt,
                        W=NULL,lambda=NULL,Q=diag(1,ncol(y),ncol(y)), #initial parameter values
                        calcNoise=T,
                        options=list())
{
  #default options
  defaults = list(
                  noiseMethod="y", #options: y, y0, scaled #you must provide y0 if you pick y0 (add to options)
                  l2W=0,
                  initLambdaZero=F, #specifically initialize lambda as 0s (i.e. exponential model)
                  lambdaOptions = list(),
                  wOptions = list(),
                  qOptions = list(),
                  skipFit=F
                  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]

  
  
  #1. fit lambda first
  if(options[["initLambdaZero"]])
  {
    warning("starting lambda at 0 ...")
    lambda = matrix(0,nrow=ncol(y),ncol=ncol(x))
    colnames(lambda) = colnames(x)
    rownames(lambda) = colnames(y)
  }
  else if(is.null(lambda))
  {
    lambda = FitLambdaEq.SF(y=y,x=x,dt=dt,Q=Q,options=options[["lambdaOptions"]])
    dlambda = lambda[["dlambda"]]
    lambda = lambda[["lambda"]]
  }
  else if (is.null(W))
  {
    warning("Lambda provided but not W, skipping lambda estimation...")
  }
  else
  {
    lambda = FitLambda.SF(y=y,x=x,dt=dt,W=W,Q=Q,options=options[["lambdaOptions"]])
    dlambda = lambda[["dlambda"]]
    lambda = lambda[["lambda"]]
  }
  
  #2. now fit W
  W = FitW.SF(y=y,x=x,dt=dt,lambda=lambda,Q=Q,options=options[["wOptions"]])
  dW = W[["dW"]]
  W = W[["W"]]

  
  #3. now estimate noise
  if(calcNoise)
  {
    sf = list(W=W,lambda=lambda)
    class(sf)="SF"
    
    if(options[["noiseMethod"]]=="y0") #use unimputed values
    {
      y0 = options[["y0"]]
      if(is.null(y0)) stop("You must provide y0 to options if you want to use noiseMethod='y0'")
      Q = FitQ.SF(y=y0,x=x,dt=dt,sf=sf,options=options[["qOptions"]])
      sigma = Q[["sigma"]]
      Q = Q[["Q"]]
    }
    else if (any(grepl("scale",tolower(options[["noiseMethod"]])))) #scale by # of imputed values (can be biased)
    {
      if(is.null(y0)) stop("You must provide y0 to options if you want to use noiseMethod='scale'")
      Q = FitQ.SF(y=y,x=x,dt=dt,sf=sf,options=options[["qOptions"]])
      sigma = Q[["sigma"]]*sum(!is.na(y))/(sum(!is.na(y0))-1)
      Q = Q[["Q"]]/sum(!is.na(y))*(sum(!is.na(y0))-1)
    }
    else if (any(grepl("y",tolower(options[["noiseMethod"]])))) #use imputed values (can be biased)
    {
      Q = FitQ.SF(y=y,x=x,dt=dt,sf=sf,options=options[["qOptions"]])    
      sigma = Q[["sigma"]]
      Q = Q[["Q"]]
    }
  }
  else
  {
    sigma = solve(Q)
  }

  #fits SF model using moments
  l = list(W=W,lambda=lambda,Q=Q,sigma=sigma,
           dW=dW,dlambda=dlambda,dQ=NULL,dsigma=NULL,options=options)
  class(l)="SF"
  return(l)
}

FitLambdaEq.SF = function(y, #individual x variable x time
                          x, #covariates
                          dt,
                          Q=diag(1,ncol(y)),
                          options=list()
)
{
  #equilibrium solution
    #does not require W
  #i.e. when dy = 0
  
  defaults = list(l2=0, #l2 penalty #warning: unvalidated
                  use = "pairwise.complete" #missingness handling; sends to cov()
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  time = dim(y)[3]-1
  
  ynt = list()
  xnt = list()
  xn = list()
  dysign = list()
  for (i in 1:time)
  {
    xmat = matrix(x[,,i],nrow=nrow(x),ncol=ncol(x))#in case x has only 1 column
    
    xt = xmat
    for (j in 1:ncol(xt)) xt[,j] = xt[,j]*sqrt(abs(dt[,i]))
    xnt[[i]] = xt
    rm(xt)
    
    xn[[i]] = xmat
    
    yt = matrix(y[,,i],nrow=nrow(y),ncol=ncol(y)) #in case y has only 1 column
    for (j in 1:ncol(yt)) yt[,j] = yt[,j]*sqrt(abs(dt[,i]))
    ynt[[i]] = yt
    rm(yt)
    
  }
  xnt = do.call(rbind,xnt)
  xn = do.call(rbind,xn)
  ynt = do.call(rbind,ynt)
  
  Cx = cov2(xnt,use=options[["use"]],center=F)
  Cxy = cov2(xnt,ynt,use=options[["use"]],center=F)
  
  l2 = options[["l2"]]
  l2 = l2 / nrow(y) / (dim(y)[3] - 1)
  
  lambda = t(solve(Cx+diag(l2,ncol(Cx),ncol(Cx)),Cxy)) #l2 assumes conjugate hyper-covariance
  
  colnames(lambda) = colnames(x)
  rownames(lambda) = colnames(y)
  
  l = list(lambda=lambda,dlambda=NULL,options=options)
  
  return(l)
}

FitLambda.SF = function(y, #individual x variable x time
                        x, #covariates, dim(x)[3] = dim(y)[3]-1 #inbetween timepoints i.e. should be evaluated at t+dt/2 
                        dt,
                        W,
                        Q = diag(1,ncol(y)),
                        options = list()
)
{
  
  defaults = list(l2=0, #l2 penalty #warning: unvalidated
                  prior_mu=matrix(0,nrow=ncol(y),ncol=ncol(x)), #dimensions are ncol(y) x ncol(x) (same as lambda) #unvalidated
                  use = "pairwise.complete", #missingness handling; sends to cov()
                  diagonal=F,  #use if W and noise are diagonal
                  error ="none",
                  nboot=100
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  l20 = options[["l2"]]
  WTQ = t(W)%*%Q
  WTQW = WTQ%*%W
  
  time = dim(y)[3]-1
  
  ynt = list()
  xnt = list() #x*sqrt(abs(dt))
  xn = list()
  dysign = list()
  for (i in 1:time)
  {
    xmat = matrix(x[,,i],nrow=nrow(x),ncol=ncol(x))#as.matrix(x[,,i]) #in case x has only 1 column
    
    xt = xmat
    for (j in 1:ncol(xt)) xt[,j] = xt[,j]*sqrt(abs(dt[,i]))
    xnt[[i]] = xt
    rm(xt)
    
    xn[[i]] = xmat
    
    yt = matrix(y[,,i],nrow=nrow(y),ncol=ncol(y)) #as.matrix(y[,,i]) #in case y has only 1 column
    for (j in 1:ncol(yt)) yt[,j] = yt[,j]*sqrt(abs(dt[,i]))
    ynt[[i]] = yt
    rm(yt)
    
    dyst = as.matrix(y[,,i+1])-as.matrix(y[,,i])
    for (j in 1:ncol(dyst)) dyst[,j] = dyst[,j]*sign(dt[,i])
    dysign[[i]] = dyst
    rm(dyst)
    
  }
  xnt = do.call(rbind,xnt) #weighted by sqrt(abs(dt))
  xn = do.call(rbind,xn)
  ynt = do.call(rbind,ynt) #weighted by sqrt(abs(dt))
  dysign = do.call(rbind,dysign) #weighted by sign(dt)
  
  Cx = cov2(xnt,use=defaults[["use"]],center=F)
  Cdyx = cov2(dysign,xn,use=defaults[["use"]],center=F)
  Cyx = cov2(ynt,xnt,use=defaults[["use"]],center=F)
  
  l2 = options[["l2"]]
  prior_mu = options[["prior_mu"]]
  if(options[["diagonal"]]) #if W and Q are both diagonal
  {
    l2 = l2 / nrow(y) / (dim(y)[3] - 1)
    warning("l2 and weights are not validated for diagonal case")
    
    #weight by W eigenvalues
    wCyx = Cyx
    for (i in 1:nrow(wCyx)) wCyx[i,] = Cyx[i,]*diag(W)
    
    #fix weighting
    for (i in 1:ncol(lambda)) lambda[,i] = lambda[,i]*diag(W)
  }
  else
  {
    l2 = diag(l2,ncol(y),ncol(y))
    l2 = l2 / nrow(y) / (dim(y)[3] - 1)
    Qx = solve(Cx) 
    
    lambda = solve(WTQW+l2,(WTQW%*%Cyx-WTQ%*%Cdyx)%*%Qx+l2%*%prior_mu)
  }
  
  colnames(lambda) = colnames(x)
  rownames(lambda) = colnames(y)
 
  #estimate errors
  dlambda = NULL 
  if(grepl("asy",options[["error"]]) | grepl("anal",options[["error"]]))
  {
    WTQWinv = solve(WTQW)
    Cxinv = solve(Cx)
    dlambda = sqrt( outer(diag(WTQWinv),diag(Cxinv))/nrow(y)/time )
    
    dimnames(dlambda) = dimnames(lambda)
  }
  else if(grepl("boot",options[["error"]]))
  {
    names = dimnames(lambda)
    lambda = list(lambda)
    for (i in 2:nboot)
    {
      ind = sample(1:nrow(y),replace=T)
      subOptions = options
      subOptions[["error"]] = "none"
      lambda[[i]] = FitSF.lambda(y=y, x=x, dt=dt, W=W, Q=Q, options=subOptions)
    }
    lambda = ListMeanSD(lambda)
    dimnames(lambda[[1]]) = names
    dimnames(lambda[[2]]) = names
  }
  l = list(lambda=lambda,dlambda=dlambda,options=options)
  return(l)
}

FitW.SF = function(y, #individual x variable x time
                   x, #covariates, dim(x)[3] = dim(y)[3]-1 #inbetween timepoints i.e. should be evaluated at t+dt/2 
                   dt,
                   lambda,
                   Q, #only needed for analytical error estimates
                   Qinv=NULL, #a little faster if pre-computed
                   options=list()
                  )
{

  #set defaults
  defaults = list(l2=0, #l2 penalty for W #can be value or matrix of same size as W #warning: unvalidated
                  prior_mu=matrix(0,ncol(y),ncol(y)), #number of matrix the size of W (ncol(y) x ncol(y)) #unvalidated
                  use = "pairwise.complete", #missingness handling; sends to cov()
                  error="none", #which error to compute? bootstrap, asymptotic or none?
                  nboot=100, #only used if error=bootstrap
                  diagonal=F #force W to be diagonal?
              )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  l2 = options[["l2"]]
  if(is.null(dim(l2))) 
  {
    l2 = diag(l2,ncol(y),ncol(y))
  }
  prior_mu = options[["prior_mu"]]
  if(is.null(dim(prior_mu))) 
  {
    if(options[["diagonal"]])
    {
      prior_mu = diag(prior_mu,ncol(y),ncol(y))  
    }
    else prior_mu = matrix(prior_mu,ncol(y),ncol(y))
  }
  
  priormu_l2 = prior_mu%*%l2 
  l2_priormu = t(priormu_l2)
  
  time = dim(y)[3]-1
  
  ynt = list() #(yn-mu)*sqrt(abs(dt))
  yntsign = list() #(yn-mu)*sign(dt)
  dy = list()
  for (i in 1:time)
  {
    mun = x[,,i]%*%t(lambda)
    
    yst = as.matrix(y[,,i]) #in case x has only 1 column
    for (j in 1:ncol(yst)) yst[,j] = (yst[,j]-mun[,j])*sign(dt[,i])
    yntsign[[i]] = yst
    rm(yst)
    
    yt = as.matrix(y[,,i])
    for (j in 1:ncol(yt)) yt[,j] = (yt[,j]-mun[,j])*sqrt(abs(dt[,i]))
    ynt[[i]] = yt
    rm(yt)
    
    dy[[i]] = as.matrix(y[,,i+1])-as.matrix(y[,,i])
  }
  ynt = do.call(rbind,ynt)
  yntsign = do.call(rbind,yntsign)
  dy = do.call(rbind,dy)
  
  if(options[["diagonal"]])
  {
    l2 = l2 / nrow(y) / time
    priormu_l2 = priormu_l2/ nrow(y) / time
    #guessing where to put l2
    Cy = apply(ynt^2,2,mean,na.rm=T)
    l = (apply(dy*yntsign,2,mean,na.rm=T)+diag(priormu_l2))/(Cy + diag(l2))
    
    W = diag(l,ncol(y),ncol(y))
  }
  else
  {
    Cy = cov2(ynt,use=defaults[["use"]],center=F)
    Cy0y1 = cov2(yntsign,dy,use=defaults[["use"]],center=F)
    
    
    l2 = l2 / nrow(y) / time
    l2_priormu = l2_priormu/ nrow(y) / time

    W = t(solve(Cy+l2,Cy0y1+l2_priormu))

  }
  
  dW = NULL
  #compute errors?
  if(grepl("asy",options[["error"]]) | grepl("anal",options[["error"]]))
  {
    if(is.null(Qinv)) sig = solve(Q)
    else sig = Qinv
    
    if(options[["diagonal"]]) Cyinv = 1/Cy
    else Cyinv = solve(Cy)
    if(options[["diagonal"]]) dW = sqrt(outer(diag(sig),Cyinv)/time/nrow(y))
    else dW = sqrt(outer(diag(sig),diag(Cyinv))/time/nrow(y))
    
    dimnames(dW) = dimnames(W)
  }
  else if(grepl("boot",options[["error"]]))
  {
    names = dimnames(W)
    W = list(W)
    for (i in 2:nboot)
    {
      ind = sample(1:nrow(y),replace=T)
      subOptions = options
      subOptions[["error"]] = "none"
      W[[i]] = FitW.SF(y=y, x=x, dt=dt, lambda=lambda, Q=Q, Qinv=Qinv, options=subOptions)
    }
    W = ListMeanSD(W,sem=F,na.rm=T)
    dimnames(W[[1]]) = names
    dimnames(W[[2]]) = names
  }
  
  l = list(W=W,dW=dW,options=options)
  return(l)
}

FitQ.SF = function(y,
                       x,
                       dt,
                       sf, #sf fit, as bare minimum list including: W, lambda, trans and transinv
                       options=list(),
                       predictOptions=list() #sent to predict.SF
)
{
  #inverse covariance of the residual (noise)
  #i.e. Q
  
  #set defaults
  defaults = list(
                  diagonal=F #force W to be diagonal?
  )
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]

  pr = predict.SF(sf=sf,y=y,x=x,dt=dt,options=predictOptions)
  residual = y[,,-1,drop=F]-pr
  
  #divide out dt so it cancels
  for (j in 1:ncol(residual)) for (k in 1:dim(residual)[3]) residual[,j,k] = residual[,j,k]/sqrt(abs(dt[,k]))
  
  if(options[["diagonal"]])
  {
    C = matrix(0,nrow=ncol(y),ncol=ncol(y)) #covariance
    rownames(C) = colnames(y)
    colnames(C) = rownames(C)
    for (i in 1:ncol(y))
    {
      C[i,i] = mean(residual[,i,]^2,na.rm=T) #probably a biased estimator
    }
    
    Q = diag(1/diag(C),nrow(C),ncol(C))
  }
  else
  {
    C = matrix(NA,nrow=ncol(y),ncol=ncol(y)) #covariance
    rownames(C) = colnames(y)
    colnames(C) = rownames(C)
    for (i in 1:ncol(y))
    {
      for (j in 1:ncol(y))
      {
        C[i,j] = mean(residual[,i,]*residual[,j,],na.rm=T) #probably a biased estimator
      }
    }
    Q = solve(C)
  }
  
  l = list(Q=Q,sigma=C,options=options)
  
  return(l)
}

predict.SF = function(sf,
                      y, #data
                      x, #covariates
                      dt, #time steps
                      options=list() #none (yet)
)
{
  #predicts next timepoint using data from current timepoint
  
  #there are no defaults
  #defaults = list()
  #for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  #if(class(sf)!="SF") stop("first argument should be SF class. if this is in error set class(sf) = 'SF' ") #never implemented
  
  W = sf$W
  lambda = sf$lambda
  
  #y is 3d array of multiple times:
  if(length(dim(y)) > 2)
  {
    timeInds = 1:(dim(y)[3]-1)
    
    pr = y[,,timeInds,drop=F]
    for (k in 1:length(timeInds))
    {
      xmat = matrix(x[,,timeInds[k]],nrow=nrow(x),ncol=ncol(x)) #as.matrix(x[,,timeInds[k]])
      mulambda = xmat%*%t(lambda)
      
      y0 = matrix(y[,,timeInds[k]],nrow=nrow(y),ncol=ncol(y)) #in case y has only 1 column #as.matrix(y[,,timeInds[k]]) 
      for (i in 1:ncol(y0)) y0[,i] = y0[,i]-mulambda[,i]
      pr[,,k] = y0%*%t(W)
      for (i in 1:ncol(pr)) pr[,i,k] = y[,i,timeInds[k]] + dt[,timeInds[k]]*pr[,i,k]
      
    }
    
  }
  else #y is flat
  {
    timeInds = 1
    if (!is.null(dim(dt))) dt = dt[,timeInds]
    
    xmat = matrix(x[,,timeInds],nrow=nrow(x),ncol=ncol(x)) #as.matrix(x[,,timeInds])
    mulambda = xmat%*%t(lambda)
    
    y0 = y
    for (i in 1:ncol(y0)) y0[,i] = y0[,i]-mulambda[,i]
    pr = y0%*%t(W)
    for (i in 1:ncol(pr)) pr[,i] = y[,i] + dt*pr[,i]
    
  }
  
  return(pr)
}

FitOP.SF = function(y,
                    x,
                    dt,
                    W = NULL,#initial guess
                    lambda = NULL,#initial guess
                    Q = diag(1,ncol(y)), #initial guess
                    calcNoise=F, #do you want to fix the noise or not?
                    options=list()
)
{
  warning("FitOP.SF is unvalidated")
  #default options
  defaults = list(
    error=F,
    method="L-BFGS-B",
    Q = diag(1,ncol(y),ncol(y)), #starting guess and/or fixed value depending on fixedQ
    initializeWith="lm", #"lm" or "0" or defaults to moments
    minVal = -1e20, #prevents -Inf
    lower=NULL, #lower limit for parameters
    upper=NULL, #upper limit for parameters
    l1W=NULL,
    addIntercept=F,
    frozenW = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient)
    frozenLambda = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient)
    frozenQ = NULL, # matrix(T,ncol(y),ncol(y)), #matrix of values to freeze (bool, T: freeze at 0 gradient)
    relax=T, #fit a second time after applying l1W penalty
    zeroCutoff=1e-3
  )
    
  for (i in 1:length(defaults)) if(is.null(options[[names(defaults)[i]]])) options[[names(defaults)[i]]] = defaults[[i]]
  
  error=options[["error"]]
  method=options[["method"]]
  Q = options[["Q"]]
  lambda0=options[["lambda0"]]
  W0=options[["W0"]]
  initializeWith=options[["initializeWith"]]
  minVal =options[["minVal"]]
  lower=options[["lower"]]
  upper=options[["upper"]]
  l1W=options[["l1W"]]
  addIntercept=options[["addIntercept"]]
  fixedQ=calcNoise
  frozenW = options[["frozenW"]]
  frozenLambda = options[["frozenLambda"]]
  frozenQ = options[["frozenQ"]]
  relax=options[["relax"]]
  zeroCutoff=options[["zeroCutoff"]]
  
    if(relax & is.null(l1W))
    {
      warning("relax will only be applid if l1W is supplied")
      relax=F
    }
    
    #add intercept to x if not already present
    if(!("mu0"%in%colnames(x)) & addIntercept)
    {
      xpp = array(1,dim=dim(x)+c(0,1,0))
      xpp[,-1,] = x
      colnames(xpp)=c("mu0",colnames(x)) 
    }
    else xpp = x
    
    #starting guesses
    if(tolower(initializeWith)=="lm")
    {
      m = FitLM.SF(y=y,x=xpp,dt=dt,calcNoise=F,options=list(diagonalW=F))
      lambda0 = m[["lambda"]]
      W = m[["W"]]
      if(is.null(W))
      {
        W = m[["W"]]
      }
      if(is.null(lambda))
      {
        lambda = m[["lambda"]]
      }
    }
    else if(tolower(initializeWith)=="zeros" | tolower(initializeWith) == "0")
    {
      if(is.null(W))
      {
        W = diag(-.25,ncol(y),ncol(y))
        colnames(W) = colnames(y)
        rownames(W) = colnames(y)
      }
      if(is.null(lambda))
      {
        lambda = matrix(0,nrow=ncol(y),ncol=ncol(xpp))
        colnames(lambda) = colnames(xpp)
        rownames(lambda) = colnames(y)    
      }
    }
    else
    {
      m = FitMoments.SF(y=y,x=xpp,dt=dt,calcNoise=F)
      if(is.null(W))
      {
        W = m[["W"]]
      }
      if(is.null(lambda))
      {
        lambda = m[["lambda"]]
      }
    }
    

    if(fixedQ)
    {
      frozenQ = matrix(T,nrow=ncol(y),ncol=ncol(y))
      rownames(frozenQ) = colnames(y)
      colnames(frozenQ) = colnames(y)
    }
  
    W0 = W
    lambda0=lambda
    Q0 = Q
    
    #number of variables
    Nw = length(W0)
    Nlambda = length(lambda0)
    Nq = length(Q0)
    
    lwrapper = function(par)
    {
      Wpar = matrix(par[1:Nw],ncol=ncol(y),nrow=ncol(y))
      lambdapar = matrix(par[1:Nlambda+Nw],nrow=ncol(y),ncol=ncol(xpp))
      Qpar = matrix(par[1:Nq+Nlambda+Nw],nrow=ncol(y),ncol=ncol(y))
      ll = loglik.SF(y=y,x=xpp,dt=dt,lambda=lambdapar,W=Wpar,Q=Qpar,skipZ=F,l1W=l1W,
                     frozenW=frozenW,frozenLambda=frozenLambda,frozenQ=frozenQ)
      #print(ll)
      if(is.na(ll) | is.nan(ll))
      {
        #print("ll:")
        #print(ll)
        #print(par)
        ll = minVal
      }
      if(ll < minVal) ll = minVal
      #print(ll)
      return(-ll)
    }
    gradwrapper = function(par)
    {
      Wpar = matrix(par[1:Nw],ncol=ncol(y),nrow=ncol(y))
      lambdapar = matrix(par[1:Nlambda+Nw],nrow=ncol(y),ncol=ncol(xpp))
      Qpar = matrix(par[1:Nq+Nlambda+Nw],nrow=ncol(y),ncol=ncol(y))
      gr = loglik.grad.SF(y=y,x=xpp,dt=dt,lambda=lambdapar,W=Wpar,Q=Qpar,skipZ=F,l1W=l1W,
                          frozenW=frozenW,frozenLambda=frozenLambda,frozenQ=frozenQ)
      
      #drop Q (for now)
      gr = gr[1:(Nlambda+Nw+Nq)]
      gr = -gr
      
      #drop NAs -- how???
      #gr[is.na(gr)] = 
      
      #drop -Inf
      logi = abs(gr) > abs(minVal)
      if(any(logi)) gr[logi] = abs(minVal)*sign(gr[logi])
      
      return(gr)
    }
    
    #first entry is intercept (W0, mu0, Q) 
    par0 = c(c(W0),c(lambda0),c(Q0))
    
    
    if (grepl("stogo",tolower(method))) #having trouble with Q
    {
      library(nloptr)
      if(is.null(lower))
      {
        q = quantile(abs(c(dt)),probs=.25,na.rm=T)
        #lower_W = -rep(2/q,Nw)
        #lower_lambda = rep(-10,Nlambda)
        #lower_Q = rep(1e-3,Nq)
        lower_W = -10*rep(max(abs(c(W0))),Nw)
        lower_lambda = -10*rep(max(abs(c(lambda0))),Nlambda)
        lower_Q = diag(min(abs(diag(Q0)))/10,ncol(Q0)) #off diagonals can be 0 but diagonals can't
        lower_Q = c(lower_Q)
        lower = c(lower_W,lower_lambda,lower_Q)
      }
      if(is.null(upper))
      {
        q = quantile(abs(c(dt)),probs=.25,na.rm=T)
        #upper_W = rep(2/q,Nw)
        #upper_lambda = rep(10,Nlambda)
        #upper_Q = rep(1e3,Nq)
        upper_W = 10*rep(max(abs(c(W0))),Nw)
        upper_lambda = 10*rep(max(abs(c(lambda0))),Nlambda)
        upper_Q = rep(max(abs(c(Q0))),Nq)*10 #off diagonals can't actually be as big but don't worry about it
        upper = c(upper_W,upper_lambda,upper_Q)
      }
      
      #print(lower)
      #print(par0)
      #print(par0 < lower)
      #print(par0 > upper)
      
      opts = opts = list("algorithm"="NLOPT_GD_STOGO") #"xtol_rel"=1.0e-8
      #opts = opts = list("algorithm"="NLOPT_GD_STOGO_RAND") #not sure how this is different...
      
      op = nloptr(x0=par0,eval_f=lwrapper,eval_grad_f=gradwrapper,lb=lower,ub=upper,opts=opts)  
      
      if(op$status != 1) warning(sprintf("FAILED TO CONVERGE (%s)",op$message))
      par = op$solution
      status = op$status
      fisher = NULL
      if(error) warning("Errors not supported for nloptr")
    }
    else if (grepl("mlsl",tolower(method))) #seems ok, so long as you initialize with lm
    {
      library(nloptr)
      if(is.null(lower))
      {
        q = quantile(abs(c(dt)),probs=.25,na.rm=T)
        #lower_W = -rep(2/q,Nw)
        #lower_lambda = rep(-10,Nlambda)
        #lower_Q = rep(1e-3,Nq)
        lower_W = -10*rep(max(abs(c(W0))),Nw)
        lower_lambda = -10*rep(max(abs(c(lambda0))),Nlambda)
        lower_Q = diag(min(abs(diag(Q0)))/10,ncol(Q0)) #off diagonals can be 0 but diagonals can't
        lower_Q = c(lower_Q)
        lower = c(lower_W,lower_lambda,lower_Q)
      }
      if(is.null(upper))
      {
        q = quantile(abs(c(dt)),probs=.25,na.rm=T)
        #upper_W = rep(2/q,Nw)
        #upper_lambda = rep(10,Nlambda)
        #upper_Q = rep(1e3,Nq)
        upper_W = 10*rep(max(abs(c(W0))),Nw)
        upper_lambda = 10*rep(max(abs(c(lambda0))),Nlambda)
        upper_Q = rep(max(abs(c(Q0))),Nq)*10 #off diagonals can't actually be as big but don't worry about it
        upper = c(upper_W,upper_lambda,upper_Q)
      }
      
      
      op = nloptr(x0=par0,eval_f=lwrapper,eval_grad_f=gradwrapper,lb=lower,ub=upper,opts=list("algorithm"="NLOPT_GD_MLSL",local_opts=list("algorithm"="NLOPT_LD_LBFGS")))
      
      if(op$status != 1) warning(sprintf("FAILED TO CONVERGE (%s)",op$message))
      par = op$solution
      status = op$status
      fisher = NULL
      if(error) warning("Errors not supported for nloptr")
    }
    else #usually sufficient
    {
      op = optim(par=par0,fn=lwrapper,gr=gradwrapper,method=method,hessian=error)
      par = op$par
      if(op$convergence != 0) warning(sprintf("FAILED TO CONVERGE (%s)",op$message))
      status = op$convergence
      #hessian is negative fisher information - but since I evaluate -loglikelihood it IS the fisher inforamtion
      fisher = op$hessian
    }
    #status: 0 is success for optim and nloptr
    
    
    W = matrix(par[1:Nw],ncol=ncol(y),nrow=ncol(y))
    colnames(W) = colnames(y)
    rownames(W) = colnames(y)
    lambda = matrix(par[1:Nlambda+Nw],nrow=ncol(y),ncol=ncol(xpp))
    colnames(lambda) = colnames(xpp)
    rownames(lambda) = colnames(y)
    
    Q = matrix(par[1:Nq+Nlambda+Nw],nrow=ncol(y),ncol=ncol(y))
    colnames(Q) = colnames(y)
    rownames(Q) = colnames(y)
    
    if(relax)
    {
      #call function a second time, but constrain W
      zeros = abs(W) < zeroCutoff #freeze these at 0
      W[zeros] = 0
      if(is.null(frozenW))
      {
        frozenW = zeros
      }
      else
      {
        frozenW = zeros | frozenW
      }
      res = FitOP.sf = function(y=y,
                                x=x,
                                dt=dt,
                                error=error,
                                method=method,
                                Q = Q, 
                                lambda0=lambda,
                                W0=W, 
                                initializeWith="default", #should already have a good starting point
                                minVal = minVal, 
                                lower=lower,
                                upper=upper, 
                                l1W=NULL,  #the whole point is to fit without this!
                                addIntercept=F, #will already be done
                                fixedQ=fixedQ, 
                                frozenW = frozenW,
                                frozenLambda = frozenLambda, 
                                frozenQ = frozenQ, 
                                relax=F, #don't get stuck in infinite loop
                                zeroCutoff=zeroCutoff)
      return(res)
    }
    
    #hessian
    dW = NA*W
    dlambda = NA*lambda   
    dQ = NA*Q 
    errorMat=NA
    if(error & !is.null(fisher))
    {
      warning("Errors are unvalidated")
      errorMat = tryCatch(solve(fisher), error=function(e) return(NA))
      
      if(all(is.na(errorMat)))
      {
        dW = NA*W
        dlambda = NA*lambda   
        dQ = NA*Q 
      }
      else
      {
        dW = matrix(sqrt(diag(errorMat)[1:Nw]),nrow=nrow(W),ncol=ncol(W))
        rownames(dW) = rownames(W)
        colnames(dW) = colnames(W)
        
        dlambda = matrix(sqrt(diag(errorMat)[1:Nlambda+Nw]),nrow=nrow(lambda),ncol=ncol(lambda))
        rownames(dlambda) = rownames(lambda)
        colnames(dlambda) = colnames(lambda)
        
        dQ = matrix(sqrt(diag(errorMat)[1:Nq+Nlambda+Nw]),nrow=nrow(Q),ncol=ncol(Q))
        rownames(dQ) = rownames(Q)
        colnames(dQ) = colnames(Q)
      }
    }
    
    sigma=tryCatch(solve(Q), error=function(e) return(NA))
    if(all(is.na(sigma))) sigma=tryCatch(solve(Q+diag(1e-6,ncol(Q))), error=function(e) return(NA)) #try again with a little boost to prevent rounding errors
    
    l = list(W=W,lambda=lambda,Q=Q,sigma=sigma,
         dW=dW,dlambda=dlambda,dQ=dQ,
         lower=lower,upper=upper,
         Q0=Q0,W0=W0,lambda0=lambda0,
         errorMat=errorMat,
         fisher=fisher,op=op,status=status,
         frozenQ=frozenQ,frozenW=frozenW,frozenLambda=frozenLambda)
    
    return(l)
  }


loglik.grad.SF = function(y,
                          x,
                          W,
                          lambda,
                          Q=diag(1,ncol(y),ncol(y)),
                          sigma=solve(Q),
                          dt=matrix(1,nrow=nrow(y),ncol=dim(y)[3]-1),
                          logdetQ=NULL,
                          skipZ=F,
                          pr0 = NULL, #optionally provide prediction (for speed)
                          ave=T, #scale by number of datum?
                          l1W = NULL, #l1 penalty (matrix or number)
                          frozenW = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient)
                          frozenLambda = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient)
                          frozenQ = NULL, # matrix(T,ncol(y),ncol(y)), #matrix of values to freeze (bool, T: freeze at 0 gradient)
                          biasCorrect=F #correct bias in Q https://stats.stackexchange.com/questions/167494/bias-in-the-mle-of-variance-component-in-a-multivariate-gaussian
)
{
  #validation: looks right (at least in 1d)
  if(is.null(logdetQ) & !skipZ) #check if we needs this for Q
  {
    logdetQ = tryCatch(determinant(Q,T),error=function(e){return(list(modulus=NA,sign=NA))})
    
    logdetQ = as.numeric(logdetQ$sign*logdetQ$modulus)
  }
  timeInds = 1:(dim(y)[3]-1)
  l = 0
  grW = rep(0,length(W))
  grLambda = rep(0,length(lambda))
  grQ = rep(0,length(Q))
  Nobs = 0 #update Aug 2023 #will be same as before so long as no NAs, otherwise smaller
  for (k in 1:length(timeInds))
  {
    xmat = matrix(x[,,timeInds[k]],nrow=nrow(x),ncol=ncol(x))
    mulambda = xmat%*%t(lambda)
    
    if(is.null(pr0)) #update Nov 24, 2022
    {
      #prediction of yn+1
      y0 = as.matrix(y[,,timeInds[k]]) #in case y has only 1 column
      for (j in 1:ncol(y0)) y0[,j] = y0[,j]-mulambda[,j]
      pr = y0%*%t(W)
      for (j in 1:ncol(pr)) pr[,j] = y[,j,timeInds[k]] + dt[,timeInds[k]]*pr[,j]     
    }
    else pr = pr0[,,timeInds[k]]
    
    dy = (y[,,timeInds[k]+1]-pr)
    Nobs = Nobs + sum(!is.na(dy))
    #for W
    dyQ = dy%*%Q
    #for lambda
    dyQW = dy%*%Q%*%W
    #probably faster to do for loop vs outer then sum
    #drop NAs
    dyQ[is.na(dyQ)] = 0
    dyQW[is.na(dyQW)] = 0
    y0[is.na(y0)] = 0
    dy[is.na(dy)] = 0
    for (i in 1:nrow(xmat)) 
    {
      #update W
      grW = grW + sign(dt[i,timeInds[k]])*c(outer(dyQ[i,],y0[i,]))
      #grW = grW + sign(dt[i,timeInds[k]])*c(outer(dyQ[i,],y[,,timeInds[k]]-c(lambda%*%xmat[i,])))
      
      #update lambda
      grLambda = grLambda - sign(dt[i,timeInds[k]])*c(outer(dyQW[i,],xmat[i,]))
      
      #update Q 
      grQ = grQ - .5*c(outer(dy[i,],dy[i,]))/abs(dt[i,timeInds[k]])
    }
    
  }
  #Z contribution
  if(biasCorrect) grQ = grQ + (Nobs-length(grQ)-length(grLambda))/2*c(sigma) #ad hoc #makes sense to be though
  else grQ = grQ + Nobs/2*c(sigma)
  
  if(!is.null(frozenW))
  {
    grW[c(frozenW)] = 0
    l1penalty = sum(l1W*sign(W[!frozenW]))
  }
  else l1penalty = sum(l1W*sign(W))
  if(!is.null(frozenLambda))
  {
    grLambda[c(frozenLambda)] = 0
  }
  if(!is.null(frozenQ))
  {
    grQ[c(frozenQ)] = 0
  }
  if(ave) 
  {
    grW = grW/Nobs
    grLambda = grLambda/Nobs
    grQ = grQ/Nobs
  }
  if(!is.null(l1W)) grW = grW - l1penalty #always after ave
  gr = c(grW,grLambda,grQ)
  
  names(gr) = c(sprintf("W%02d",1:length(grW)),c(outer(colnames(y),colnames(x),paste,sep="_")),sprintf("Q%02d",1:length(grQ)))
  
  return(gr)
}


loglik.SF = function(y,
                     x,
                     W,
                     lambda,
                     mu0=rep(0,ncol(y)),
                     Q=diag(1,ncol(y),ncol(y)),
                     dt=matrix(1,nrow=nrow(y),ncol=dim(y)[3]-1),
                     logdetQ=NULL,
                     skipZ=F,
                     pr0 = NULL, #optionally provide prediction (for speed)
                     ave=T, #scale by number of datum?
                     l1W = NULL, #l1 penalty (matrix or number)
                     frozenW = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient) #only this one matters
                     frozenLambda = NULL, #matrix of values to freeze (bool, T: freeze at 0 gradient) #not used
                     frozenQ = NULL, #matrix(T,ncol(y),ncol(y)), #matrix of values to freeze (bool, T: freeze at 0 gradient) #not used
                     biasCorrect=F #correct bias in Q https://stats.stackexchange.com/questions/167494/bias-in-the-mle-of-variance-component-in-a-multivariate-gaussian
)
{
  if(is.null(logdetQ) & !skipZ)
  {
    logdetQ = tryCatch(determinant(Q,T),error=function(e){return(list(modulus=NA,sign=NA))})
    
    logdetQ = as.numeric(logdetQ$sign*logdetQ$modulus)
    if(is.na(logdetQ)) logdetQ = -1e100 #make likelihood very small
  }
  timeInds = 1:(dim(y)[3]-1)
  l = 0
  Nobs = 0 #update Aug 2023 #running sum
  if(biasCorrect) Nobs_est = sum(!is.na(y[,,-1]) & !is.na(y[,,-dim(y)[3]])) #need this for Z  #update Aug 2023
  for (k in 1:length(timeInds))
  {
    xmat = matrix(x[,,timeInds[k]],nrow=nrow(x),ncol=ncol(x))
    mulambda = xmat%*%t(lambda)
    
    if(is.null(pr0)) #update Nov 24, 2022
    {
      #prediction of yn+1
      y0 = as.matrix(y[,,timeInds[k]]) #in case y has only 1 column
      for (j in 1:ncol(y0)) y0[,j] = y0[,j]-mu0[j]-mulambda[,j]
      pr = y0%*%t(W)
      for (j in 1:ncol(pr)) pr[,j] = y[,j,timeInds[k]] + dt[,timeInds[k]]*pr[,j]     
    }
    else pr = pr0[,,timeInds[k]]
    
    dy = y[,,timeInds[k]+1]-pr
    obsLogi = !is.na(dy)
    Nobsn = sum(obsLogi)
    Nobs = Nobs + Nobsn
    
    if(skipZ) lnZ = 0
    else 
    {
      if(biasCorrect) #ad hoc implementation #https://stats.stackexchange.com/questions/167494/bias-in-the-mle-of-variance-component-in-a-multivariate-gaussian
      {
        lnZ = -.5*Nobsn*logdetQ*(Nobs_est-length(W)-length(lambda))/Nobs_est 
        for (j in 1:ncol(y)) lnZ = lnZ + 1/2*sum(log(2*pi*abs(dt[obsLogi[,j],timeInds[k]])),na.rm=T)
      }
      else
      {
        #lnZ = ncol(y)/2*sum(log(2*pi*abs(dt[,timeInds[k]])),na.rm=T)-.5*nrow(y)*logdetQ
        lnZ = -.5*Nobsn*logdetQ #aug 2023 update: accounts for dropped NAs
        for (j in 1:ncol(y)) lnZ = lnZ + 1/2*sum(log(2*pi*abs(dt[obsLogi[,j],timeInds[k]])),na.rm=T)
        #print("lnz:")
        #print(lnZ)
      }
    }
    #print("l")
    #print(- .5*sum(dy*(dy%*%Q),na.rm=T))
    
    #don't forget dt!    
    dysdt = dy
    for (j in 1:ncol(dysdt)) dysdt[,j] = dy[,j]/sqrt(abs(dt[,timeInds[k]]))
    #update likelihood
    l = l - lnZ 
    #l = l -.5*sum(t(y[,,timeInds[k]+1]-pr)%*%Q%*%(y[,,timeInds[k]+1]-pr),na.rm=T) #not sure about this...
    l = l - .5*sum(dysdt*(dysdt%*%Q),na.rm=T)
  }
  if(ave) l = l/Nobs
  if(!is.null(l1W)) 
  {
    if(!is.null(frozenW))
    {
      l1penalty = sum(l1W*abs(W[!frozenW]))
    }
    else l1penalty = sum(l1W*abs(W))
    l = l - l1penalty #always after ave
  }
  
  return(l)
}


FitSF.wrapper = function(l)
{
  #l: list of info, must contain:
  #y
  #x
  #dt ('should')
  
  y = l$y
  x = l$x

  #set defaults if no option provided --- skip conditionals for now
  defaults = list()
  defaults$dt = matrix(1,nrow=nrow(y),ncol=dim(y)[3]-1)
  defaults$iter=5
  defaults$FitFun=FitMoments.SF
  defaults$fitOptions = list()
  defaults$ErrorFun=AymErrorMoments.SF
  defaults$errorFunOptions = list()
  defaults$estimateErrors=T
  defaults$ImputeFirstFun = Impute
  defaults$imputeFirstFunOptions = list()
  defaults$imputeFirst=T              
  defaults$ImputeEMFun = ImputeMean.SF
  defaults$imputeEMFunOptions = list()
  defaults$imputeMean=T
  defaults$imputeAgeCols = c("Age","age","t","time")
  defaults$imputeDead=T
  defaults$imputeCensored=T
  defaults$age = NULL
  defaults$s = NULL
  defaults$where = is.na(y)
  defaults$Npc = NULL
  defaults$imputeInPCSpace = F
  defaults$pcaIndex=1
  defaults$center=!is.null(l$Npc)
  defaults$scale=F
  defaults$addIntercept=T
  defaults$fixedQ = diag(1,ncol(y),ncol(y))
  defaults$diagonalW = F 
  
  #add missing options to l
  missing = setdiff(names(defaults),names(l))
  for (i in 1:length(missing)) l[[missing[i]]] = defaults[[missing[i]]]
  
  #add conditionals
  if(is.null(l[['imputeAge']])) l[['imputeAge']] = (l[["imputeMean"]] | l[["imputeFirst"]])
  
  #you must index by name because of alias bullshit
  fit = FitSF(y=y,
              x=x, #covariates e.g. sex and age
              dt = l[["dt"]],
              iter = l[["iter"]], 
              FitFun = l[["FitFun"]],
              fitOptions = l[["fitOptions"]], 
              ErrorFun = l[["ErrorFun"]],
              errorFunOptions = l[["errorFunOptions"]],
              estimateErrors=l[["estimateErrors"]],
              ImputeFirstFun = l[["ImputeFirstFun"]], 
              imputeFirstFunOptions = l[["imputeFirstFunOptions"]],
              imputeFirst = l[["imputeFirst"]],
              ImputeEMFun = l[["ImputeEMFun"]],
              imputeEMFunOptions = l[["imputeEMFunOptions"]],
              imputeMean=l[["imputeMean"]],
              imputeAge=l[["imputeAge"]],
              imputeAgeCols = l[["imputeAgeCols"]],
              imputeDead=l[["imputeDead"]],
              imputeCensored=l[["imputeCensored"]],
              age = l[["age"]],
              s = l[["s"]],
              where = l[["where"]],
              Npc = l[["Npc"]],
              imputeInPCSpace = l[["imputeInPCSpace"]],
              pcaIndex=l[["pcaIndex"]],
              center=l[["center"]],
              scale=l[["scale"]],
              addIntercept=l[["addIntercept"]],
              fixedQ = l[["fixedQ"]],
              diagonalW = l[["diagonalW"]]
  )
  
  return(fit)
  
}


Natural.SF = function(sf,realOnly=T) #converts fit to natural variables
{
  
  eig = eigen(sf$W)
  l = sort.list(eig[[1]],decreasing=T)
  eig[[1]] = eig[[1]][l]
  eig[[2]] = eig[[2]][,l]
  #sort eigenvalues
  
  P = eig[[2]]
  colnames(P) = sprintf("z%02d",1:nrow(P))
  rownames(P) = colnames(sf$y)
  Pinv = solve(P)
  rownames(Pinv) = sprintf("z%02d",1:nrow(Pinv))
  colnames(Pinv) = colnames(sf$y)
  lambdaz = Pinv%*%sf$lambda
  dlambdaz = sqrt((Pinv)^2%*%(sf$dlambda^2))
  Qz = Pinv%*%sf$Q%*%t(Pinv)
  Wz = Pinv%*%sf$W%*%P
  dWz = sqrt(Pinv^2%*%sf$dW^2%*%P^2)

  z = sf$yimp ######## y or yimp? ######
  for (kk in 1:dim(z)[3]) z[,,kk] = z[,,kk]%*%t(Pinv)
  colnames(z) = sprintf("z%d",1:ncol(z))
  
  
  if(is.null(sf$yimpsd)) zsd = NULL
  else
  {
    zsd = sf$yimpsd
    for (kk in 1:dim(zsd)[3]) zsd[,,kk] = sqrt(zsd[,,kk]^2%*%t(Pinv^2))
    colnames(zsd) = sprintf("z%d",1:ncol(zsd))
  }
  
  if(realOnly)
  {
    z = Re(z)
    if(!is.null(zsd)) zsd = Re(zsd)
    eig[[1]] = Re(eig[[1]])
    eig[[2]] = Re(eig[[2]])
    P = Re(P)
    Pinv = Re(Pinv)
    Wz = Re(Wz)
    lambdaz = Re(lambdaz)
    Qz = Re(Qz)
    dWz = Re(dWz)
    dlambdaz = Re(dlambdaz)
  }
  
  nat = list(z=z,zsd=zsd,eig=eig,P=P,Pinv=Pinv,W=Wz,lambda=lambdaz,Q=Qz,
             dW=dWz,dlambda=dlambdaz)
  return(nat)
}

Natural = function(bs,realOnly=T)
{
  return(Natural.BS(bs,realOnly=realOnly))
}

Natural.BS = function(bs,realOnly=T) #converts fit to natural variables
{
  
  W = bs[["meansd"]][["W"]]
  lambda = bs[["meansd"]][["lambda"]]
  Q = bs[["meansd"]][["Q"]]
  eig = eigen(W[[1]])
  l = sort.list(eig[[1]],decreasing=T)
  eig[[1]] = eig[[1]][l]
  eig[[2]] = eig[[2]][,l]
  #sort eigenvalues
  
  P = eig[[2]]
  colnames(P) = sprintf("z%02d",1:nrow(P))
  rownames(P) = colnames(sf$y)
  Pinv = solve(P)
  rownames(Pinv) = sprintf("z%02d",1:nrow(Pinv))
  colnames(Pinv) = colnames(sf$y)
  lambdaz = Pinv%*%lambda[[1]]
  dlambdaz = sqrt((Pinv)^2%*%(lambda[[2]]^2))
  Qz = Pinv%*%Q[[1]]%*%t(Pinv)
  dQz = sqrt(Pinv^2%*%Q[[2]]%*%t(Pinv)^2)
  Wz = Pinv%*%W[[1]]%*%P
  dWz = sqrt(Pinv^2%*%W[[2]]^2%*%P^2)
  
  
  z = sf$yimp ######## y or yimp? ######
  for (kk in 1:dim(z)[3]) z[,,kk] = z[,,kk]%*%t(Pinv)
  colnames(z) = sprintf("z%d",1:ncol(z))
  
  
  if(is.null(sf$yimpsd)) zsd = NULL
  else
  {
    zsd = sf$yimpsd
    for (kk in 1:dim(zsd)[3]) zsd[,,kk] = sqrt(zsd[,,kk]^2%*%t(Pinv^2))
    colnames(zsd) = sprintf("z%d",1:ncol(zsd))   
  }
  
  
  if(realOnly)
  {
    z = Re(z)
    if(!is.null(zsd)) zsd = Re(zsd)
    eig[[1]] = Re(eig[[1]])
    eig[[2]] = Re(eig[[2]])
    P = Re(P)
    Pinv = Re(Pinv)
    Wz = Re(Wz)
    lambdaz = Re(lambdaz)
    Qz = Re(Qz)
    dWz = Re(dWz)
    dlambdaz = Re(dlambdaz)
  }
  
  
  
  nat = list(z=z,zsd=zsd,eig=eig,P=P,Pinv=Pinv,W=Wz,lambda=lambdaz,Q=Qz,
             dW=dWz,dlambda=dlambdaz,dQ=dQz)
  return(nat)
}

Sim.SF = function(N=1000,
                 W=matrix(rnorm(10*10,0,.1),nrow=10,ncol=10),
                 times = 8, #number of timepoints
                 mu0 = rep(0,ncol(W)),
                 dt = outer(rep(1,N),rep(1,times)),
                 x = array(0,dim=c(N,1,times)), #timepoints are inbetween t+dt/2 #covariates
                 lambda = matrix(rnorm(ncol(W)*ncol(x)),nrow=ncol(W),ncol=ncol(x)), #covariate coefficients
                 rate = rep(1,N), #individual rate
                 mu0i = rep(0,N), #individual mean/offset
                 noiseCov = diag(1,10,10),
                 signal_to_noise_ratio=10,
                 y0=NULL,
                 noNoise=F,
                 qterm=NULL,
                 transform=T, #are you working in the diagonal space or do you need to transform to get there? (quadratic term only)
                 trans=NULL, #needed for qterm
                 transinv=NULL,
                 interventionTime=-1, #add at this time index (always adds after everything else)
                 intervention=rep(0,nrow(W)) #add this to y
)
{
  #generate data using model:
  # y_{n+1}-y_{n} = dt*W(y_n-mu_{in})+eps
  #mu_{in} = mu_0 + x%*%lambda
  #lambda is covariate regression coefficients
  
  
  p = ncol(W)
  
  if (p != ncol(noiseCov)) stop("noise covariance and W not same size!")
  if (ncol(noiseCov)!=nrow(noiseCov)) stop("noise covariance must be square")
  if(times < 2) stop("Need at least 2 timepoints!")
  
  y = array(NA,dim=c(N,p,times))
  if(is.null(y0))
  {
    if(noNoise) eps = 0
    else
    {
      eps = mvrnorm(N,mu=rep(0,p),Sigma = noiseCov)
      for (j in 1:ncol(eps)) eps[,j] = eps[,j]
    }
    y[,,1] = eps*signal_to_noise_ratio
  }
  else 
  {
    colnames(y) = colnames(y0)
    y[,,1] = y0
  }
  
  if(interventionTime==1) for (j in 1:ncol(y)) y[,j,1] = y[,j,1] + intervention[j]
  
  for (i in 2:times)
  {
    odtrate = outer(dt[,i-1]*rate,rep(1,p))
    
    if(noNoise) eps = 0
    else
    {
      eps = mvrnorm(N,mu=rep(0,p),Sigma = noiseCov)
      for (j in 1:ncol(eps)) eps[,j] = eps[,j]*sqrt(abs(dt[,i-1]))
    }
    
    yc = y[,,i-1] -  x[,,i-1]%*%t(lambda)
    for (j in 1:ncol(yc)) yc[,j] = yc[,j] - mu0[j] - mu0i   
    
    y[,,i] = y[,,i-1]+odtrate*(yc%*%t(W)) + eps
    
    if(!is.null(qterm)) 
    {
      if(transform)
      {
        z2 = (as.matrix(y[,,i-1])%*%t(trans))^2
        y2 = z2
      }
      else y2 = y[,,i-1]^2
      y2 = y2%*%t(qterm)
      y[,,i] = y[,,i]+odtrate^2*y2
    }
    
    if(interventionTime==i) for (j in 1:ncol(y)) y[,j,i] = y[,j,i] + intervention[j]
  }
  if(is.null(colnames(y))) colnames(y)=sprintf("y%d",1:ncol(y))
  
  results = list(y=y,x=x,mu0=mu0,p=p,N=N,W=W,dt=dt,times=times,noiseCov=noiseCov,lambda=lambda,mu0i=mu0i)
  
  return(results)
}

################################# helper functions ######################################
################################# helper functions ######################################
################################# helper functions ######################################
Preprocess = function(y,
                      x,
                      normVar=colnames(y), #variables to normalize
                      sexSpecific=T, #do you want to use a sex-specific mean?
                      normIndex=1, #which index do you want to use for normalization (default: first/baseline)
                      sexCol="sex" #sex column name
    )
{
  if(!(sexCol%in%colnames(x)))
  {
    stop("Sex column not found in x, either change column name parameter (sexCol) or do not preprocess by sex (sexSpecific=F).")
  }
  
  data_sd = rep(NA,length(normVar))
  names(data_sd) = normVar
  mu_sex = matrix(NA,nrow=2,ncol=length(normVar))
  rownames(mu_sex) = c("sex 0","sex 1")
  colnames(mu_sex) = normVar

  ypp = y
  for (nv in normVar)
  {
    #lm for sex-specific mean?
    if(sexSpecific)
    {
      mdata = data.frame(sex=c(x[,sexCol,normIndex]),y=c(y[,nv,normIndex]))
      m = lm(y~sex,mdata)
      mu_sex[1,nv] = predict(m,newdata=data.frame(sex=0,y=NA))
      mu_sex[2,nv] = predict(m,newdata=data.frame(sex=1,y=NA))
      data_sd[nv] = sd(y[,nv,normIndex]-predict(m,newdata=mdata),na.rm=T)
      
      newdata = data.frame(sex=c(x[,sexCol,]),y=c(ypp[,nv,]))
      mu = predict(m,newdata=newdata)
      mu = array(mu,dim=c(nrow(ypp),ncol=dim(ypp)[3]))
      ypp[,nv,] = (ypp[,nv,]-mu)
      ypp[,nv,] = ypp[,nv,]/data_sd[nv]
      rm(mu)
    }
    else 
    {
      mu_sex[,nv] = mean(y[,nv,normIndex],na.rm=T)
      mu = array(mean(y[,nv,normIndex],na.rm=T),dim=c(nrow(ypp),ncol=dim(ypp)[3]))
      data_sd[nv] = sd(y[,nv,normIndex],na.rm=T)
      ypp[,nv,] = (ypp[,nv,]-data_mu[nv])/data_sd[nv]
    }
  }
  
  return(list(ypp=ypp,mu=mu_sex,sd=data_sd))
  
}


cov2 = function(x,
                y=NULL,
                use="pairwise.complete",
                center=T,
                cond = NULL, #variables to condition on, should be matrix/dataframe
                ...
)
{

  if (grepl("complete",use)) na.rm=T
  else na.rm=F
  
  if(!is.null(cond))
  {
    #matrix cookbook notation
    #A: col(x)
    #B: col(cond)
    warning("Conditional covariance needs validation. Consider doing a binary-specific solution.")
    if(!is.null(y)) stop("y cannot be supplied when conditioning.")
    if(!center)    
    {
      mua = 0 #apply(x,2,mean,na.rm=na.rm) #not centering, so assuming mean of a (x) is 0
      xb = cond
      for (i in 1:ncol(xb)) xb[,i] = cond[,i]-mean(cond[,i],na.rm=na.rm)
      
      Cb = cov(cond,use=use,...)
      Cc = cov(x,cond,use=use,...)
      Qb = solve(Cb)
      
      xa = x
      for (i in 1:nrow(xa)) 
      {
        xa[i,] = xa[i,] - (mua + Cc%*%Qb%*%xb[i,])[,1]
      }
      
      Ccond = cov2(xa,use=use,center=F,cond=NULL)
      
      rownames(Ccond) = colnames(x)
      colnames(Ccond) = colnames(x)
      
      return(Ccond) 
    }
    else #ignore mean since we're implicitly going to subtract it anyways
    {
      if(!is.null(y)) stop("no cond + y implemented")
      Ca = cov(x,use=use,...)
      Cb = cov(cond,use=use,...)
      Cc = cov(x,cond,use=use,...)
      
      
      Qb = solve(Cb)
      
      Ccond = Ca - Cc%*%Qb%*%t(Cc)
      
      rownames(Ccond) = colnames(x)
      colnames(Ccond) = colnames(x)
      
      
      return(Ccond) 
    }
  }
  else if (center)
  {
    return(cov(x,y,use=use,...))
  }
  if(is.null(y)) y = x
  
  
  C = cov(x,y,use=use,...)
  
  if(!center)
  {
    mux = apply(x,2,mean,na.rm=na.rm)
    muy = apply(y,2,mean,na.rm=na.rm)
    C = C + outer(mux,muy)
  }

  return(C)
}

prcomp2 = function(x,
                   center=T,
                   scale.=F,
                   use = "pairwise.complete.obs", #sent to cov
                   eps=0, #prevents negative eigenvalues
                   meanImpute = T, #don't change this unless you're sure what you're doing
                   ... #send to cov #cond can go there
)
{
  #purpose is to do prcomp() but allow available case analysis (i.e. missing values are okay)
  #prcomp uses svd (supposedly better numerics)
  #here I use cov (handles missing data 'better')
  #uses the eigendecomposition of the covariance matrix instead of SVD
  
  for (i in 1:ncol(x))
  {
    if(center) x[,i] = x[,i]-mean(x[,i],na.rm=T)
    if(scale.) x[,i] = x[,i]/sd(x[,i],na.rm=T)
  }
  
  if(ncol(x)==1)
  {
    warning("prcomp2 Only one variable provided!")
    rotation = matrix(1,1,1)
    colnames(rotation)=sprintf("PC%d",1:ncol(rotation))
    rownames(rotation)=colnames(x)
    pc = as.matrix(x)
    
    l = list(sdev=sqrt(sd(x[,1])+eps),rotation=rotation,center=center,scale=scale.,x=pc)
    return(l)
  }
  
  
  if(!center) C = cov2(x,use=use,center=F,...)  #already centered/or not
  else C = cov2(x,use=use,...)
  e = eigen(C,symmetric=T)
  rotation = e[[2]]
  colnames(rotation)=sprintf("PC%d",1:ncol(rotation))
  rownames(rotation)=colnames(x)
  
  
  #imp mean for x
  ximp = x
  if(meanImpute)
  {
    for (i in 1:ncol(x))
    {
      ximp[is.na(x[,i]),i] = mean(x[,i],na.rm=T)
    }
    if(sum(is.na(ximp))>0) stop("NAs found in ximp")
  }
  
  pc = as.matrix(ximp)%*%rotation

  
  l = list(sdev=sqrt(e[[1]]+eps),rotation=rotation,center=center,scale=scale.,x=pc,eigen=e)
  return(l)
}

ArrayToStartStop = function(age, #2d array of ages, dimensions individuals x times
                            y, #3d array of time-dependent covariates, individuals x variables x times
                            s, #survival information: start, stop and event
                            dropNAMeasurement=T #drops all of the info from start to first measurement (IF there's a gap)
)
{
  #convert 3d array, y, into start-stop format
  #age, y and s show all be aligned
  
  #doesn't work for time-to-event data that doesn't cause censorship e.g. dementia
  
  #notes:
  #doesn't extrapolate backwards
  #e.g. paquid data has start times < first timepoint, so there are no values for the period from start time to first measurement
  
  if(nrow(age)!=nrow(y)) stop("age and y must have same number of individuals")
  if(nrow(age)!=nrow(s)) stop("age and s must have same number of individuals")
  if(ncol(age)!=dim(y)[3]) 
  {
    print("age:")
    print(str(age))
    print("y:")
    print(str(y))
    stop("age and y must have same number of timepoints")
  }
  
  #steps:
  #1. convert to long format
  #2. add identifier for individual and measurement number
  #3. add to tmerge
  
  ydf = list()
  for (i in 1:ncol(age))
  {
    df = data.frame(y[,,i])
    colnames(df) = colnames(y)
    df[,"age"] = age[,i]
    df[,"id"] = 1:nrow(age)
    df[,"measurement"] = i
    df[,"start"] = s[,"start"]
    df[,"death"] = s[,"stop"]
    df[,"status"] = s[,"status"]
    
    ydf[[i]] = df
  }
  ydf = do.call(rbind,ydf)
  
  #print(head(ydf))
  #print(ydf[sort.list(ydf[,"id"]),][1:10,])
  
  #start-stop format
  stst = tmerge(subset(ydf,measurement==1)[,c("id"),drop=F],
                data2=subset(ydf,measurement==1),tstop=death,tstart=start,id=id)
  
  stst = tmerge(stst,ydf,id=id,status=event(death,status),measurement=tdc(age,measurement))
  
  
  #add predictors/covariates on
  stst[,colnames(y)] = NA
  for (i in 1:nrow(stst))
  {
    logi = (ydf[,"id"] == stst[i,"id"]) & (ydf[,"measurement"] == stst[i,"measurement"])
    logi[is.na(logi)] = F
    
    if(sum(logi) < 1) next #can happen if we have survival info from before first measurement
    if(sum(logi) > 1) stop("multiple measurements with same id and measurement")
    stst[i,colnames(y)] = ydf[logi,colnames(y)]
  }
  
  if(dropNAMeasurement) stst = stst[!is.na(stst[,"measurement"]),]
  
  return(stst)  
}

ListMeanSD = function(l,sem=F,negToZero=F,na.rm=F,skipsd=F)
{
  #calculates SD of list of objects of same size (e.g. matrices/arrays)
  #preserves dimensions!
  #returns mean and sd or sem
  #negToZero: force negative variance to 0 (happens if mu^2 ~= v)
  
  if(length(l)==1) 
  {
    if(is.null(dim(l[[1]]))) 
    {
      mu = l[[1]]
      s = rep(0,length(l))
      if(!is.null(names(mu)))  names(s) = names(mu)
      return(list(mean=mu,sd=s))
    }
    else 
    {
      mu = l[[1]]
      s = array(0,dim(l[[1]]))
      if(!is.null(colnames(mu))) colnames(s) = colnames(mu)
      return(list(mean=mu,sd=s))
    }
  }
  else if (length(l)<1) return(list(mean=NA,sd=NA))
  else if(na.rm)
  {
    lnotNA = l
    for (i in 1:length(l)) 
    {
      lnotNA[[i]] = !is.na(l[[i]])
      l[[i]][is.na(l[[i]])] = 0
    }
    numNotNA = Reduce("+",lnotNA)
    
    mu = Reduce("+",l)/numNotNA
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/numNotNA
    v = v - mu^2
    v = v*numNotNA/(numNotNA-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(numNotNA)
    return(list(mean=mu,sd=s))
  }
  else
  {
    mu = Reduce("+",l)/length(l)
    l = lapply(l,function(x){return(x^2)})
    
    if(skipsd) return(list(mean=mu,sd=NULL))
    
    v = Reduce("+",l)/length(l)
    v = v - mu^2
    v = v*length(l)/(length(l)-1) #convert to unbiased estimate
    
    s = sqrt(v)
    
    if(negToZero)
    {
      v[is.nan(v)]=NA
      s[v<0] = 0
    }
    
    if(sem) s = s/sqrt(length(l))
    return(list(mean=mu,sd=s))
  }
  
}

PWLinearInterp = function(age, #vector of times #should be sorted
                          minAge = NA,
                          maxAge = NA
)
{
  #puts upper and lower bounds on each element then interpolates using a broken-stick style model
  #that is, piecewise linear between lower and upper bound (linear function of number of steps)
  
  nonNA = sum(!is.na(c(age,minAge,maxAge)))
  if(nonNA < 1.9) 
  {
    print(nonNA)
    warning("Not enough information to impute")
    return(age) #nothing to do! 
  }
  #where do we get values before/after this age?
  lb = age[1]
  if(is.na(lb)) lb = minAge
  #figure out upper and lower bound for each element
  j = 1
  missing = integer()
  while (j <= length(age))
  {
    #print(j)
    #print(missing)
    if(is.na(age[j])) #increment missing list
    {
      missing = c(missing,j)
    }
    else #interpolate and raise lower bound
    {
      #interpolate if there is a gap
      if(length(missing)> 0)
      {
        ub = age[j] 
        age[missing] = lb + (ub-lb)/(length(missing)+1)*(missing-missing[1]+1)
      }
      #increment lower bound
      lb = age[j]
      
      #empty missing list in preparation for next
      missing = integer()
    }
    j = j + 1
  }
  
  #interpolate if there is a gap
  if(length(missing)> 0)
  {
    ub = maxAge
    age[missing] = lb + (ub-lb)/(length(missing)+1)*(missing-missing[1]+1)
  }
  
  return(age)
}

#reshape into 3d array
Reshape = function(data, #long format dataframe (multiple entries for each individual)
                   ycols, #predictor columns to put into y
                   xcols, #covariate columns to put into x
                   doSurv=F, #also create survival object?
                   startCol = "baseline_age", #set to NULL to default to 0 #age of start time column (for survival)
                   stopCol = "death_age", 
                   statusCol = "status",
                   idCol="id", #unique identifier
                   ageCol="age", #name of age column
                   stationary=c("sex"), # will be hard imputed across all timepoints; must be in xcols
                   preprocess=T, #apply preprocessing to y after reshaping?
                   ... #sent to Preprocess()
)
{    
  #converts from long (?) to 3d arrays
  #rearrange so id column is first
  data = data[,c(idCol,setdiff(colnames(data),idCol))]
  
  #uses id frequency to figure out max number of times (measurements)
  tab = table(data[,idCol])
  #print(tab)
  n = length(tab)
  print(sprintf("N: %d",n))
  m = max(tab)
  print(sprintf("timepoints: %d",m))
  #print(tab)
  y = array(NA,dim=c(n,ncol(data),m))
  dimnames(y) = list(names(tab),colnames(data),sprintf("t%02d",1:dim(y)[3]))
  
  #print("filling y...")
  temp = data
  temp[,"timepoint"] = 1
  un = names(tab)
  for (i in 1:length(un)) #loop through each unique id
  {
    logi = data[,idCol] == un[i] 
    #sort.list = rank in many cases e.g. pre-sorted
    temp[logi,"timepoint"] = rank(data[logi,ageCol]) 
    for (j in 1:sum(logi)) #go through each row
    {
      tn = temp[logi,"timepoint"][j]
      #print("converting")
      y[i,,tn] = as.numeric(data[logi,][j,])
    }
  }
  rm(temp)    
  
  data_arr = y
  
  #sort
  #print("sorting")
  #print(tab[is.na(as.character(sort(as.numeric(names(tab)))))])
  #print(as.character(sort(as.numeric(rownames(data_arr))))[is.na(as.character(sort(as.numeric(rownames(data_arr)))))])
  #print(dim(data_arr))
  #print(colnames(data_arr))
  data_arr = data_arr[as.character(sort(as.numeric(rownames(data_arr)))),,]
  #print(dim(data_arr))
  
  #print("orphans:")
  #print(setdiff(ycols,colnames(data_arr)))
  
  y = data_arr[,ycols,]
  
  x = data_arr[,xcols,,drop=F]
  
  #survival
  #all entry age are 0 ( I think)
  if(doSurv) 
  {
    if(is.null(startCol)) s = Surv(time=rep(0,nrow(data_arr)),time2=data_arr[,stopCol,1],event=data_arr[,statusCol,1])
    else s = Surv(time=data_arr[,startCol,1],time2=data_arr[,stopCol,1],event=data_arr[,statusCol,1])
  }
  else s = NULL
  
  #hard imputations / carry forward
  for (j in 1:length(stationary))
  {
    sanity = apply(x[,stationary[j],],1,sd,na.rm=T) > 0
    sanity[is.na(sanity)] = F
    if(any(sanity))
    {
      print(sprintf("found non-stationary values in %s",stationary[j]))
      print(x[sanity,stationary[j],])
      
      print("rows violating stationarity:")
      print(1:nrow(x)[sanity])
      
      stop("sanity check failed")
    }
    else
    {
      for (k in 2:dim(x)[3]) x[is.na(x[,stationary[j],k]),stationary[j],k] = x[is.na(x[,stationary[j],k]),stationary[j],1]
    }
  }
  
  print("random spot checks (ages should be increasing):")
  print(x[sample(1:nrow(x),1),ageCol,])
  print(x[sample(1:nrow(x),1),ageCol,])
  print(x[sample(1:nrow(x),1),ageCol,])
  
  
  #imputing dt:
  #this can be tricky if we have lots of irregularely spaced data or one individual has more timepoints than the others
  #e.g. suppose an individual is measured at two closely-spaced timepoints, then an imputed age could be larger than the age at the next timepoint
  #interval censorship is hard to deal with
  #e.g. panel 0 isn't measured for an individual then they're censored periodically
  #hard impute
  #new version:
  minAge = min(c(x[,ageCol,]),na.rm=T)*.95    #add a bit of space to prevent dt = 0
  maxAge = max(c(x[,ageCol,]),na.rm=T)*1.05   #add a bit of space to prevent dt = 0
  age = t(apply(x[,ageCol,],1,PWLinearInterp,minAge=minAge,maxAge=maxAge))
  print("sanity check on imputed ages, should be 0:")
  print(sum(abs(x[,ageCol,]-age),na.rm=T))
  if(sum(abs(x[,ageCol,]-age),na.rm=T)>0.01) stop("Age sanity check failed")
  x[,ageCol,] = age
  dt = x[,ageCol,-1]-x[,ageCol,-dim(x)[3]]
  
  #hard imput dt # depricated
  #for (i in 1:ncol(dt)) dt[is.na(dt[,i]),i] = median(dt[,i],na.rm=T)
  
  
  #now hard impute age
  #forward impute:
  for (k in 2:dim(x)[3]) x[is.na(x[,ageCol,k]),ageCol,k] = x[is.na(x[,ageCol,k]),ageCol,k-1] + dt[is.na(x[,ageCol,k]),k-1] 
  #backwards impute: #shouldn't be necessary
  #for (k in 1:(dim(x)[3]-1)) x[is.na(x[,ageCol,k]),ageCol,k-1] = x[is.na(x[,ageCol,k]),ageCol,k] - dt[is.na(x[,ageCol,k]),k] 
  
  #sanity check that imputation is self-consistent
  if(!isTRUE(all.equal(dt,x[,ageCol,-1]-x[,ageCol,-dim(y)[3]]))) 
  {
    print(all.equal(dt,x[,ageCol,-1]-x[,ageCol,-dim(y)[3]]))
    print("different:")
    print(apply(abs(dt-x[,ageCol,-1]+x[,ageCol,-dim(y)[3]]) > 0.001,2,sum,na.rm=T))
    print(data_arr[apply(abs(dt-x[,ageCol,-1]+x[,ageCol,-dim(y)[3]]) > 0.001,1,any,na.rm=T),,1])
    print("panels:")
    print(table(data_arr[apply(abs(dt-x[,ageCol,-1]+x[,ageCol,-dim(y)[3]]) > 0.001,1,any,na.rm=T),"panel",1]))
    print("NAs:")
    print(sum(is.na(dt)))
    print(sum(is.na(x[,ageCol,])))
    print(x[apply(is.na(x[,ageCol,]),1,any),,1])
    print(y[apply(is.na(x[,ageCol,]),1,any),,1])
    print(which(apply(is.na(x[,ageCol,]),1,any)))
    stop("problem with age and dt")
  }
  
  if(preprocess)
  {
    ypp = Preprocess(y=y,x=x,...)
  }
  else ypp = list()
  
  l = list(y=y,x=x,dt=dt,s=s)
  l = c(l,ypp)
  
  return(l)
}

################################# overloaded functions ######################################
################################# overloaded functions ######################################
################################# overloaded functions ######################################
#giving me an error
#note: use setMethod and methods to add new methods (overloads existing functions like print, summary and plot)
#you must first add the new class to the class list
#setOldClass("SF") #add SF to class list
#print.SF = function(x)
#{
#  print(x)
#}
#setMethod(print,signature="SF",definition=print.SF)
#summary.SF = function(object)
#{
#  print(x[["iter"]])
#}
#setMethod(summary,signature=signature(object="SF"),definition=summary.SF)
#plot.SF = function(x)
#{
#  library(ggplot)
#  TilePlot(x[["W"]])
#}
#setMethod(plot,signature(x="SF"),definition=plot.SF)
