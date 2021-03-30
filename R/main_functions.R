#' @export
get_orthoMOD = function(X,siglev=1e-04,numSearch=100) {
  # X = p by n matrix
  X.pca = pca.hy(X,subt.mean=F)
  init.d = length(which(X.pca$eigenval>1e-10))
  cutoff = sqrt(qchisq(p=(1-siglev),df=init.d))
  MODout = matrix(0,nrow=nrow(X),ncol=init.d)

  curdim = init.d
  curX = X
  while (curdim > 1) {
    MODres = find_MOD(X=curX,siglev=siglev,numSearch=numSearch,minmax=F)

    if (length(MODres$outliers)>0) {
      MODout[,(init.d - curdim + 1)] = MODres$MOdir
      curX = rm_onedim(X=curX,u=MODres$MOdir)$resd.X
      curdim = curdim - 1
    } else {
      MODout[,(init.d - curdim + 1):init.d] = MODres$MOdir
      curdim = 0
    }
    rm(MODres)
    cat((init.d - curdim),"\n")
  }
  if (curdim==1) {
    MODout[,init.d] = pca.hy(curX,subt.mean=F)$dirmat[,1]
  }
  MODprojscore = t(apply(t(MODout)%*%X,1,FUN=function(x){pd.rate.hy(x,qrsc=T)}))
  rownames(MODprojscore) = paste("orthoMOD",1:nrow(MODprojscore))
  colnames(MODout) = paste("orthoMOD",1:nrow(MODprojscore))
  outliers = which(apply(abs(MODprojscore),2,max)>cutoff)
  return(list(orthoMOD=MODout,MOscore=MODprojscore,outliers=outliers,cutoff=cutoff))
}

## One-layer function
#' @export
find_MOD = function(X,siglev=1e-04,numSearch=300,qrsc=T,minmax=F) {
  # One-layer function for finding MOD
  # Input: X = PC scores (p by n)

  X.pca = pca.hy(X,subt.mean=F)
  init.d = length(which(X.pca$eigenval>1e-10))
  cutoff = sqrt(qchisq(p=(1-siglev),df=init.d))

  if (init.d == 1) {
    MOdirection = X.pca$dirmat[,1]
    MOoutliers = which(abs(pd.rate.hy(x=X.pca$projmat[1,],qrsc=qrsc))>cutoff)
  } else {
    original.X = X; rm(X);
    original.PC = X.pca$dirmat[,1:init.d]

    X0 = X.pca$projmat[1:init.d,]
    stdX0 = t(apply(X0,1,FUN=function(x){pd.rate.hy(x,qrsc=qrsc)}))
    X0.pca = pca.hy(X0,subt.mean=F)

    X0.dirs = get_outdirs(X=X0,numSearch=numSearch)
    X0.POmat = t(apply(X0.dirs%*%X0,MARGIN=1,FUN=function(x){pd.rate.hy(x,qrsc=qrsc)})) # n by length(indir)
    indexmax0 = which.max(apply(abs(X0.POmat),1,max))
    PO1maxout = which.max(abs(X0.POmat[indexmax0,]))
    # X.POmat = t(apply(X=A%*%X,MARGIN=1,get_stdvar)) # n by length(indir)
    # X.ADstat = apply(X.POmat,1,ADstatWins.hy) #
    X0index = which(abs(stdX0[,PO1maxout])>qnorm(0.8))
    X = X0
    X[which(! c(1:nrow(X0)) %in% X0index),] = 0
    X.dirs = get_outdirs(X=X,numSearch=numSearch)
    X.POmat = t(apply(X.dirs%*%X,MARGIN=1,FUN=function(x){pd.rate.hy(x,qrsc=qrsc)}))
    indexmax0 = which.max(abs(X.POmat[,PO1maxout]))
    # max(apply(abs(X.POmat),1,max))

    if (abs(X.POmat[indexmax0,PO1maxout])<cutoff) {
      cat("No outlier signals left in data")
      MOdirection = original.PC%*%X0.pca$dirmat[,1:init.d]
      MOoutliers = c()
    } else {
      indexset = which(abs(X.POmat)[,PO1maxout]>cutoff)
      length(indexset)
      # plot_pileup(Pileup=ScissorOutput$logData,cases=PO1maxout,
      #             Ranges=Ranges,col.pileup="black")
      # plot_pileup(Pileup=PCdir%*%X,cases=PO1maxout,
      #             main=paste("low-dim"),
      #             Ranges=Ranges,col.pileup="black")
      # plot_pileup(Pileup=PCdir%*%original.PC%*%X.dirs[indexmax0,],ylim=c(-0.08,0.08),
      #             main=paste("original MOD"),
      #             Ranges=Ranges,col.pileup="black")
      # kdeplot.hy(X.POmat[indexmax0,],indlist=which(abs(X.POmat[indexmax0,])>cutoff),text=T)
      # abline(v=c(-cutoff,cutoff),col="red")

      # Generates orthogonal directions
      nindex = length(indexset)
      if (nindex==1) {
        cat("Only one PC selected","\n")
        MOdirection = original.PC%*%X.dirs[indexmax0,]
        MOoutliers = which(abs(X.POmat[indexmax0,])>cutoff)
      } else {
        X1.dirs = X.dirs[indexset,]
        X1.POmat = X.POmat[indexset,]
        X2.POmaxmean.res = apply(X1.dirs,1,FUN=function(x){get_X2POfromX1maxout_v2(X1dir=x,X=X,PO1maxout=PO1maxout,
                                                                              numSearch=200,cutoff=cutoff,qrsc=qrsc)})
        X2.POmaxmean = X2.POmaxmean.res[1,]
        X2.POmaxmax = X2.POmaxmean.res[2,]
        min(X2.POmaxmean)
        min(X2.POmaxmax)

        if (min(X2.POmaxmax)>cutoff) {
          cat("At least two outlier signals detected","\n")
          # Both directions are outlier signals
          indexmax = which.min(X2.POmaxmean)
          MOdirection = original.PC%*%X1.dirs[indexmax,]
          MOoutliers = which(abs(X1.POmat[indexmax,])>cutoff)

          # plot_pileup(Pileup=PCdir%*%original.PC%*%X1.dirs[indexmax,],ylim=c(-0.08,0.08),
          #             main=paste("orthogonal MOD"),
          #             Ranges=Ranges,col.pileup="black")
          # kdeplot.hy(X1.POmat[indexmax,],indlist=MOoutliers,text=T)
          # abline(v=c(-cutoff,cutoff),col="red")
          # for (case in which(abs(X1.POmat[indexmax,])>cutoff)) {
          #   plot_pileup(Pileup=PCdir%*%X,cases=case,
          #               Ranges=Ranges,col.pileup="black")
          # }
        } else {
          cat("Only one outlier signal detected","\n")
          # Direction 2 is a main signal
          indexset2 = which(X2.POmaxmax<cutoff)
          indexmax = indexset2[which.max(abs(X1.POmat[indexset2,PO1maxout]))]
          MOdirection = original.PC%*%X1.dirs[indexmax,]
          MOoutliers = which(abs(X1.POmat[indexmax,])>cutoff)

          # plot_pileup(Pileup=PCdir%*%X1.dirs[indexmax,],ylim=c(-0.08,0.08),
          #             main=paste(curdim,"| orthogonal MOD"),
          #             Ranges=Ranges,col.pileup="black")
          # kdeplot.hy(X1.POmat[indexmax,],indlist=MOoutliers,text=T)
          # abline(v=c(-ScissorOutput$GSCout$cutoff,ScissorOutput$GSCout$cutoff),col="red")
        }
      }
    }
  }
  return(list(MOdir=as.vector(MOdirection),outliers=MOoutliers))
}

## Inner functions
#' @export
get_X2POfromX1maxout_v1 = function(X1dir,X,PO1maxout,cutoff=4,numSearch=100,qrsc=T) {
  # Get the second outlying direction where X1 max outlier is the most outlying.
  # Save only PO
  X_temp = rm_onedim(X=X,u=X1dir)$resd.X
  A_temp = get_outdirs(X=X_temp,numSearch=numSearch)
  if (dim(A_temp)[1]==1) {
    POmat_temp = pd.rate.hy(as.vector(A_temp%*%X_temp),qrsc=T)
    return(POmat_temp)
  } else {
    POmat_temp = t(apply(X=A_temp%*%X_temp,MARGIN=1,FUN=function(x){pd.rate.hy(x,qrsc=qrsc)})) # n by length(indir)
    # ADstat_temp = apply(POmat_temp,1,ADstatWins.hy)
    ## Strategy 1: only PO1maxout
    index_temp = which.max(abs(POmat_temp[,PO1maxout]))
    ## Strategy 2: mean of outliers
    # inner_tmp = function(x) {
    #   tmp_cases = which(x>cutoff)
    #   if (length(tmp_cases)>0) {
    #     return(mean(x[tmp_cases]))
    #   } else {
    #     return(x[PO1maxout])
    #   }
    # }
    # index_temp = which.max(apply(abs(POmat_temp),1,inner_tmp))
    return(POmat_temp[index_temp,])
  }
}

#' @export
get_X2POfromX1maxout_v2 = function(X1dir,X,PO1maxout,cutoff=4,numSearch=100,qrsc=T) {
  # Consider PO2 max mean of the PO1 outliers
  # Save only PO
  PO1std = abs(pd.rate.hy(x=as.vector(t(X1dir)%*%X),qrsc=qrsc))
  PO1out = which(PO1std>cutoff)

  X_temp = rm_onedim(X=X,u=X1dir)$resd.X
  A_temp = get_outdirs(X=X_temp,numSearch=numSearch)
  # plot_pileup(Pileup=PCdir%*%original.PC%*%X_temp,case=152,
  #             Ranges=Ranges,col.pileup="black")
  if (dim(A_temp)[1]==1) {
    POmat_temp = pd.rate.hy(as.vector(A_temp%*%X_temp),qrsc=qrsc)
    max.output = max(abs(POmat_temp))
    mean.output = mean(abs(POmat_temp[PO1out]))
  } else {
    POmat_temp = t(apply(X=A_temp%*%X_temp,MARGIN=1,FUN=function(x){pd.rate.hy(x,qrsc=qrsc)})) # n by length(indir)
    # ADstat_temp = apply(POmat_temp,1,ADstatWins.hy)
    ## Strategy 1: only PO1maxout
    # index_temp = which.max(abs(POmat_temp[,PO1maxout]))
    ## Strategy 2: mean of outliers
    if (length(PO1out)>1) {
      # index_temp = which.max(apply(abs(POmat_temp[,PO1out]),1,mean))
      max.output = max(apply(abs(POmat_temp),2,max))
      mean.output = mean(apply(abs(POmat_temp[,PO1out]),2,max))
    } else {
      max.output = max(abs(POmat_temp))
      mean.output = max(abs(POmat_temp[,PO1out]))
      # index_temp = which.max(abs(POmat_temp[,PO1out]))
    }
  }
  # max.output = max from all data points
  # mean.output = mean from PO1 outliers
  return(c(mean.output,max.output))
}


#' @export
get_outdirs = function(X,numSearch=300) {
  n = ncol(X); M = nrow(X);

  X.pca = pca.hy(X,subt.mean=F)
  d = length(which(X.pca$eigenval>1e-10))
  ndir = numSearch*d
  if (d==1) {
    A = matrix(X.pca$dirmat[,1],nrow=1)
  } else {
    if (d < M) {
      A0 = generdir(t(X.pca$projmat[1:d,]),ndir=ndir)
      A = t(X.pca$dirmat[,1:d]%*%t(A0))
    } else {
      A = generdir(t(X),ndir=ndir) # generates `ndir' directions (ndir by M)
      ## Add more potential directions
      Scov = (X%*%t(X))/n
      B0 = solve(Scov)%*%X   # M by n
      B0 = B0[,-which(apply(B0,2,FUN=function(x){sqrt(sum(x^2))})<1e-10)]
      B1 = sweep(B0,2,apply(B0,2,FUN=function(x){sqrt(sum(x^2))}),"/")  # normalized M by n
      A = rbind(A,t(B1))  # ndir by M
    }
  }
  return(A)
}

#' @export
get_PO2max = function(PO1,PO2,siglev=1e-04,rm.case=NULL,qrsc=T) {
  cutoff = qnorm(1-siglev)

  # i = 555
  # PO1 = out_temp1[,i]
  # PO2 = out_temp2[,i]
  PO1std = abs(pd.rate.hy(x=PO1,qrsc=qrsc))
  PO2std = abs(pd.rate.hy(x=PO2,qrsc=qrsc))
  PO1out = which(PO1std>cutoff)
  # plot(x=X[1,],y=X[2,],col="black",xlim=c(-8,15),ylim=c(-8,15))
  # abline(h=0,v=0,lty=2)
  # points(x=X[1,PO1out],y=X[2,PO1out],col="green",pch=19)
  # abline(a=0,b=A[i,2]/A[i,1],col="skyblue")
  # abline(a=0,b=B[i,2]/B[i,1],col="pink")

  if (length(PO1out)>0) {
    if (length(which(!PO1out %in% rm.case))>0) {
      PO1out = PO1out[which(!PO1out %in% rm.case)]
    }
    if (length(which(PO2std>cutoff))>0) {
      # This is the case when two directions are both outlier signals
      return(max(PO2std[PO1out]))
    } else {
      # Direction 2 is a main signal
      return(0)
    }
  } else {
    # Direction 1 is a main signal
    return(0)
  }
}

#' @export
get_PO2mean = function(PO1,PO2,siglev=1e-04,rm.case=NULL,qrsc=T) {
  cutoff = qnorm(1-siglev)

  # i = 555
  # PO1 = out_temp1[,i]
  # PO2 = out_temp2[,i]
  PO1std = abs(pd.rate.hy(x=PO1,qrsc=qrsc))
  PO2std = abs(pd.rate.hy(x=PO2,qrsc=qrsc))
  PO1out = which(PO1std>cutoff)

  # plot(x=X[1,],y=X[2,],col="black",xlim=c(-8,15),ylim=c(-8,15))
  # abline(h=0,v=0,lty=2)
  # points(x=X[1,PO1out],y=X[2,PO1out],col="green",pch=19)
  # abline(a=0,b=A[i,2]/A[i,1],col="skyblue")
  # abline(a=0,b=B[i,2]/B[i,1],col="pink")

  if (length(PO1out)>0) {
    if (length(which(!PO1out %in% rm.case))>0) {
      PO1out = PO1out[which(!PO1out %in% rm.case)]
    }
    if (length(which(PO2std>cutoff))>0) {
      # This is the case when two directions are both outlier signals
      return(mean(PO2std[PO1out]))
    } else {
      # Direction 2 is a main signal
      return(0)
    }
  } else {
    # Direction 1 is a main signal
    return(0)
  }
}

#' @export
get_POscore = function(X,u) {
  # X: d by n data matrix
  # u: d by 1 vector
  u = u/sqrt(sum(u^2))
  score = t(X)%*%u
  std.score = (score - median(score))/mad(score)
  return(list(score=score,std.score=std.score))
}

#' @export
get_unitdir = function(u) {
  if (sum(u^2)>1e-10) {
    u/sqrt(sum(u^2))
  } else {
    u
  }
}

#' @export
get_stdvar = function(x,robust=T) {
  if (mad(x)>1e-10) {
    (x-median(x))/mad(x)
  } else {
    rep(0,length(x))
  }
}

#' @export
rm_onedim = function(X,u) {
  # X = d by n data
  # u = d by 1 vector
  u = get_unitdir(u)
  pscore = t(u)%*%X
  return(list(appr.X = u%*%pscore,resd.X = (X - u%*%pscore), pscore = as.vector(pscore)))
}

#' @export
rm_subspace = function(X,basis) {
  # X = d by n data
  # basis = d by K basis or d by 1 vector
  if (is.null(dim(basis))) {
    K = 1
    return(rm_onedim(X=X,u=basis)$resd.X)
  } else {
    K = dim(basis)[2]
    for (k in 1:K) {
      X = rm_onedim(X=X,u=basis[,k])$resd.X
    }
    return(X)
  }
}

