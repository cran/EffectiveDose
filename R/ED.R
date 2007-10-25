
setOldClass("locfit")
setClass("locpoly", representation(x="numeric", y="numeric", call="call"))
setClass("monofit", representation(x="numeric",y="numeric",z="matrix"))
setMethod("plot", signature(x="monofit",y="missing"), function(x,y,...){
plot(x@x,x@y,...)
})

setMethod("lines", signature(x="monofit"), function(x,...){
lines(x@x,x@y,...)
})

locpoly<-function(...){
	call<-match.call(KernSmooth:::locpoly)
	ans<-KernSmooth:::locpoly(...)
	ret<-new("locpoly", x=ans$x, y=ans$y, call=call)
return(ret)
}

setMethod("plot", signature(x="locpoly",y="missing"), function(x,y,...){
plot(x@x,x@y,...)
})


setMethod("lines", signature(x="locpoly"), function(x,...){
lines(x@x,x@y,...)
})

"mono.1d"<-function(fit, bandwidth, xx, kernel="epanech", mono1)
{	
	n<-length(fit[[1]])
	h<-bandwidth
	mono1<-match.arg(mono1, c("increasing", "decreasing"))
	if(missing(xx)||class(xx)!="numeric") xx<-fit[[1]]
	xl<-(xx-min(fit[[1]]))/(max(fit[[1]])-min(fit[[1]]))
	yl<-fit[[2]]
	N<-length(xl)
	r<-numeric(N)
	monotoneinverse<-function(tn,kernel="epanech")
	{
		xx<-(yl-tn)/h
		g<-rep(1,n)
		if (kernel == "epanech"){
		if(mono1=="increasing") (sum(1/2-3/4*(xx[abs(xx)<1]-1/3*xx[abs(xx)<1]^3))+sum(g[xx<(-1)]))/n
		else if(mono1=="decreasing") (sum(1/2+3/4*(xx[abs(xx)<1]-1/3*xx[abs(xx)<1]^3))+sum(g[xx>1]))/n
		}
	}
	for(l in 1:N)
	{
		tmax  = ifelse(mono1 == "increasing", max(yl), min(yl))
		tmin  = ifelse(mono1 == "increasing", min(yl), max(yl))
		wertmin<-monotoneinverse(tmin)
		wertmax<-monotoneinverse(tmax)
		if(wertmax==wertmin) r[l]<-tmin
		wertt<-xl[l]
		for(i in 1:150)
		{
			tneu<-tmin+(wertt-wertmin)*(tmax-tmin)/(wertmax-wertmin)
			wertneu<-monotoneinverse(tneu)
			if(abs(wertt-wertneu)<0.0001) break;
			if(tneu<(min(yl)-0.1)) break;
			if(tneu>(max(yl)+0.1)) break;
			if(wertneu >wertt)
			{
				tmax<- tneu
				wertmax<-wertneu
				tmin <- tmin
				wertmin <-wertmin
			} else 
			{
				tmin <- tneu
				wertmin <- wertneu
				wertmax <-wertmax
				tmax<-tmax
			}
		}
	
		r[l]<-(tmin + (wertt-wertmin)*(tmax - tmin)/(wertmax- wertmin))
	}
	fit<-new("monofit",x = xx, y = r)
	return(fit)
}

setClassUnion("fitprob", c("list","locfit", "locpoly"))

setClass("fitED", representation(alpha="numeric", ED="numeric", call="call", fitold="fitprob"))
setClass("ED.Boot.CI", representation(CI="numeric",R="numeric", call="call"))

setGeneric("ED", function(fitprob, alpha, ...)standardGeneric("ED"))
setGeneric("Boot.CI", function(fitprob, alpha, ...)standardGeneric("Boot.CI"))

setGeneric("aic.ED.locfit", function(fit,...)standardGeneric("aic.ED.locfit"))

setMethod("ED", signature(fitprob="list", alpha="missing"), function(fitprob, alpha,bandwidth,N=101, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)
setMethod("ED", signature(fitprob="locfit", alpha="missing"), function(fitprob, alpha,bandwidth,N=101, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)
setMethod("ED", signature(fitprob="locpoly", alpha="missing"), function(fitprob, alpha,bandwidth,N=101, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)


setMethod("ED", signature(fitprob="list", alpha="numeric"), function(fitprob, alpha, bandwidth, N=101, mono="increasing", type="cont",...){
	call<-match.call(call=sys.call(sys.parent()))
	locfit=locfit.raw(fitprob[[1]], fitprob[[2]], deg=1)
	xx=lfmarg(locfit, m=N)
	xmin=min(xx[[1]])
	xmax=max(xx[[1]])
	yy=predict(locfit, newdata=xx)
	if(missing(bandwidth)) hd=sd(fitprob[[1]])/length(fitprob[[1]])^(1/5)
	else hd=bandwidth
	n1=length(alpha)
	if(type=="cont") {alpha1=alpha}
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	if(mono=="decreasing") alpha1=1-alpha1
	dose=numeric(n1)
	for(l in 1:n1){
	W=numeric(N)
	v=numeric(N)
	for(i in 1:N){
	v[i]=(yy[i]-alpha1[l])/hd
	if(mono=="increasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2-3/4*v[i]+1/4*v[i]^3 
	else if(v[i]>1) W[i]=0
	else W[i]=1
	}
	if(mono=="decreasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2+3/4*v[i]-1/4*v[i]^3 
	else if(v[i]>1) W[i]=1
	else W[i]=0
	}

	}	
	dose[l]=xmin +(xmax-xmin)*sum(W)/(N+1)
	}
	dose
	fited<-new("fitED", alpha=alpha, ED=dose, call=call, fitold=fitprob)
	return(fited)
})

setMethod("ED", signature(fitprob="locfit", alpha="numeric"), function(fitprob, alpha, bandwidth, N=101, mono="increasing",type="cont",...){
	call<-match.call(call=sys.call(sys.parent()))
	xx=lfmarg(fitprob, m=N)
	xmin=min(xx[[1]])
	xmax=max(xx[[1]])
	yy=predict(fitprob, newdata=xx)
	
	if(missing(bandwidth)){
	 data <-{ if (is.null(fit$call$data)){ 
       	sys.frame(sys.parent())}
       else eval(fit$call$data)}
       m<-locfit.matrix(fit, data = data)
	  hd=sd(m$x)/length(m$x)^(1/5)
	}
	else hd=bandwidth
	n1=length(alpha)
	if(type=="cont") alpha1=alpha
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	if(mono=="decreasing") alpha1=1-alpha1
	dose=numeric(n1)
	for(l in 1:n1){
	W=numeric(N)
	v=numeric(N)
	for(i in 1:N){
	v[i]=(yy[i]-alpha1[l])/hd
	if(mono=="increasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2-3/4*v[i]+1/4*v[i]^3 
	else if(v[i]>1) W[i]=0
	else W[i]=1
	}
	if(mono=="decreasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2+3/4*v[i]-1/4*v[i]^3 
	else if(v[i]>1) W[i]=1
	else W[i]=0
	}

	}	
	dose[l]=xmin+ (xmax-xmin)*sum(W)/(N+1)
	}
	dose
fited<-new("fitED", alpha=alpha, ED=dose, call=call,  fitold=fitprob)
return(fited)
})


setMethod("ED", signature(fitprob="locpoly", alpha="numeric"), function(fitprob, alpha, bandwidth, N=101, mono="increasing", type="cont",...){
	call<-match.call(call=sys.call(sys.parent()))
	fitprob@call$gridsize<-N
	pfitnew<-eval(fitprob@call)
	xx=pfitnew@x
	xmin=min(xx)
	xmax=max(xx)
	yy=pfitnew@y
	if(missing(bandwidth)) { 
		x=eval(fitprob@call[[2]])
		 hd=sd(x)/length(x)^(1/5)
		}
	else hd=bandwidth
	n1=length(alpha)
	dose=numeric(n1)
	if(type=="cont") alpha1=alpha
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	if(mono=="decreasing") alpha1=1-alpha1
	for(l in 1:n1){
	W=numeric(N)
	v=numeric(N)
	for(i in 1:N){
	v[i]=(yy[i]-alpha1[l])/hd
	if(mono=="increasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2-3/4*v[i]+1/4*v[i]^3 
	else if(v[i]>1) W[i]=0
	else W[i]=1
	}
	if(mono=="decreasing"){
	if(v[i]<=1 && v[i]>=-1) W[i]=1/2+3/4*v[i]-1/4*v[i]^3 
	else if(v[i]>1) W[i]=1
	else W[i]=0
	}

	}	
	dose[l]=xmin+ (xmax-xmin)*sum(W)/(N+1)
	}
	dose
	fited<-new("fitED", alpha=alpha, ED=dose, call=call,  fitold=fitprob)
return(fited)
})

setMethod("aic.ED.locfit", signature(fit="fitED"), function(fit, pen=2, ...){
	locfitob=fit@fitold
	if(class(locfitob)=="locfit"){
	df=locfitob[['dp']][['df1']]	
	data <-{ if (is.null(locfitob$call$data)){ 
       	sys.frame(sys.parent())}
       else eval(locfitob$call$data)}
       m<-locfit.matrix(locfitob, data = data)
       aic=aic(locfitob, pen=pen)['aic']
       LK=aic(locfitob, pen=pen)['lik']
       sy<-m$y
       sx=as.numeric(m$x)
	   if(is.null(fit@call$bandwidth)) bandwidth=sd(m$x)/length(m$x)^(1/5)
		else bandwidth=fit@call$bandwidth
      if(is.null(fit@call$N)) N=101
		else N=fit@call$N
       if(is.null(fit@call$mono)) mono="increasing"
    else mono=fit@call$mono	
	xx=lfmarg(locfitob, m=N)
	yy=predict(locfitob, newdata=xx)
	fits=list(xx[[1]],yy)
	}
	if(class(locfitob)=="list"){
	if(is.null(fit@call$bandwidth)) bandwidth=sd(locfitob[[1]])/length(locfitob[[1]])^(1/5)
	else bandwidth=fit@call$bandwidth
	if(is.null(fit@call$N)) N=101
	else N=fit@call$N
	if(is.null(fit@call$mono)) mono="increasing"
    else mono=fit@call$mono	
	locfit=locfit.raw(locfitob[[1]], locfitob[[2]], deg=1)
	df=locfit[['dp']][['df1']]
	LK=locfit[['dp']][['lk']]
	aic=-2 * LK + pen * df
	xx=lfmarg(locfit, m=N)
	xmin=min(xx[[1]])
	xmax=max(xx[[1]])
	yy=predict(locfit, newdata=xx)
	sy<-locfitob[[2]]
	sx<-locfitob[[1]]
	fits=list(xx[[1]],yy)
	}
	rz=mono.1d(fits, bandwidth=bandwidth, xx=sx, mono1=mono)
	lk=-sum((sy-rz@y)^2)/2
	AIC=-2*lk+pen*df
	r=c(df,LK, lk, aic, AIC)
	names(r)=c("df", "lk LL", "lk ED", "AIC LL", "AIC ED")
	return(r)
	})


setMethod("Boot.CI", signature(fitprob="list", alpha="missing"), function(fitprob, alpha,bandwidth,level=0.95, N=101, R=100, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)
setMethod("Boot.CI", signature(fitprob="locfit", alpha="missing"), function(fitprob, alpha, bandwidth,level=0.95, N=101, R=100, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)
setMethod("Boot.CI", signature(fitprob="locpoly", alpha="missing"), function(fitprob, alpha, bandwidth,level=0.95, N=101, R=100, mono="increasing",...){
cat("the alpha value has to be specified!\n")
}
)


setMethod("Boot.CI", signature(fitprob="list", alpha="numeric"), function(fitprob, alpha, bandwidth, level=0.95, N=101, R=100, mono="increasing",type="cont",...){
	if(length(alpha)!=1) stop("Please use only one alpha-value per time!\n")
	call<-match.call(call=sys.call(sys.parent()))
	x=fitprob[[1]]
	y=fitprob[[2]]
	if(type=="cont") alpha1=alpha
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	if(missing(bandwidth)) hd=sd(x)/length(x)^(1/5)
	else hd=bandwidth
	n=length(fitprob[[1]])
	res=locfit.raw(x,y, deg=1)
	res1=locfit.raw(x,y,deg=1, deriv=1)
	if(abs(predict(res1, newdata=alpha1))<0.1) warning("WARNING: The bootstrap interval might be inaccurate\n")
	xx=lfmarg(res, m=N)
	yy=predict(res, newdata=xx)
	fit=list(xx[[1]],yy)
	rz=mono.1d(fit, bandwidth=hd, xx=x, mono1=mono)
	epsilon=y-rz@y
	phatinv=ED(res, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
	Tstern1=numeric(R)
	for(k in 1:R){
		epsilonstern=(-1)^rbinom(n,1,0.5)*epsilon
		ystern=rz@y+epsilonstern
		mstern1=locfit.raw(x,ystern, deg=1)
		phatsterninv=ED(mstern1, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
		Tstern1[k]=abs(phatsterninv-phatinv)
		}
		tstern1=quantile(Tstern1, level)
		CI=c(phatinv-tstern1, phatinv+tstern1)
		names(CI)=c(paste(100*(1-level)/2, "%"), paste(100- 100*(1-level)/2, "%"))
		ret<-new("ED.Boot.CI", CI=CI, R=R, call=call)	
	return(ret)
	})


setMethod("Boot.CI", signature(fitprob="locfit", alpha="numeric"), function(fitprob, alpha, bandwidth, level=0.95, N=101, R=100, mono="increasing",type="cont",...){
	if(length(alpha)!=1) stop("Please use only one alpha-value per time!\n")
	call<-match.call(call=sys.call(sys.parent()))
	data <-{ if (is.null(fitprob$call$data)){ 
          	sys.frame(sys.parent())}
       	else eval(fitprob$call$data)}
       	m<-locfit.matrix(fitprob, data = data)
	if(missing(bandwidth)) hd=sd(m$x)/length(m$x)^(1/5)
	else hd=bandwidth	
	xx=lfmarg(fitprob, m=N)
	yy=predict(fitprob, newdata=xx)
	fit=list(xx[[1]],yy)
	n=fitprob$mi["n"]
	if(type=="cont") alpha1=alpha
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	rz=mono.1d(fit, bandwidth=hd, xx=as.numeric(m$x), mono1=mono)
	epsilon=as.numeric(m$y)-rz@y
	phatinv=ED(fitprob, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
	Tstern1=numeric(R)
	for(k in 1:R){
		epsilonstern=(-1)^rbinom(n,1,0.5)*epsilon
		ystern=rz@y+epsilonstern
		mstern1=locfit.raw(as.numeric(m$x),ystern, deg=1)
		phatsterninv=ED(mstern1, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
		Tstern1[k]=abs(phatsterninv-phatinv)
		}
		tstern1=quantile(Tstern1, level)
		CI=c(phatinv-tstern1, phatinv+tstern1)
		names(CI)=c(paste(100*(1-level)/2,"%"), paste(100- 100*(1-level)/2,"%"))
		fitprob$call$deriv<-1
		deriv=eval(fitprob$call)
		if(abs(predict(deriv, newdata=alpha1))<0.1) warning("WARNING: The bootstrap interval might be inaccurate\n")

		ret<-new("ED.Boot.CI", CI=CI,R=R, call=call)	
	return(ret)
	})
	

setMethod("Boot.CI", signature(fitprob="locpoly", alpha="numeric"), function(fitprob, alpha, bandwidth, level=0.95, N=101, R=100, mono="increasing", type="cont",...){
	if(length(alpha)!=1) stop("Please use only one alpha-value per time!\n")
	call<-match.call(call=sys.call(sys.parent()))
	x=eval(fitprob@call[[2]])
	y=eval(fitprob@call[[3]])
	if(type=="cont") alpha1=alpha
	if(type=="prob") {
		if(any(alpha>1)&any(alpha<0)) stop("Please use type=cont\n")
		rr=range(tapply(fitprob[[2]],fitprob[[1]], mean))
		alpha1=rr[1]*alpha+rr[2]*(1-alpha)
		}
	n=length(x)
	fitprob@call$gridsize<-N
	pfitnew<-eval(fitprob@call)
	if(missing(bandwidth)) hd=sd(x)/length(x)^(1/5)
	else hd=bandwidth
	xx=pfitnew@x
	xmin=min(xx)
	xmax=max(xx)
	yy=pfitnew@y
	fit=list(xx,yy)
	rz=mono.1d(fit, bandwidth=hd, xx=x, mono1=mono)
	epsilon=y-rz@y
	phatinv=ED(fitprob, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
	Tstern1=numeric(R)
	for(k in 1:R){
		epsilonstern=(-1)^rbinom(n,1,0.5)*epsilon
		ystern=rz@y+epsilonstern
		fitprob@call$y<-ystern
		mstern1=eval(fitprob@call)
		phatsterninv=ED(mstern1, alpha=alpha1, bandwidth=hd, N=N, mono=mono)@ED
		Tstern1[k]=abs(phatsterninv-phatinv)
		}
		tstern1=quantile(Tstern1, level)
		CI=c(phatinv-tstern1, phatinv+tstern1)
		names(CI)=c(paste(100*(1-level)/2,"%"), paste(100- 100*(1-level)/2,"%"))
		fitprob@call$drv<-1
		fitprob@call$y<-y
		deriv=eval(fitprob@call)
		nn=which.min(abs(deriv@x-alpha1))
		if(abs(deriv@y[nn])<0.1) warning("WARNING: The bootstrap interval might be inaccurate\n")
		ret<-new("ED.Boot.CI", CI=CI,R=R, call=call)	
	return(ret)
	})

setMethod("plot", signature(x="fitED", y="missing"), function(x, xlab="alpha", ylab="effective dose", ...){
	plot(x@alpha,x@ED, type="l", xlab="alpha", ylab="effective dose",...)
	})

setMethod("lines", signature(x="fitED"), function(x, ...){
	lines(x@alpha,x@ED, ...)
	})

setMethod("print", signature(x="fitED"), function(x,...){
	cat("Effective Dose level ED\n")
	cat("CALL:\n")
	print(x@call)
	cat("alpha-values:\n")
	print(x@alpha) 
	cat("effective dose levels:\n")
	print(x@ED)
	})
	
setMethod("print", signature(x="ED.Boot.CI"), function(x,...){
	cat("BOOTSTRAP CONFIDENCE INTERVAL for EFFECTIVE DOSE \n")
	cat("based on", x@R, "bootstrap replications\n")
	cat("CALL:\n")
	print(x@call)
	cat("Interval:\n")
	print(x@CI)
	})
	
setMethod("show", signature(object="fitED"), function(object){
	cat("Effective Dose level ED\n")
	cat("CALL:\n")
	print(object@call)
	cat("alpha-values:\n")
	print(object@alpha) 
	cat("effective dose levels:\n")
	print(object@ED)
	})
	
setMethod("show", signature(object="ED.Boot.CI"), function(object){
	cat("BOOTSTRAP CONFIDENCE INTERVAL for EFFECTIVE DOSE \n")
	cat("based on", object@R, "bootstrap replications\n")
	cat("CALL:\n")
	print(object@call)
	cat("Interval:\n")
	print(object@CI)
	})
	
