
### Below is an R function for computing/plotting the median residual time estimator (using right-censored data)
### along with bootstrap-based pointwise confidence intervals.
###
### Input:   'survobj' = survival object
###          'times' = vectors of times at which the median residual time estimator will be computed
###          'confint' = TRUE if confidence intervals requested, FALSE otherwise
###          'conflevel' = confidence level of confidence intervals
###          'bootruns' = number of replications performed in bootstrap scheme
###          'plot' = TRUE if plot requested, FALSE otherwise
###          'xrange' = vector containing the lower and upper ends of the x-axis
###          'yrange' = vector containing the lower and upper ends of the y-axis
###          'xtitle' = title to be displayed under x-axis
###          'ytitle' = title to be displayed next to y-axis
###
### Output:  list with elements 'times' = vector of times at which the median residual time estimator was calculated;
###                       	    'ci.lower' = vector of lower endpoints for each confidence interval computed;
###                             'estimates' = vector of estimated median residual times;
###                             'ci.upper' = vector of upper endpoints for each confidence interval computed.

library(survival)

getmedianres = function(survobj,times,confint=FALSE,conflevel=0.95,bootruns=2000,
					    plot=FALSE,xrange=NULL,yrange=NULL,xtitle=NULL,ytitle=NULL){

	fit = survfit(survobj~1);
	sumfit = summary(fit);
	t = sumfit$time; s = sumfit$surv;
	km = stepfun(t,c(1,s));

	medrl = function(x){
		p = 1-0.5*km(x);
		return(as.numeric(quantile(fit,probs=p,conf.int=FALSE))-x);
	}
	medrl = Vectorize(medrl);
	medrlvalues = medrl(times);

	if(plot){
		
		if(is.null(xrange)) xrange = c(min(times),max(times));
		curve(medrl(x),from=xrange[1],to=xrange[2],ylim=yrange,lwd=1.5,cex.axis=1,cex.lab=1,xlab=xtitle,ylab=ytitle,n=10000);
	
	}

	if(confint){

		y = survobj[,1]; delta = survobj[,2];
	
		makemedrl.vec = function(dat,v){
			fit_s = survfit(Surv(dat$y,dat$delta)~1);
			sumfit_s = summary(fit_s);
			t_s = sumfit_s$time; s_s = sumfit_s$surv;
			km_s = stepfun(t_s,c(1,s_s));
			medrl_s = function(x){
				p_s = 1-0.5*km_s(x);
				return(as.numeric(quantile(fit_s,probs=p_s,conf.int=FALSE))-x);
			}
			medrl_s = Vectorize(medrl_s);
			return(medrl_s(v));
		}

		boottable = matrix(0,nrow=bootruns,ncol=length(times));

		for(i in 1:bootruns){
			
			ind = sample(1:length(y),replace=T);
			newdat = NULL;
			newdat$y = y[ind]; newdat$delta = delta[ind];
			boottable[i,] = makemedrl.vec(newdat,times);
			
		}
		
		int = matrix(0,ncol=2,nrow=length(times));
		lower = (1-conflevel)/2; higher = lower+conflevel;
		for(i in 1:length(times)){
			int[i,] = quantile(boottable[,i],c(0.025,0.975),na.rm=T)
		}

		if(plot) for(i in 1:length(times)) points(rep(times[i],2),int[i,],type='l',lty=3);

		return(list("times"=times,"ci.lower"=int[,1],"estimates"=medrlvalues,"ci.upper"=int[,2]));

	}
	
	return(list("times"=times,"estimate"=medrlvalues));

}