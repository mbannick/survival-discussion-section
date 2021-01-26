
### This updated version was created on 01/24/2018 by Marco Carone
## JW added a couple of other useful functions: 
##      get.nelson.aalen
##      get.weighted.logrank.estimates
##      get.parametric.estimates

fitparametric = function(survobj,dist,feature=NULL,pi=0.5,t=NULL,t0=NULL,g=NULL,conflevel=0.95){
    	
	f = flexsurvreg(survobj~1, dist=dist)
    alpha = 1-conflevel; q = qnorm(1-alpha/2);
    loglik = round(f$loglik,2); n.parameters = f$npars; n.total = f$N; n.events = f$events;
  
  	modelout = list("fit"=f,"loglik"=loglik,"n.parameters"=n.parameters,"n.total"=n.total,"n.events"=n.events,"alpha"=alpha)
    
    if(dist=="exp"){
    	
    		lambda.est = f$res[1,1]; lambda.se = f$res[1,4];
    		lambda.lower = lambda.est-q*lambda.se;
    		lambda.upper = lambda.est+q*lambda.se;
    		
    		coeff = matrix(0,ncol=4,nrow=1)
    		coeff[1,1] = lambda.est; coeff[1,2] = lambda.lower; coeff[1,3] = lambda.upper; coeff[1,4] = lambda.se;
    		colnames(coeff) = c("estimate","ci.lower","ci.upper","se")
    		parnames = c("lambda");
    		rownames(coeff) = parnames;
    		
    		varcov = matrix(lambda.se,nrow=1,ncol=1);
    		rownames(varcov) = parnames; colnames(varcov) = parnames;
    		
    		familyfit = "### SUMMARY OF FITTED EXPONENTIAL MODEL ###"

    		if(!is.null(feature)||!is.null(g)){
    		
    			if(!is.null(feature)){
	    			if(feature=="mean") g=function(u) 1/u[1]; 
    				if(feature=="quantile") g=function(u) (1/u[1])*(-log(1-pi));
    				if(feature=="survival"&&!is.null(t)) g=function(u) exp(-u[1]*t);
				if(feature=="condsurvival"&&!is.null(t)&&!is.null(t0)) g=function(u) exp(-u[1]*t+u[1]*t0);
			}
			
			g.est = g(c(lambda.est));
			der = grad(g,c(lambda.est));
			g.se = abs(der)*lambda.se;
     		g.lower = g.est-q*g.se;
     		g.upper = g.est+q*g.se;
     	
     	}
     	
     	coeffout = list("coeff"=coeff,"varcov"=varcov);  		
    		
    }
    
    if(dist=="weibull"){
    	
    		shape.est = f$res[1,1]; shape.se = f$res[1,4];
		scale.est = f$res[2,1]; scale.se = f$res[2,4];

		sigma = f$cov;
		transf.matrix = matrix(c(0,-1/scale.est,shape.est,0),byrow=TRUE,nrow=2);
		newsigma = transf.matrix%*%sigma%*%t(transf.matrix)

		lambda.est = 1/scale.est; lambda.se = sqrt(newsigma[1,1]);
    		lambda.lower = lambda.est-q*lambda.se;
    		lambda.upper = lambda.est+q*lambda.se;
    	
    		p.est = shape.est; p.se = sqrt(newsigma[2,2]);
    		p.lower = p.est-q*p.se;
    		p.upper = p.est+q*p.se;
    		
    		coeff = matrix(0,ncol=4,nrow=2)
    		coeff[1,1] = lambda.est; coeff[1,2] = lambda.lower; coeff[1,3] = lambda.upper; coeff[1,4] = lambda.se;
    		coeff[2,1] = p.est; coeff[2,2] = p.lower; coeff[2,3] = p.upper; coeff[2,4] = p.se;
    		colnames(coeff) = c("estimate","ci.lower","ci.upper","se")
    		parnames = c("lambda","p");
   		rownames(coeff) = parnames;	
   		
    		varcov = newsigma;
    		rownames(varcov) = parnames; colnames(varcov) = parnames;
    		
    		familyfit = "### SUMMARY OF FITTED WEIBULL MODEL ###"    		
    		
    		if(!is.null(feature)||!is.null(g)){ ### first argument = lambda; second argument = p

    			if(!is.null(feature)){    		
	    			if(feature=="mean") g=function(u) gamma(1+1/u[2])/u[1]; 
    				if(feature=="quantile") g=function(u) (1/u[1])*(-log(1-pi))^(1/u[2]);
    				if(feature=="survival"&&!is.null(t)) g=function(u) exp(-(u[1]*t)^u[2]);
				if(feature=="condsurvival"&&!is.null(t)&&!is.null(t0)) g=function(u) exp(-(u[1]*t)^u[2]+(u[1]*t0)^u[2]);
			}
			
			g.est = g(c(lambda.est,p.est));
			der = t(as.matrix(grad(g,c(lambda.est,p.est))));
			g.se = sqrt(as.numeric(der%*%newsigma%*%t(der)));
     		g.lower = g.est-q*g.se;
     		g.upper = g.est+q*g.se;
     	
     	}
     	
		coeffout = list("coeff"=coeff,"varcov"=varcov);   		
			
    	}
    	
    	if(dist=="gamma"){
 
    		alpha.est = f$res[1,1]; alpha.se = f$res[1,4];
		lambda.est = f$res[2,1]; lambda.se = f$res[2,4];

		alpha.lower = alpha.est-q*alpha.se;
		alpha.upper = alpha.est+q*alpha.se;
		
		lambda.lower = lambda.est-q*lambda.se;
		lambda.upper = lambda.est+q*lambda.se;
		
		sigma = f$cov;
		transf.matrix = matrix(c(0,lambda.est,alpha.est,0),byrow=TRUE,nrow=2);
		newsigma = transf.matrix%*%sigma%*%t(transf.matrix);
		
		coeff = matrix(0,ncol=4,nrow=2)
    		coeff[1,1] = lambda.est; coeff[1,2] = lambda.lower; coeff[1,3] = lambda.upper; coeff[1,4] = lambda.se;
    		coeff[2,1] = alpha.est; coeff[2,2] = alpha.lower; coeff[2,3] = alpha.upper; coeff[2,4] = alpha.se;
    		colnames(coeff) = c("estimate","ci.lower","ci.upper","se")		
    		parnames = c("lambda","alpha");
   		rownames(coeff) = parnames;	
    		
    		varcov = newsigma;
    		rownames(varcov) = parnames; colnames(varcov) = parnames;
    	
    		familyfit = "### SUMMARY OF FITTED GAMMA MODEL ###"    		

    		if(!is.null(feature)||!is.null(g)){ ### first argument = lambda; second argument = alpha
    		
    			if(!is.null(feature)){
	    			if(feature=="mean") g=function(u) u[2]/u[1]; 
    				if(feature=="quantile") g=function(u) qgamma(pi,rate=u[1],shape=u[2]);
    				if(feature=="survival"&&!is.null(t)) g=function(u) pgamma(t,rate=u[1],shape=u[2],lower.tail=FALSE);
				if(feature=="condsurvival"&&!is.null(t)&&!is.null(t0)){
					g=function(u) pgamma(t,rate=u[1],shape=u[2],lower.tail=FALSE)/
								  pgamma(t0,rate=u[1],shape=u[2],lower.tail=FALSE);
				}
			}
			
			g.est = g(c(lambda.est,alpha.est));
			der = t(as.matrix(grad(g,c(lambda.est,alpha.est))));
			g.se = sqrt(as.numeric(der%*%newsigma%*%t(der)));
     		g.lower = g.est-q*g.se;
     		g.upper = g.est+q*g.se;
     	
     	}
     	
		coeffout = list("coeff"=coeff,"varcov"=varcov);
		    		
    	}
    	
    	if(dist=="gengamma"){
    		
    		mu.est = f$res[1,1]; mu.se = f$res[1,4];
    		sigma.est = f$res[2,1]; sigma.se = f$res[2,4];
    		Q.est = f$res[3,1]; Q.se = f$res[3,4];
    		
    		mu.lower = mu.est-q*mu.se; mu.upper = mu.est+q*mu.se;
    		sigma.lower = sigma.est-q*sigma.se; sigma.upper = sigma.est+q*sigma.se;
    		Q.lower = Q.est-q*Q.se; Q.upper = Q.est+q*Q.se;
    		    		
    		sigma = f$cov;
		transf.matrix = matrix(c(1,0,0,0,sigma.est,0,0,0,1),byrow=TRUE,nrow=3);
		newsigma = transf.matrix%*%sigma%*%t(transf.matrix);
		
		
		coeff = matrix(0,ncol=4,nrow=3)
    		coeff[1,1] = mu.est; coeff[1,2] = mu.lower; coeff[1,3] = mu.upper; coeff[1,4] = mu.se;
    		coeff[2,1] = sigma.est; coeff[2,2] = sigma.lower; coeff[2,3] = sigma.upper; coeff[2,4] = sigma.se;
    		coeff[3,1] = Q.est; coeff[3,2] = Q.lower; coeff[3,3] = Q.upper; coeff[3,4] = Q.se;
    		colnames(coeff) = c("estimate","ci.lower","ci.upper","se")		
    		parnames = c("mu","sigma","Q");
   		rownames(coeff) = parnames;	

	    	varcov = newsigma;
    		rownames(varcov) = parnames; colnames(varcov) = parnames;
    	
    		familyfit = "### SUMMARY OF FITTED GENERALIZED GAMMA MODEL ###"    		

    		if(!is.null(feature)||!is.null(g)){ ### first argument = lambda; second argument = alpha
    		
    			if(!is.null(feature)){
	    			if(feature=="mean"){
	    				upperlim = qgengamma(1-1e-10,mu=mu.est,sigma=sigma.est,Q=Q.est);
	    				g=function(u){
	    					integrate(function(w) pgengamma(w,mu=u[1],sigma=u[2],Q=u[3],lower.tail=FALSE),
	    						  	  lower=0,upper=upperlim)$value;
	    				}
	    			}
    				if(feature=="quantile") g=function(u) qgengamma(pi,mu=u[1],sigma=u[2],Q=u[3]);
    				if(feature=="survival"&&!is.null(t)){
    					g=function(u) pgengamma(t,mu=u[1],sigma=u[2],Q=u[3],lower.tail=FALSE);
    				}
				if(feature=="condsurvival"&&!is.null(t)&&!is.null(t0)){
					g=function(u) pgengamma(t,mu=u[1],sigma=u[2],Q=u[3],lower.tail=FALSE)/
								  pgengamma(t0,mu=u[1],sigma=u[2],Q=u[3],lower.tail=FALSE);
				}
			}
			
			g.est = g(c(mu.est,sigma.est,Q.est));
			der = t(as.matrix(grad(g,c(mu.est,sigma.est,Q.est))));
			g.se = sqrt(as.numeric(der%*%newsigma%*%t(der)));
     		g.lower = g.est-q*g.se;
     		g.upper = g.est+q*g.se;
     	
     	}
     	
		coeffout = list("coeff"=coeff,"varcov"=varcov);
		     	
    	}

	print(noquote(familyfit));    	
	print(noquote(""));        
	print(noquote("MODEL FIT SUMMARIES"));
	print(noquote(""));
	print(noquote(paste("Total number of observations:",n.total," ")));
	print(noquote(paste("Number of events observed:",n.events," ")));
	print(noquote(paste("Number of model parameters:",n.parameters," ")));
	print(noquote(paste("Maximized loglikelihood value:",loglik," ")));
	print(noquote(""));
	print(noquote("INFERENCE ON MODEL COEFFICIENTS"));
	print(noquote(""));
	print(round(coeff,5));
	
	gout = NULL;
	
	if(!is.null(feature)||!is.null(g)){

		if(!is.null(feature)){
			if(feature=="mean") summarydesc = "mean";
    		if(feature=="quantile") summarydesc = paste(pi,"-quantile",sep="");
    		if(feature=="survival"&&!is.null(t)) summarydesc = paste("probability of survival beyond t=",t,sep="");
			if(feature=="condsurvival"&&!is.null(t)&&!is.null(t0)){
				summarydesc = paste("probability of survival beyond t=",t," given survival to t=",t0,sep="");
			}
		} else { summarydesc = "user-defined summary"; }
		
		summaryval = matrix(0,ncol=4,nrow=1);
		summaryval[1,1] = g.est; summaryval[1,2] = g.lower; summaryval[1,3] = g.upper; summaryval[1,4] = g.se;
		colnames(summaryval) = c("estimate","ci.lower","ci.upper","se")		
   		rownames(summaryval) = c("summary");	

		print(noquote(""));   		
		print(noquote("INFERENCE ON SUMMARY"));
		print(noquote(""));
		print(noquote(paste("Summary = ",summarydesc,sep="")));
		print(noquote(""));
		print(round(summaryval,5));
		
		gout = list("feature"=summaryval);
		
	}
	
	invisible(c(modelout,coeffout,gout))
    
}

# Helper functions
## Function to return Nelson-Aalen estimator and standard error
get.nelson.aalen = function(survfit.object, all.times, times.to.show=NULL, alpha=0.05, CI.type=c("loglog, standard")) {
  if (class(survfit.object) != "survfit") {
    stop("need to provide a survfit object")
  }
  CI.type = match.arg(CI.type)
  # want to evaluate at all integer times from t=min(survobject.time) to t=max(survobject.time)
  n.times = length(all.times)
  idx.nonmissing.time = which(all.times %in% survfit.object$time)
  results = matrix(NA, nrow=n.times, ncol=5)
  colnames(results) = c("Time", "Estimate", "SE", "CI_LB", "CI_UB")
  results[, "Time"] = all.times
  results[idx.nonmissing.time, "Estimate"] = cumsum(survfit.object$n.event/survfit.object$n.risk)
  results[idx.nonmissing.time, "SE"] = sqrt(cumsum(survfit.object$n.event/(survfit.object$n.risk * (survfit.object$n.risk - survfit.object$n.event))))
  results[, "Estimate"] = na.locf(results[, "Estimate"])
  results[, "SE"] = na.locf(results[, "SE"])
  if (CI.type == "standard") {
    results[, "CI_LB"] = results[, "Estimate"] - qnorm(p=1-alpha/2) * results[, "SE"]
    results[, "CI_UB"] = results[, "Estimate"] + qnorm(p=1-alpha/2) * results[, "SE"]
  } else {
    results[, "CI_LB"] = results[, "Estimate"] * exp(- qnorm(p=1-alpha/2) * results[, "SE"] / results[, "Estimate"])
    results[, "CI_UB"] = results[, "Estimate"] * exp(qnorm(p=1-alpha/2) * results[, "SE"] / results[, "Estimate"])
  }
  
  if (is.null(times.to.show) == TRUE) {
    return(as.data.frame(results))
  } else {
    #idx.to.show = which(results[, "Time"] %in% times.to.show)
    times.to.show = times.to.show[times.to.show != 0]
    idx.to.show = sapply(times.to.show, function(x) findInterval(x, results[, "Time"]))
    # remove a leading 0 if it's included
    results.to.show = as.data.frame(results)[idx.to.show, ]
    results.to.show[, "Time"] = times.to.show
    return(results.to.show)
  }
}



get.parametric.estimates = function(Surv.object, dist=c("exponential", "weibull", "gengamma.orig")) {
  dist = match.arg(dist)
  survreg.result = flexsurvreg(Surv.object ~ 1, dist=dist)
  if (dist == "exponential") {
    mat.results = matrix(NA, nrow=1, ncol=6)
    colnames(mat.results) = c("Estimate", "SE", "z", "pval", "CI_LB", "CI_UB")
    rownames(mat.results) = "lambda"
    mat.results["lambda", c("Estimate", "SE", "CI_LB", "CI_UB")] = survreg.result$res["rate", c("est", "se", "L95%", "U95%")]
    mat.results["lambda", "z"] = mat.results["lambda", "Estimate"] / mat.results["lambda", "SE"]
    mat.results[, "pval"] = 2 * pnorm(-abs(mat.results["lambda", "z"]))
  } else if (dist == "weibull") {
    mat.results = matrix(NA, nrow=2, ncol=6)
    colnames(mat.results) = c("Estimate", "SE", "z", "pval", "CI_LB", "CI_UB")
    rownames(mat.results) = c("lambda", "p")
    mat.results["lambda", "Estimate"] = 1 / survreg.result$res["scale", "est"]
    mat.results["lambda", "SE"] = deltamethod(g=~1/x1, mean=survreg.result$res["scale", "est"], cov=survreg.result$res["scale", "se"]^2, ses=TRUE)
    mat.results["lambda", c("CI_LB", "CI_UB")] = c(mat.results["lambda", "Estimate"] - 1.96 * mat.results["lambda", "SE"],
                                                   mat.results["lambda", "Estimate"] + 1.96 * mat.results["lambda", "SE"])
    mat.results["p", c("Estimate", "SE", "CI_LB", "CI_UB")] = survreg.result$res["shape", c("est", "se", "L95%", "U95%")]
    mat.results[, "z"] = mat.results[, "Estimate"] / mat.results[, "SE"]
    mat.results[, "pval"] = 2 * pnorm(-abs(mat.results[, "z"]))
  } else if (dist == "gengamma.orig") {
    mat.results = matrix(NA, nrow=3, ncol=6)
    colnames(mat.results) = c("Estimate", "SE", "z", "pval", "CI_LB", "CI_UB")
    rownames(mat.results) = c("lambda", "p", "alpha")
    mat.results["lambda", "Estimate"] = 1 / survreg.result$res["scale", "est"]
    mat.results["lambda", "SE"] = deltamethod(g=~1/x1, mean=survreg.result$res["scale", "est"], cov=survreg.result$res["scale", "se"]^2, ses=TRUE)
    mat.results["lambda", c("CI_LB", "CI_UB")] = c(mat.results["lambda", "Estimate"] - 1.96 * mat.results["lambda", "SE"],
                                                   mat.results["lambda", "Estimate"] + 1.96 * mat.results["lambda", "SE"])
    mat.results["p", c("Estimate", "SE", "CI_LB", "CI_UB")] = survreg.result$res["shape", c("est", "se", "L95%", "U95%")]
    mat.results["alpha", c("Estimate", "SE", "CI_LB", "CI_UB")] = survreg.result$res["k", c("est", "se", "L95%", "U95%")]
    mat.results[, "z"] = mat.results[, "Estimate"] / mat.results[, "SE"]
    mat.results[, "pval"] = 2 * pnorm(-abs(mat.results[, "z"]))
  }
  return(mat.results)
}

get.weighted.logrank.estimates = function(Surv.object, group.variable) {
  group.variable.factor = as.factor(group.variable)
  df = length(levels(group.variable.factor)) - 1
  ## get test statistics and p-values
  test.statistic.tware = statistic(logrank_test(Surv.object ~ group.variable.factor, type="Tarone-Ware")) ^ 2
  pvalue.tware = 1 - pchisq(test.statistic.tware, df=df)
  test.statistic.breslow = statistic(logrank_test(Surv.object ~ group.variable.factor, type="Gehan-Breslow")) ^ 2
  pvalue.breslow = 1 - pchisq(test.statistic.breslow, df=df)
  ## store test statistics and p-values
  mat.weighted.logrank = matrix(NA, nrow=2, ncol=2)
  colnames(mat.weighted.logrank) = c("test statistic", "pvalue")
  rownames(mat.weighted.logrank) = c("tware", "breslow")
  mat.weighted.logrank["tware", ] = c(test.statistic.tware, pvalue.tware)
  mat.weighted.logrank["breslow", ] = c(test.statistic.breslow, pvalue.breslow)
  return(mat.weighted.logrank)
}




