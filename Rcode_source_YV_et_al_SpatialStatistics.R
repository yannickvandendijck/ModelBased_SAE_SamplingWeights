#-----------------------------------------------#
#								#
# Model-Based Inference for Small Area		#
# Estimation with Sampling Weights			#
#								#
# Y. Vandenijck et al. (2016)				#
# 								#
# Spatial Statistics					#
#								#
# Source R-code to perform analysis			#
#								#
#-----------------------------------------------#


### unweighted data
unwMean <- function(simdata, arronddata){
	y.help = tapply(simdata$y, simdata$arrond, sum)
	n.help = tapply(simdata$n, simdata$arrond, sum)
	x = y.help/n.help

	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2 = simdata2[order(simdata2$arrond),]
	simdata2$wn=1
	y.ext = c(rep(rep(1,length(simdata2$y)),times=simdata2$y) , rep(rep(0,length(simdata2$y)),times=(simdata2$n-simdata2$y)))
	arrond.ext = c(rep(simdata2$arrond,times=simdata2$y) , rep(simdata2$arrond,times=(simdata2$n-simdata2$y)))
	wn.ext = c(rep(simdata2$wn,times=simdata2$y) , rep(simdata2$wn,times=(simdata2$n-simdata2$y)))
	popsize.ext = c(rep(simdata2$popsize,times=simdata2$y) , rep(simdata2$popsize,times=(simdata2$n-simdata2$y)))
	extended.data = data.frame(y=y.ext, arrond=arrond.ext, wn=wn.ext, popsize=popsize.ext)
	my.svydesign = svydesign(id=~1, strata=~arrond, weights=~wn, data=extended.data, fpc=~popsize)
	design.consistent.estimates = svyby(~y, ~arrond, svymean, design=my.svydesign)
	return(list(x, design.consistent.estimates$se))
}


### weighted data
weightedMean <- function(simdata, arronddata){
	y.help = tapply(simdata$y*simdata$wn, simdata$arrond, sum)
	n.help = tapply(simdata$n*simdata$wn, simdata$arrond, sum)
	x = y.help/n.help

	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2 = simdata2[order(simdata2$arrond),]
	y.ext = c(rep(rep(1,length(simdata2$y)),times=simdata2$y) , rep(rep(0,length(simdata2$y)),times=(simdata2$n-simdata2$y)))
	arrond.ext = c(rep(simdata2$arrond,times=simdata2$y) , rep(simdata2$arrond,times=(simdata2$n-simdata2$y)))
	wn.ext = c(rep(simdata2$wn,times=simdata2$y) , rep(simdata2$wn,times=(simdata2$n-simdata2$y)))
	popsize.ext = c(rep(simdata2$popsize,times=simdata2$y) , rep(simdata2$popsize,times=(simdata2$n-simdata2$y)))
	extended.data = data.frame(y=y.ext, arrond=arrond.ext, wn=wn.ext, popsize=popsize.ext)
	my.svydesign = svydesign(id=~1, strata=~arrond, weights=~wn, data=extended.data, fpc=~popsize)
	design.consistent.estimates = svyby(~y, ~arrond, svymean, design=my.svydesign)
	return(list(x, design.consistent.estimates$se))
}


### Unadjusted binomial model
unadjusted.binomial <- function(simdata,arronddata,graph.loc){
	sumy = tapply(simdata$y, simdata$arrond, sum)
	help = data.frame(arrond=as.numeric(names(sumy)), sumy=sumy)
	help2 = merge(arronddata, help, by="arrond", all.x=TRUE)
	help2$region.unstruct = help2$index
	help2$region.struct = help2$index

	formula = sumy ~ f(region.struct,model="besag", graph=graph.loc, param=c(0.5,0.008)) + 
		f(region.unstruct,model="iid", param=c(0.5,0.008)) 
	mod <- inla(formula, family = "binomial", data = help2, Ntrials=samplesize, control.predictor=list(compute=TRUE,link=1))
	p.est <- mod$summary.fitted.values[, "0.5quant"]
	p.ll  <- mod$summary.fitted.values[, "0.025quant"]
	p.ul  <- mod$summary.fitted.values[, "0.975quant"]
	return(list(p.est, p.ll, p.ul))
}


### Logit normal model
logit.normal <- function(simdata,arronddata,graph.loc){

	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2 = simdata2[order(simdata2$arrond),]
	y.ext = c(rep(rep(1,length(simdata2$y)),times=simdata2$y) , rep(rep(0,length(simdata2$y)),times=(simdata2$n-simdata2$y)))
	arrond.ext = c(rep(simdata2$arrond,times=simdata2$y) , rep(simdata2$arrond,times=(simdata2$n-simdata2$y)))
	wn.ext = c(rep(simdata2$wn,times=simdata2$y) , rep(simdata2$wn,times=(simdata2$n-simdata2$y)))
	popsize.ext = c(rep(simdata2$popsize,times=simdata2$y) , rep(simdata2$popsize,times=(simdata2$n-simdata2$y)))
	extended.data = data.frame(y=y.ext, arrond=arrond.ext, wn=wn.ext, popsize=popsize.ext)
	my.svydesign = svydesign(id=~1, strata=~arrond, weights=~wn, data=extended.data, fpc=~popsize)
	design.consistent.estimates = svyby(~y, ~arrond, svymean, design=my.svydesign)

		# empirical logistics transform
	var.est = design.consistent.estimates$se^2
	logit.p = log(design.consistent.estimates$y/(1-design.consistent.estimates$y))
	logit.prec = 1/var.est * (design.consistent.estimates$y*(1-design.consistent.estimates$y))^2

	simdata4 = data.frame(arrond=design.consistent.estimates$arrond, logit.p=logit.p, logit.prec=logit.prec)
	simdata5 = merge(arronddata, simdata4, by="arrond", all.x=TRUE)
	simdata5$region.unstruct = simdata5$index
	simdata5$region.struct = simdata5$index

	formula = logit.p ~ 1 + f(region.unstruct,model="iid",param=c(0.5,0.008)) +
		f(region.struct,model="besag", graph=graph.loc, param=c(0.5,0.008))
	mod <- inla(formula, family = "gaussian", data =simdata5, 
		control.predictor=list(compute=TRUE,link=1),
		control.family=list(hyper=list(prec=list(initial=log(1),fixed=TRUE))),
		scale=logit.prec)

	# backtransform results
	p.est = p.ll = p.ul = NULL
	for (i in 1:43){
		m = inla.tmarginal( function(x) exp(x)/(1 +exp(x)), mod$marginals.fitted.values[[i]])
		p.est[i]= inla.qmarginal(c(0.025,0.50,0.975), m)[2]
		p.ll[i] = inla.qmarginal(c(0.025,0.50,0.975), m)[1]
		p.ul[i] = inla.qmarginal(c(0.025,0.50,0.975), m)[3]
	}
	return(list(p.est, p.ll, p.ul, mod))
}


### arcsin normal model
arcsin.normal <- function(simdata,arronddata,graph.loc){

	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2 = simdata2[order(simdata2$arrond),]
	y.ext = c(rep(rep(1,length(simdata$y)),times=simdata$y) , rep(rep(0,length(simdata$y)),times=(simdata$n-simdata$y)))
	arrond.ext = c(rep(simdata$arrond,times=simdata$y) , rep(simdata$arrond,times=(simdata$n-simdata$y)))
	wn.ext = c(rep(simdata$wn,times=simdata$y) , rep(simdata$wn,times=(simdata$n-simdata$y)))
	popsize.ext = c(rep(simdata2$popsize,times=simdata$y) , rep(simdata2$popsize,times=(simdata$n-simdata$y)))
	extended.data = data.frame(y=y.ext, arrond=arrond.ext, wn=wn.ext, popsize=popsize.ext)
	my.svydesign = svydesign(id=~1, strata=~arrond, weights=~wn, data=extended.data, fpc=~popsize)
	design.consistent.estimates = svyby(~y, ~arrond, svymean, design=my.svydesign)

		# arcsine transformation
	var.est = design.consistent.estimates$se^2
	arcsin.p = asin(sqrt(design.consistent.estimates$y))
	eff.sample.size = design.consistent.estimates$y*(1-design.consistent.estimates$y) / var.est
	arcsin.prec = 4*eff.sample.size

	simdata4 = data.frame(arrond=design.consistent.estimates$arrond, arcsin.p=arcsin.p, arcsin.prec=arcsin.prec)
	simdata5 = merge(arronddata, simdata4, by="arrond", all.x=TRUE)
	simdata5$region.unstruct = simdata5$index
	simdata5$region.struct = simdata5$index

	formula = arcsin.p~ 1 + f(region.unstruct,model="iid",param=c(0.5,0.008)) +
		f(region.struct,model="besag", graph=graph.loc, param=c(0.5,0.008))
	mod <- inla(formula, family = "gaussian", data =simdata5, 
		control.predictor=list(compute=TRUE,link=1),
		control.family=list(hyper=list(prec=list(initial=log(1),fixed=TRUE))),
		scale=arcsin.prec)

	# backtransform results
	p.est = p.ll = p.ul = NULL
	for (i in 1:43){
		m = inla.tmarginal( function(x) sin(x)^2 , mod$marginals.fitted.values[[i]])
		p.est[i]= inla.qmarginal(c(0.025,0.50,0.975), m)[2]
		p.ll[i] = inla.qmarginal(c(0.025,0.50,0.975), m)[1]
		p.ul[i] = inla.qmarginal(c(0.025,0.50,0.975), m)[3]
	}
	return(list(p.est, p.ll, p.ul))
}


### pseudo likelihood
pseudo.likelihood.binomial <- function(simdata, arronddata, graph.loc){
	y.pseudo.lik = NULL
	k=1
	for (i in unique(simdata$arrond)){
		set = simdata[simdata$arrond==i,]
		y.pseudo.lik[k] = sum(set$wn*set$y)
		k=k+1
	}

	simdata2 = arronddata
	simdata2$y.pseudo.lik = NA
	simdata2$y.pseudo.lik[!is.na(simdata2$samplesize)] = y.pseudo.lik
	simdata2$region.unstruct = simdata2$index
	simdata2$region.struct = simdata2$index

	formula = y.pseudo.lik ~ f(region.unstruct,model="iid", param=c(0.5,0.008)) + 
		f(region.struct,model="besag", graph=graph.loc,param=c(0.5,0.008))
	mod <- inla(formula, family = "binomial", data = simdata2, 
		Ntrials=samplesize, control.predictor=list(compute=TRUE,link=1))
	p.est <- mod$summary.fitted.values[, "0.5quant"]
	p.ll  <- mod$summary.fitted.values[, "0.025quant"]
	p.ul  <- mod$summary.fitted.values[, "0.975quant"]
	return(list(p.est, p.ll, p.ul))
}


### Effective sample size method
effective.samplesize.binomial <- function(simdata, arronddata, graph.loc){

	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2 = simdata2[order(simdata2$arrond),]
	y.ext = c(rep(rep(1,length(simdata$y)),times=simdata$y) , rep(rep(0,length(simdata$y)),times=(simdata$n-simdata$y)))
	arrond.ext = c(rep(simdata$arrond,times=simdata$y) , rep(simdata$arrond,times=(simdata$n-simdata$y)))
	wf.ext = c(rep(simdata$wf,times=simdata$y) , rep(simdata$wf,times=(simdata$n-simdata$y)))
	popsize.ext = c(rep(simdata2$popsize,times=simdata$y) , rep(simdata2$popsize,times=(simdata$n-simdata$y)))
	extended.data = data.frame(y=y.ext, arrond=arrond.ext, wf=wf.ext, popsize=popsize.ext)
	my.svydesign = svydesign(id=~1, strata=~arrond, weights=~wf, data=extended.data, fpc=~popsize)
	design.consistent.estimates = svyby(~y, ~arrond, svymean, design=my.svydesign)

		# effective sample size transformation
	var.est = design.consistent.estimates$se^2
	eff.sample.size = design.consistent.estimates$y*(1-design.consistent.estimates$y) / var.est
	y.eff.sample.size = design.consistent.estimates$y*eff.sample.size

	simdata4 = data.frame(arrond=design.consistent.estimates$arrond, y.eff.sample.size=y.eff.sample.size, eff.sample.size=eff.sample.size)
	simdata5 = merge(arronddata, simdata4, by="arrond", all.x=TRUE)
	simdata5$region.unstruct = simdata5$index
	simdata5$region.struct = simdata5$index

	formula = y.eff.sample.size ~ f(region.unstruct,model="iid", param=c(0.5,0.008)) +
		f(region.struct,model="besag", graph=graph.loc,param=c(0.5,0.008))
  	mod = inla(formula, family = "binomial", data = simdata5, Ntrials=eff.sample.size, control.predictor=list(compute=TRUE,link=1))
	p.est <- mod$summary.fitted.values[, "0.5quant"]
	p.ll  <- mod$summary.fitted.values[, "0.025quant"]
	p.ul  <- mod$summary.fitted.values[, "0.975quant"]
	return(list(p.est, p.ll, p.ul))
}


### proposed method model 1
model.based.method1 <- function(simdata, arronddata, graph.loc,B){
		# create data for model predictions
	pred.data1 = merge(simdata, arronddata, by="arrond")
	pred.data1$region.struct = pred.data1$index
	pred.data1$region.unstruct = pred.data1$index
	pred.data2 = pred.data1
	
	# model predictions using inla
	formula = (y ~ f(wn, model="rw1", scale.model=TRUE, hyper = list(prec = list(prior="loggamma",param=c(1,0.01)))) + 
		f(region.unstruct,model="iid", param=c(0.5,0.008)) + 
		f(region.struct,model="besag", graph=graph.loc, param=c(0.5,0.008)) )
	mod.pred <- inla(formula, family = "binomial", data = pred.data2, Ntrials=n, 
		control.predictor=list(compute=TRUE,link=1),
		control.compute=list(config = TRUE), num.threads=1)
	pred.data2$pred = mod.pred$summary.fitted.values[, "0.5quant"]

		# population size predictions
	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2$idx=simdata2$arrond
	simdata2$jdx=1:nrow(simdata)
	simdata2$x = 1

	formula = n ~ -1 +
	     f(idx, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
	     f(jdx, x, model="iid", constr = FALSE, hyper = list(prec = list(initial = log(0.0001),fixed = TRUE)))
	mod.N = inla(formula, family = "poisson", data = simdata2, control.predictor=list(compute=TRUE,link=1), control.compute = list(config=TRUE))

	norm.N = vector(); k=1
	for (i in arronddata$arrond){
		eta = exp(mod.N$summary.random$jdx$mean[simdata2$arrond==i]) * simdata2$wf[simdata2$arrond==i]
		norm.N = c(norm.N, eta/sum(eta)*arronddata$popsize[arronddata$arrond==i])
		k=k+1
	}
	simdata2$norm.N = norm.N

		# obtain model-based prediction
	model.based.pred = vector()
	k=1
	for (i in arronddata$arrond){
		help.y = simdata2$y[simdata2$arrond==i]	
		if ( length(help.y) > 0 ){
			help.N = simdata2$norm.N[simdata2$arrond==i]
			help.n = simdata2$n[simdata2$arrond==i]
			help.pred = pred.data2$pred[pred.data2$arrond==i]
			term2 = sum((help.N-help.n)*help.pred)
			term1 = sum(help.y)
			model.based.pred[k] = sum(term1+term2) / sum(help.N)
		}
		if (length(help.y)==0){
			beta0 = mod.pred$summary.fixed[1,1]
			rw.coef = merge(data.frame(ID=pred.data1$wn), data.frame(mod.pred$summary.random$wn), by="ID", sort=FALSE, all=TRUE)[,2]
			v.k = mod.pred$summary.random$region.struct[,2]
			index = arronddata$index[arronddata$arrond==i]
			pred.off.sample = expit( beta0 + rw.coef + v.k[index] )
			model.based.pred[k] = 1/sum(simdata2$norm.N) * sum(simdata2$norm.N * pred.off.sample)
		}
		k=k+1
	}

		# obtain posterior variance
	model.based.posterior = matrix(NA, B, 43)
	posterior.pred.y = inla.posterior.sample(B, mod.pred)
	posterior.pred.N = inla.posterior.sample(B, mod.N)

	set.seed(123)
	for (b in 1:B){
		posterior.pred = expit( posterior.pred.y[[b]]$latent[1:nrow(simdata)] )
		i1 = nrow(simdata) + length(unique(simdata2$arrond)) + 1
		i2 = 2*nrow(simdata) + length(unique(simdata2$arrond))
		posterior.N = exp( posterior.pred.N[[b]]$latent[ i1:i2 ] ) * simdata2$wf

		help1=tapply(posterior.N, simdata2$arrond, sum)
		help2=tapply(posterior.N, simdata2$arrond, length)
		norm.N = posterior.N / rep(help1,help2) * simdata2$popsize

		k=1
		for (i in arronddata$arrond){
			help.y = simdata2$y[simdata2$arrond==i]	
			term1 = sum(help.y)
			if (length(help.y) > 0){
				help.N = norm.N[simdata2$arrond==i]
				help.n = simdata2$n[simdata2$arrond==i]
				help.pred = posterior.pred[pred.data2$arrond==i]
				term2 = sum((help.N-help.n)*help.pred)
				model.based.posterior[b,k] = sum(term1+term2) / sum(help.N)
			}
			if (length(help.y)==0){
				n.diff.rw = nrow(mod.pred$summary.random$wn)
				n.diff.uk = nrow(mod.pred$summary.random$region.unstruct)
				n.diff.vk = nrow(mod.pred$summary.random$region.struct)
				posterior.rw.coef = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+1):(nrow(pred.data2)+n.diff.rw),]
				mod.pred$summary.random$wn$post.wn = posterior.rw.coef
				posterior.rw.coef.all = merge(data.frame(ID=pred.data1$wn), data.frame(mod.pred$summary.random$wn), by="ID", sort=FALSE, all=TRUE)[,9]
				posterior.v.k = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+n.diff.rw+n.diff.uk+1):(nrow(pred.data2)+n.diff.rw+n.diff.uk+n.diff.vk),]
				posterior.beta0 = posterior.pred.y[[b]]$latent[nrow(pred.data2)+n.diff.rw+n.diff.uk+n.diff.vk+1,]

				index = arronddata$index[arronddata$arrond==i]
				posterior.pred.off.sample = expit( posterior.beta0 + posterior.rw.coef.all + posterior.v.k[index] )
				term2 = sum( norm.N * posterior.pred.off.sample )
				model.based.posterior[b,k] = sum(term2) / sum(norm.N)
			}
			k=k+1
		}
	}

	x =list(est.p = model.based.pred, posterior.p=model.based.posterior, model.fit=mod.pred)
	return(x)
}


### proposed method model 2
model.based.method2 <- function(simdata, arronddata, graph.loc,B){
		# create data for model predictions
	simdata$pi = 1/simdata$wn
	pred.data1 = merge(simdata, arronddata, by="arrond")
	pred.data1$region.struct = pred.data1$index
	pred.data1$region.unstruct = pred.data1$index
	pred.data2 = pred.data1
	
	# model predictions using inla
	formula = (y ~ f(pi, model="rw1", scale.model=TRUE, hyper = list(prec = list(prior="loggamma",param=c(1,0.01)))) + 
		f(region.unstruct,model="iid", param=c(0.5,0.008)) + 
		f(region.struct,model="besag", graph=graph.loc, param=c(0.5,0.008)) )
	mod.pred <- inla(formula, family = "binomial", data = pred.data2, Ntrials=n, 
		control.predictor=list(compute=TRUE,link=1),
		control.compute=list(config = TRUE), num.threads=1)
	pred.data2$pred = mod.pred$summary.fitted.values[, "0.5quant"]

		# population size predictions
	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2$idx=simdata2$arrond
	simdata2$jdx=1:nrow(simdata)
	simdata2$x = 1

	formula = n ~ -1 +
	     f(idx, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
	     f(jdx, x, model="iid", constr = FALSE, hyper = list(prec = list(initial = log(0.0001),fixed = TRUE)))
	mod.N = inla(formula, family = "poisson", data = simdata2, control.predictor=list(compute=TRUE,link=1), control.compute = list(config=TRUE))

	norm.N = vector(); k=1
	for (i in arronddata$arrond){
		eta = exp(mod.N$summary.random$jdx$mean[simdata2$arrond==i]) * simdata2$wf[simdata2$arrond==i]
		norm.N = c(norm.N, eta/sum(eta)*arronddata$popsize[arronddata$arrond==i])
		k=k+1
	}
	simdata2$norm.N = norm.N

		# obtain model-based prediction
	model.based.pred = vector()
	k=1
	for (i in arronddata$arrond){
		help.y = simdata2$y[simdata2$arrond==i]	
		if ( length(help.y) > 0 ){
			help.N = simdata2$norm.N[simdata2$arrond==i]
			help.n = simdata2$n[simdata2$arrond==i]
			help.pred = pred.data2$pred[pred.data2$arrond==i]
			term2 = sum((help.N-help.n)*help.pred)
			term1 = sum(help.y)
			model.based.pred[k] = sum(term1+term2) / sum(help.N)
		}
		if (length(help.y)==0){
			beta0 = mod.pred$summary.fixed[1,1]
			rw.coef = merge(data.frame(ID=pred.data1$pi), data.frame(mod.pred$summary.random$pi), by="ID", sort=FALSE, all=TRUE)[,2]
			v.k = mod.pred$summary.random$region.struct[,2]
			index = arronddata$index[arronddata$arrond==i]
			pred.off.sample = expit( beta0 + rw.coef + v.k[index] )
			model.based.pred[k] = 1/sum(simdata2$norm.N) * sum(simdata2$norm.N * pred.off.sample)
		}
		k=k+1
	}

		# obtain posterior variance
	model.based.posterior = matrix(NA, B, 43)
	posterior.pred.y = inla.posterior.sample(B, mod.pred)
	posterior.pred.N = inla.posterior.sample(B, mod.N)

	set.seed(123)
	for (b in 1:B){
		posterior.pred = expit( posterior.pred.y[[b]]$latent[1:nrow(simdata)] )
		i1 = nrow(simdata) + length(unique(simdata2$arrond)) + 1
		i2 = 2*nrow(simdata) + length(unique(simdata2$arrond))
		posterior.N = exp( posterior.pred.N[[b]]$latent[ i1:i2 ] ) * simdata2$wf

		help1=tapply(posterior.N, simdata2$arrond, sum)
		help2=tapply(posterior.N, simdata2$arrond, length)
		norm.N = posterior.N / rep(help1,help2) * simdata2$popsize

		k=1
		for (i in arronddata$arrond){
			help.y = simdata2$y[simdata2$arrond==i]	
			term1 = sum(help.y)
			if (length(help.y) > 0){
				help.N = norm.N[simdata2$arrond==i]
				help.n = simdata2$n[simdata2$arrond==i]
				help.pred = posterior.pred[pred.data2$arrond==i]
				term2 = sum((help.N-help.n)*help.pred)
				model.based.posterior[b,k] = sum(term1+term2) / sum(help.N)
			}
			if (length(help.y)==0){
				n.diff.rw = nrow(mod.pred$summary.random$pi)
				n.diff.uk = nrow(mod.pred$summary.random$region.unstruct)
				n.diff.vk = nrow(mod.pred$summary.random$region.struct)
				posterior.rw.coef = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+1):(nrow(pred.data2)+n.diff.rw),]
				mod.pred$summary.random$pi$post.pi = posterior.rw.coef
				posterior.rw.coef.all = merge(data.frame(ID=pred.data1$pi), data.frame(mod.pred$summary.random$pi), by="ID", sort=FALSE, all=TRUE)[,9]
				posterior.v.k = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+n.diff.rw+n.diff.uk+1):(nrow(pred.data2)+n.diff.rw+n.diff.uk+n.diff.vk),]
				posterior.beta0 = posterior.pred.y[[b]]$latent[nrow(pred.data2)+n.diff.rw+n.diff.uk+n.diff.vk+1,]

				index = arronddata$index[arronddata$arrond==i]
				posterior.pred.off.sample = expit( posterior.beta0 + posterior.rw.coef.all + posterior.v.k[index] )
				term2 = sum( norm.N * posterior.pred.off.sample )
				model.based.posterior[b,k] = sum(term2) / sum(norm.N)
			}
			k=k+1
		}
	}

	x =list(est.p = model.based.pred, posterior.p=model.based.posterior, model.fit=mod.pred)
	return(x)
}


### proposed method model 3
model.based.method3 <- function(simdata, arronddata, graph.loc,B){
		# create data for model predictions
	pred.data1 = merge(simdata, arronddata, by="arrond")
	pred.data1$region.struct = pred.data1$index
	pred.data1$region.unstruct = pred.data1$index
	pred.data2 = pred.data1

		# create spline basis functions
	a=range(pred.data1$wn)[1];b=range(pred.data1$wn)[2];
	intKnots <- quantile(unique(pred.data1$wn),
                     seq(0,1,length=(17+2))[-c(1,(17+2))])
	intKnots = as.vector(intKnots)
	B.basis = bs(pred.data2$wn, knots=intKnots, degree=2, Boundary.knots=c(a,b), intercept=TRUE)
	P = diff(diag(20), diff = 2)
	K = t(P)%*% P
	eigOmega <- eigen(K)
	eigOmega$values
	indsZ <- 1:18
	UZ <- eigOmega$vectors[,indsZ]
	LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
	ZSpline <- B.basis%*%LZ
	pred.data2$idnum = 1:nrow(pred.data2)

		# model predictions using inla
	formula = (y ~ wn + f(idnum, model="z", Z=ZSpline, initial=5, param=c(1,0.01)) +
		f(region.unstruct,model="iid", param=c(0.5,0.008)) + 
		f(region.struct,model="besag", constr=TRUE, graph=graph.loc, param=c(0.5,0.008)) )
	mod.pred.spline <- inla(formula, family = "binomial", data = pred.data2, Ntrials=n, 
		control.predictor=list(compute=TRUE,link=1),
		control.compute=list(config = TRUE), num.threads=1)
	pred.data2$pred.spline = mod.pred.spline$summary.fitted.values[,"0.5quant"]

		# population size predictions
	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2$idx=simdata2$arrond
	simdata2$jdx=1:nrow(simdata)
	simdata2$x = 1

	formula = n ~ -1 +
	     f(idx, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
	     f(jdx, x, model="iid", constr = FALSE, hyper = list(prec = list(initial = log(0.0001),fixed = TRUE)))
	mod.N = inla(formula, family = "poisson", data = simdata2, control.predictor=list(compute=TRUE,link=1), control.compute = list(config=TRUE))

	norm.N = vector(); k=1
	for (i in arronddata$arrond){
		eta = exp(mod.N$summary.random$jdx$mean[simdata2$arrond==i]) * simdata2$wf[simdata2$arrond==i]
		norm.N = c(norm.N, eta/sum(eta)*arronddata$popsize[arronddata$arrond==i])
		k=k+1
	}
	simdata2$norm.N = norm.N

		# obtain model-based prediction
	model.based.pred = vector()
	k=1
	for (i in arronddata$arrond){
		help.y = simdata2$y[simdata2$arrond==i]	
		if ( length(help.y) > 0 ){
			help.N = simdata2$norm.N[simdata2$arrond==i]
			help.n = simdata2$n[simdata2$arrond==i]
			help.pred = pred.data2$pred.spline[pred.data2$arrond==i]
			term2 = sum((help.N-help.n)*help.pred)
			term1 = sum(help.y)
			model.based.pred[k] = sum(term1+term2) / sum(help.N)
		}
		if (length(help.y)==0){
			beta0 = mod.pred.spline$summary.fixed[1,1]
			beta1 = mod.pred.spline$summary.fixed[2,1]
			ps.coef = mod.pred.spline$summary.random$idnum[1:nrow(simdata),2]
			v.k = mod.pred.spline$summary.random$region.struct[,2]
			index = arronddata$index[arronddata$arrond==i]
			pred.off.sample = expit( beta0 + beta1*simdata$wn + ps.coef + v.k[index] )
			model.based.pred[k] = 1/sum(simdata2$norm.N) * sum(simdata2$norm.N * pred.off.sample)
		}
		k=k+1
	}

		# obtain posterior variance
	model.based.posterior = matrix(NA, B, 43)
	posterior.pred.y = inla.posterior.sample(B, mod.pred.spline)
	posterior.pred.N = inla.posterior.sample(B, mod.N)

	set.seed(123)
	for (b in 1:B){
		posterior.pred = expit( posterior.pred.y[[b]]$latent[1:nrow(simdata)] )
		i1 = nrow(simdata) + length(unique(simdata2$arrond)) + 1
		i2 = 2*nrow(simdata) + length(unique(simdata2$arrond))
		posterior.N = exp( posterior.pred.N[[b]]$latent[ i1:i2 ] ) * simdata2$wf

		help1=tapply(posterior.N, simdata2$arrond, sum)
		help2=tapply(posterior.N, simdata2$arrond, length)
		norm.N = posterior.N / rep(help1,help2) * simdata2$popsize

		k=1
		for (i in arronddata$arrond){
			help.y = simdata2$y[simdata2$arrond==i]	
			term1 = sum(help.y)
			if (length(help.y) > 0){
				help.N = norm.N[simdata2$arrond==i]
				help.n = simdata2$n[simdata2$arrond==i]
				help.pred = posterior.pred[pred.data2$arrond==i]
				term2 = sum((help.N-help.n)*help.pred)
				model.based.posterior[b,k] = sum(term1+term2) / sum(help.N)
			}
			if (length(help.y)==0){
				n.diff.uk = nrow(mod.pred.spline$summary.random$region.unstruct)
				n.diff.vk = nrow(mod.pred.spline$summary.random$region.struct)
				posterior.ps.coef = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+1):(nrow(pred.data2)+nrow(simdata)),]
				posterior.v.k = posterior.pred.y[[b]]$latent[(2*nrow(pred.data2)+18+n.diff.uk+1):(2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk),]
				posterior.beta0 = posterior.pred.y[[b]]$latent[2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk+1,]
				posterior.beta1 = posterior.pred.y[[b]]$latent[2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk+2,]

				index = arronddata$index[arronddata$arrond==i]
				posterior.pred.off.sample = expit( posterior.beta0 + posterior.beta1*simdata$wn + posterior.ps.coef + posterior.v.k[index] )
				term2 = sum( norm.N * posterior.pred.off.sample )
				model.based.posterior[b,k] = sum(term2) / sum(norm.N)
			}
			k=k+1
		}
	}

	x =list(est.p = model.based.pred, posterior.p=model.based.posterior, model.fit=mod.pred.spline)
	return(x)
}


### proposed method model 4
model.based.method4 <- function(simdata, arronddata, graph.loc,B){
		# create data for model predictions
	simdata$pi = 1/simdata$wn
	pred.data1 = merge(simdata, arronddata, by="arrond")
	pred.data1$region.struct = pred.data1$index
	pred.data1$region.unstruct = pred.data1$index
	pred.data2 = pred.data1
	
		# create spline basis functions
	a=range(pred.data1$pi)[1];b=range(pred.data1$pi)[2];
	intKnots <- quantile(unique(pred.data1$pi),
                     seq(0,1,length=(17+2))[-c(1,(17+2))])
	intKnots = as.vector(intKnots)
	B.basis = bs(pred.data2$pi, knots=intKnots, degree=2, Boundary.knots=c(a,b), intercept=TRUE)
	P = diff(diag(20), diff = 2)
	K = t(P)%*% P
	eigOmega <- eigen(K)
	eigOmega$values
	indsZ <- 1:18
	UZ <- eigOmega$vectors[,indsZ]
	LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))     
	ZSpline <- B.basis%*%LZ
	pred.data2$idnum = 1:nrow(pred.data2)

		# model predictions using inla
	formula = (y ~ pi + f(idnum, model="z", Z=ZSpline, initial=5, param=c(1,0.01)) +
		f(region.unstruct,model="iid", param=c(0.5,0.008)) + 
		f(region.struct,model="besag", constr=TRUE, graph=graph.loc, param=c(0.5,0.008)) )
	mod.pred.spline <- inla(formula, family = "binomial", data = pred.data2, Ntrials=n, 
		control.predictor=list(compute=TRUE,link=1),
		control.compute=list(config = TRUE), num.threads=1)
	pred.data2$pred.spline = mod.pred.spline$summary.fitted.values[, "0.5quant"]

		# population size predictions
	simdata2 = merge(simdata, arronddata, by="arrond")
	simdata2$idx=simdata2$arrond
	simdata2$jdx=1:nrow(simdata)
	simdata2$x = 1

	formula = n ~ -1 +
	     f(idx, model="iid", hyper = list(prec = list(initial = log(0.0001),fixed = TRUE))) +
	     f(jdx, x, model="iid", constr = FALSE, hyper = list(prec = list(initial = log(0.0001),fixed = TRUE)))
	mod.N = inla(formula, family = "poisson", data = simdata2, control.predictor=list(compute=TRUE,link=1), control.compute = list(config=TRUE))

	norm.N = vector(); k=1
	for (i in arronddata$arrond){
		eta = exp(mod.N$summary.random$jdx$mean[simdata2$arrond==i]) * simdata2$wf[simdata2$arrond==i]
		norm.N = c(norm.N, eta/sum(eta)*arronddata$popsize[arronddata$arrond==i])
		k=k+1
	}
	simdata2$norm.N = norm.N

		# obtain model-based prediction
	model.based.pred = vector()
	k=1
	for (i in arronddata$arrond){
		help.y = simdata2$y[simdata2$arrond==i]	
		if ( length(help.y) > 0 ){
			help.N = simdata2$norm.N[simdata2$arrond==i]
			help.n = simdata2$n[simdata2$arrond==i]
			help.pred = pred.data2$pred.spline[pred.data2$arrond==i]
			term2 = sum((help.N-help.n)*help.pred)
			term1 = sum(help.y)
			model.based.pred[k] = sum(term1+term2) / sum(help.N)
		}
		if (length(help.y)==0){
			beta0 = mod.pred.spline$summary.fixed[1,1]
			beta1 = mod.pred.spline$summary.fixed[2,1]
			ps.coef = mod.pred.spline$summary.random$idnum[1:nrow(simdata),2]
			v.k = mod.pred.spline$summary.random$region.struct[,2]
			index = arronddata$index[arronddata$arrond==i]
			pred.off.sample = expit( beta0 + beta1*simdata$pi + ps.coef + v.k[index] )
			model.based.pred[k] = 1/sum(simdata2$norm.N) * sum(simdata2$norm.N * pred.off.sample)
		}
		k=k+1
	}

		# obtain posterior variance
	model.based.posterior = matrix(NA, B, 43)
	posterior.pred.y = inla.posterior.sample(B, mod.pred.spline)
	posterior.pred.N = inla.posterior.sample(B, mod.N)

	set.seed(123)
	for (b in 1:B){
		posterior.pred = expit( posterior.pred.y[[b]]$latent[1:nrow(simdata)] )
		i1 = nrow(simdata) + length(unique(simdata2$arrond)) + 1
		i2 = 2*nrow(simdata) + length(unique(simdata2$arrond))
		posterior.N = exp( posterior.pred.N[[b]]$latent[ i1:i2 ] ) * simdata2$wf

		help1=tapply(posterior.N, simdata2$arrond, sum)
		help2=tapply(posterior.N, simdata2$arrond, length)
		norm.N = posterior.N / rep(help1,help2) * simdata2$popsize

		k=1
		for (i in arronddata$arrond){
			help.y = simdata2$y[simdata2$arrond==i]	
			term1 = sum(help.y)
			if (length(help.y) > 0){
				help.N = norm.N[simdata2$arrond==i]
				help.n = simdata2$n[simdata2$arrond==i]
				help.pred = posterior.pred[pred.data2$arrond==i]
				term2 = sum((help.N-help.n)*help.pred)
				model.based.posterior[b,k] = sum(term1+term2) / sum(help.N)
			}
			if (length(help.y)==0){
				n.diff.uk = nrow(mod.pred.spline$summary.random$region.unstruct)
				n.diff.vk = nrow(mod.pred.spline$summary.random$region.struct)
				posterior.ps.coef = posterior.pred.y[[b]]$latent[(nrow(pred.data2)+1):(nrow(pred.data2)+nrow(simdata)),]
				posterior.v.k = posterior.pred.y[[b]]$latent[(2*nrow(pred.data2)+18+n.diff.uk+1):(2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk),]
				posterior.beta0 = posterior.pred.y[[b]]$latent[2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk+1,]
				posterior.beta1 = posterior.pred.y[[b]]$latent[2*nrow(pred.data2)+18+n.diff.uk+n.diff.vk+2,]

				index = arronddata$index[arronddata$arrond==i]
				posterior.pred.off.sample = expit( posterior.beta0 + posterior.beta1*simdata$pi + posterior.ps.coef + posterior.v.k[index] )
				term2 = sum( norm.N * posterior.pred.off.sample )
				model.based.posterior[b,k] = sum(term2) / sum(norm.N)
			}
			k=k+1
		}
	}
	x =list(est.p = model.based.pred, posterior.p=model.based.posterior, model.fit=mod.pred.spline)
	return(x)
}




