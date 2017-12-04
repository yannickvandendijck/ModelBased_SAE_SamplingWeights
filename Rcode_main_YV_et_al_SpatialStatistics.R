#-----------------------------------------------#
#								#
# Model-Based Inference for Small Area		#
# Estimation with Sampling Weights			#
#								#
# Y. Vandenijck et al. (2016)				#
# 								#
# Spatial Statistics					#
#								#
# Main R-code to perform analysis			#
#								#
#-----------------------------------------------#


	### Libraries
	#------------
library(INLA)
library(survey)
library(VGAM)
library(dummies)
library(rje)
library(TeachingSampling)
library(maptools)
library(maps)
library(RColorBrewer)
library(shapefiles)
library(classInt)


	### Import Data
	#--------------
survey.data = read.table(".../survey_data_YV_et_al_SpatialStatistics.txt",header=TRUE)
arrond.data = read.table(".../arrond_data_YV_et_al_SpatialStatistics.txt",header=TRUE)
str(survey.data)
str(arrond.data)


	### Some Summaries and Plots
	#---------------------------
sum(survey.data$wf * survey.data$n)
sum(arrond.data$popsize)

sum(survey.data$wn * survey.data$n)

plot.hist.data = rep(survey.data$wn, survey.data$n)
hist(plot.hist.data, nclass=100, xlim=c(range(plot.hist.data)[1],range(plot.hist.data)[2]),
	xlab="Normalized weights", ylab="Frequency", cex.lab=1.6, cex.axis=1.6,
	main="Histogram of normalized weights for scenario (S3)")
range(plot.hist.data)

tapply(survey.data$n, survey.data$arrond, sum)
tapply(survey.data$y, survey.data$arrond, sum)


	### Obtain Estimated Prevalences
	#-------------------------------
source(".../Rcode_source_YV_et_al_SpatialStatistics_run.R")
graph.loc = ".../graphloc.graph"
results = arrond.data

 # Unweighted mean (UM)
unw.mean = unwMean(survey.data, arrond.data)
results$est.mean1[!is.na(arrond.data$samplesize)] = unw.mean[[1]]
results$se.mean1[!is.na(arrond.data$samplesize)] = unw.mean[[2]]

 # Weighted mean (HT)
ht.mean = weightedMean(survey.data, arrond.data)
results$est.mean2[!is.na(arrond.data$samplesize)] = ht.mean[[1]]
results$se.mean2[!is.na(arrond.data$samplesize)] = ht.mean[[2]]

 # Naive binomial (NB)
unadj.mean = unadjusted.binomial(survey.data, arrond.data, graph.loc)
results$est.mean3 = unadj.mean[[1]]
results$ll.mean3 = unadj.mean[[2]]
results$ul.mean3 = unadj.mean[[3]]

 # Logit normal (LN)
logit.mean = logit.normal(survey.data, arrond.data, graph.loc)
results$est.mean4 = logit.mean[[1]]
results$ll.mean4 = logit.mean[[2]]
results$ul.mean4 = logit.mean[[3]]

 # Arcsin normal (AS)
arcsin.mean = arcsin.normal(survey.data, arrond.data, graph.loc)
results$est.mean5 = arcsin.mean[[1]]
results$ll.mean5 = arcsin.mean[[2]]
results$ul.mean5 = arcsin.mean[[3]]

 # Pseudo likelihood (PL)
pseudo.mean = pseudo.likelihood.binomial(survey.data, arrond.data, graph.loc)
results$est.mean6 = pseudo.mean[[1]]
results$ll.mean6 = pseudo.mean[[2]]
results$ul.mean6 = pseudo.mean[[3]]

 # Effective sample size method (ES)
eff.mean = effective.samplesize.binomial(survey.data, arrond.data, graph.loc)
results$est.mean7 = eff.mean[[1]]
results$ll.mean7 = eff.mean[[2]]
results$ul.mean7 = eff.mean[[3]]

 # Proposed method model 1 (M1 - RW1)
mb.m1 = model.based.method1(survey.data, arrond.data, graph.loc, 2500)
results$est.mean8 = apply(mb.m1$posterior.p, 2, quantile, 0.50)
results$ll.mean8 = apply(mb.m1$posterior.p, 2, quantile, 0.025)
results$ul.mean8 = apply(mb.m1$posterior.p, 2, quantile, 0.975)

 # Proposed method model 2 (M2 - RW1)
mb.m2 = model.based.method2(survey.data, arrond.data, graph.loc, 2500)
results$est.mean9 = apply(mb.m2$posterior.p, 2, quantile, 0.50)
results$ll.mean9 = apply(mb.m2$posterior.p, 2, quantile, 0.025)
results$ul.mean9 = apply(mb.m2$posterior.p, 2, quantile, 0.975)

 # Proposed method model 3 (M1 - PS)
mb.m3 = model.based.method3(survey.data, arrond.data, graph.loc, 2500)
results$est.mean10 = apply(mb.m3$posterior.p, 2, quantile, 0.50)
results$ll.mean10 = apply(mb.m3$posterior.p, 2, quantile, 0.025)
results$ul.mean10 = apply(mb.m3$posterior.p, 2, quantile, 0.975)

 # Proposed method model 4 (M2 - PS)
mb.m4 = model.based.method4(survey.data, arrond.data, graph.loc, 2500)
results$est.mean11 = apply(mb.m4$posterior.p, 2, quantile, 0.50)
results$ll.mean11 = apply(mb.m4$posterior.p, 2, quantile, 0.025)
results$ul.mean11 = apply(mb.m4$posterior.p, 2, quantile, 0.975)

	### Plot Results
	#---------------
boxplot(results$est.mean1, results$est.mean2, results$est.mean3, results$est.mean4,
	results$est.mean5, results$est.mean6, results$est.mean7, results$est.mean8, 
	results$est.mean9, results$est.mean10, results$est.mean11,
	cex.axis=0.8, col="lightgrey",
	names=c("UM","HT","NB","LN","AS","PL","ES","M1-RW","M2-RW","M1-PS","M2-PS"))


map <- readShapePoly(".../BT_districts_2007.shp")
map <- map[map@data$CMCNCD=="BE",]
plot(map)

results2 = results[order(results$index),]
plotvar <- results2$est.mean10
nclr <- 8
plotclr <- brewer.pal(nclr,"RdYlGn")
class <- classIntervals(plotvar, nclr, style="fixed",fixedBreaks=seq(0.10,0.50,0.05))
colcode <- findColours(class, plotclr)
plot(map,xlim=c(3560000, 3760000), ylim=c(2540000,2760000))
plot(map, col=colcode, add=T)
legend(3740000, 2780210, legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), 
	cex=2.15, bty="n")

