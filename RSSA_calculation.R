# RSSA_calculation.R
# George Allen, March 14, 2017

######################################################
# script description:
# coming soon.


######################################################
# load required packages:
require(foreign)
require(zyp)
require(RColorBrewer)


######################################################
# functions:

# calculate quantized river surface area at each river observation:
distCalc <- function(csv){
  # calc along stream length at each width:
  nRow = nrow(csv)
  eDiff = abs(csv$utm_east[-1]-csv$utm_east[-nRow])
  nDiff = abs(csv$utm_north[-1]-csv$utm_north[-nRow])
  
  # set length values of large centerline jumps to 30 m:
  eDiff[eDiff>30] = 30
  nDiff[nDiff>30] = 30
  
  # determine length based on summing the x & y position of each vertex:
  diagDist = 30*sqrt(2)
  d = eDiff+nDiff
  d[d == 60] = diagDist
  # add a single length value of 30 to set the length measurements 
  # equal to the number of width measurements:
  if (length(d)<nRow){d = c(d, 30)}
  
  return(d)
}


# error function:
erf = function(A){ 2 * pnorm(A * sqrt(2)) - 1 }


# lognormal maximum likelihood estimation:
lnSLL = function(xL, xU, mu, extrapMode, Amin){
  # erf of lower bin break:
  erf_xu = as.double(erf( (xU-mu) / ( 2*(mu-extrapMode) )^0.5 ))
  
  # erf of upper bin break:
  erf_xl = as.double(erf( (xL-mu) / ( 2*(mu-extrapMode) )^0.5 ))
  
  # erf of Amin:
  erf_xmin = erf( (Amin-mu) / ( 2*(mu-extrapMode) )^0.5 )
  
  # Calculate likelihood:
  likelihood = 0.5 * (erf_xu - erf_xl) / (1 - (0.5 + 0.5*erf_xmin) )
  likelihood = likelihood[-which(likelihood == 0)]
  
  # should sum  to 1:
  #print(sum(unique(likelihood))) 
  
  # return sum log likelihood:
  return(sum(log(likelihood)))
}


# pareto maximum likelihood estimation:
pSLL = function(xL, xU, a, extrapMin, Amin){
  #  prob of upper bin break:
  p_xU = 1-(extrapMin/xU)^a
  
  #  prob of lower bin break:
  p_xL = 1-(extrapMin/xL)^a
  
  #  prob of Amin bin break:
  p_wMin = 1-(extrapMin/Amin)^a
  
  # Calculate likelihood:
  likelihood = (p_xU-p_xL)/(1-p_wMin)
  likelihood = likelihood[-which(likelihood == 0)]
  
  # should sum  to almost 1:
  #print(sum(unique(likelihood))) 
  
  # return sum log likelihood:
  return(sum(log(likelihood)))
  
}

# specialized functions for the pareto fit:
# code based on Elvis's reply on stackexhange forum:
# http://stats.stackexchange.com/questions/78168/how-to-know-if-my-data-fits-pareto-distribution

# distribution, cdf, quantile and random functions for Pareto distributions:
dpareto <- function(x, xm, a) ifelse(x > 0, a*xm**a/(x**(a+1)), 0)
ppareto <- function(q, xm, a) ifelse(q > 0, 1 - (xm/q)**a, 0 )
qpareto <- function(p, xm, a) ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/a))
rpareto <- function(n, xm, a) qpareto(runif(n), xm, a)

# fit pareto distribution to data with MLE:
pareto.mle <- function(x, xm=min(x)){
  alpha = length(x)/(sum(log(x))-length(x)*log(xm))
  return( list(xm = xm, alpha = alpha))
}


# get uncertainty of MLE fit by using an MLE optimization
# algorithm to claculate the hessian (rate at which fit
# falls off when moving away from optimum. Inverse of 
# hessian is the variance-covariance matrix, from which 
# standard deviation may be calculated.
mleOptimizer <- function(p, x){
  -sum(log(p*min(x)^p/(x^(p+1))))
}


# calculate goodness of fit for the distribution fits:
GOF = function(fitParams, h, jA){
  
  GOF = list(X2=vector(), X2_p=vector(), KS_D=vector(), KS_p=vector())
  
  # Pearson's Chi-squared test for count data:
  # low X2 and high p-value signifies a good fit:
  f = dpareto(h$mids, fitParams[[1]], fitParams[[2]][2])
  X2 = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=10)
  GOF[[1]] = X2$statistic
  GOF[[2]] = X2$p.value
  
  # Two sided One sample KS GOF test:
  # low D and high p-value signifies a good fit:
  KS = ks.test(jA, "ppareto",  fitParams[[1]], fitParams[[2]][2], 
               alternative="two.sided")
  GOF[[3]] = KS$statistic
  GOF[[4]] = KS$p.value
  
  return(GOF)
}


extrapSA_calculator = function(pFit, fOaMeans, Amin, Amax, sumA){
  
  pSA = list(fOWsd=vector(), pFitSD=vector())
  
  a = gpFit$alpha # use global mean alpha values. ##pFit$alpha
  xm = fOaMeans[2]
  x1 = fOaMeans
  x2 = Amin
  x3 = Amax
  
  # calc. area under curve with different 1st order stream widths:
  obsIntegral = -(a[2]*x3*(xm/x3)^a[2])/(a[2]-1) + (a[2]*x2*(xm/x2)^a[2])/(a[2]-1)
  extrapIntegral = -(a[2]*x2*(xm/x2)^a[2])/(a[2]-1) + (a[2]*x1*(xm/x1)^a[2])/(a[2]-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  
  # add observed to estimated and convert to km2:
  pSA$fOWsd = (sumA*extrap2ObsRatio+sumA)*reducer*1e-6
  
  # calc. area under curve with different pareto MLE fits:
  obsIntegral = -(a*x3*(xm/x3)^a)/(a-1) + (a*x2*(xm/x2)^a)/(a-1)
  extrapIntegral = -(a*x2*(xm/x2)^a)/(a-1) + (a*x1[2]*(xm/x1[2])^a)/(a-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  
  # add observed to estimated and convert to km2:
  pSA$pFitSD = (sumA*extrap2ObsRatio+sumA)*reducer*1e-6
  
  return(pSA)
}


# plot histogram and fit for each basin (Fig 2b):
histPlot = function(fOaMeans, Amax, pFit, gpFit, sumA, breaks, h, Alen, int){
  # calculate total area of the bins: 
  pTotA = Alen*int
  
  
  # plot empirical frequency histogram:
  xlim = c(fOaMeans[1], Amax)
  ylim = c(1, max(dpareto(fOaMeans[1], pFit[[1]], pFit[[2]][3])*pTotA,
                  dpareto(fOaMeans[1], gpFit[[1]], gpFit[[2]][2])*pTotA))
  
  title = paste0('i: ',i, '   basin: ',fN[i])
  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab=paste("Area (m2)"), 
       ylab="N measurements",
       main=title,
       type='n', 
       log='xy',
       xaxt='n')
  xAxis = axis(1, labels=F)
  axis(1, at=xAxis, labels=xAxis*reducer)
  
  mtext(paste("Nobs:", Alen, "   %SAobs: ", round(obs_prcnt, 4)))
  rect(breaks[-length(breaks)], 1, 
       breaks[-1], h$counts, 
       border=NA, col="gray")
  
  
  
  # basin-fit pareto curve:
  
  # 1st order width uncertainty:
  x1 = fOaMeans[1]
  x2 = fOaMeans[3]
  y1 = dpareto(fOaMeans[1], pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[3], pFit[[1]], pFit[[2]][2])*pTotA
  polygon(x = c(x1, x1, x2, x2),
          y = c(1, y1, y2, 1),
          col=rgb(1,0,0,0.2), border=NA)
  
  # Pareto MLE fit uncertainty: 
  y1 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][1])*pTotA
  y2 = dpareto(Amax, pFit[[1]], pFit[[2]][1])*pTotA
  y3 = dpareto(Amax, pFit[[1]], pFit[[2]][3])*pTotA
  y4 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][3])*pTotA
  polygon(x = c(fOaMeans[2], Amax, Amax, fOaMeans[2]),
          y = c(y1, y2, y3, y4),
          col=rgb(1,0,0,0.2), border=NA)
  
  # estimated area polygon:
  y1 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(Amin, pFit[[1]], pFit[[2]][2])*pTotA
  polygon(x = c(fOaMeans[2], Amin, Amin, fOaMeans[2]),
          y = c(y1, y2, 1, 1),
          col=rgb(0,0.2,1,0.2), border=NA)
  
  # Pareto fit over observed data:
  y1 = dpareto(Amin, pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(Amax, pFit[[1]], pFit[[2]][2])*pTotA
  segments(Amin, y1, Amax, y2, col=2, lwd=1.2)
  
  # extrapolated Pareto fit:
  y1 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(Amin, pFit[[1]], pFit[[2]][2])*pTotA
  segments(fOaMean, y1, Amin, y2, col=2, lwd=1.2, lty=2)
  
  # first order uncertainty bars:
  y1 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[3], pFit[[1]], pFit[[2]][2])*pTotA
  arrows(fOaMeans[2], y1, fOaMeans[3], y2, 0.05, 90, col=2, lwd=1.2)
  y1 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[1], pFit[[1]], pFit[[2]][2])*pTotA
  arrows(fOaMeans[2], y1, fOaMeans[1], y2, 0.05, 90, col=2, lwd=1.2)
  
  y2 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][1])*pTotA
  arrows(fOaMeans[2], 1, fOaMeans[1], 1, 0.05, 90, col=2, lwd=1.2)
  y2 = dpareto(fOaMeans[2], pFit[[1]], pFit[[2]][3])*pTotA
  arrows(fOaMeans[2], 1, fOaMeans[3], 1, 0.05, 90, col=2, lwd=1.2)
  
  # global-fit pareto curve:
  
  # 1st order width uncertainty:
  x1 = fOaMeans[1]
  x2 = fOaMeans[3]
  y1 = dpareto(fOaMeans[1], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[3], gpFit[[1]], gpFit[[2]][2])*pTotA
  polygon(x = c(x1, x1, x2, x2),
          y = c(1, y1, y2, 1),
          col=rgb(1,1,0,0.2), border=NA)
  
  # Pareto MLE fit uncertainty: 
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][1])*pTotA
  y2 = dpareto(Amax, gpFit[[1]], gpFit[[2]][1])*pTotA
  y3 = dpareto(Amax, gpFit[[1]], gpFit[[2]][3])*pTotA
  y4 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][3])*pTotA
  polygon(x = c(fOaMeans[2], Amax, Amax, fOaMeans[2]),
          y = c(y1, y2, y3, y4),
          col=rgb(1,0,0,0.2), border=NA)
  
  # estimated area polygon:
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(Amin, gpFit[[1]], gpFit[[2]][2])*pTotA
  polygon(x = c(fOaMeans[2], Amin, Amin, fOaMeans[2]),
          y = c(y1, y2, 1, 1),
          col=rgb(0,0.2,1,0.2), border=NA)
  
  # Pareto fit over observed data:
  y1 = dpareto(Amin, gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(Amax, gpFit[[1]], gpFit[[2]][2])*pTotA
  segments(Amin, y1, Amax, y2, col=4, lwd=1.2)
  
  # extrapolated Pareto fit:
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(Amin, gpFit[[1]], gpFit[[2]][2])*pTotA
  segments(fOaMean, y1, Amin, y2, col=4, lwd=1.2, lty=2)
  
  # first order uncertainty bars:
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[3], gpFit[[1]], gpFit[[2]][2])*pTotA
  arrows(fOaMeans[2], y1, fOaMeans[3], y2, 0.05, 90, col=4, lwd=1.2)
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(fOaMeans[1], gpFit[[1]], gpFit[[2]][2])*pTotA
  arrows(fOaMeans[2], y1, fOaMeans[1], y2, 0.05, 90, col=4, lwd=1.2)
  
  y2 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][1])*pTotA
  arrows(fOaMeans[2], 1, fOaMeans[1], 1, 0.05, 90, col=4, lwd=1.2)
  y2 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][3])*pTotA
  arrows(fOaMeans[2], 1, fOaMeans[3], 1, 0.05, 90, col=4, lwd=1.2)
  
  # add pareto statistics:
  legend("topright", paste0(
    "SA: ", round(pSA$totSD[2], 1), " km2 (", round(pSA$totSD[1],1),  
      "-", round(pSA$totSD[3],1), " km2)\n", 
    "%SA: ", round(pSA_prcnt[2], 2), "% (", round(pSA_prcnt[3], 2),  
      "-", round(pSA_prcnt[1], 2), "%)\n", 
    "alpha: ", round(pFit$alpha[2], 3), "(", round(pFit$alpha[3], 3),  
      "-", round(pFit$alpha[1], 3), ")\n",  "\n",
    "X2: ", round(pGOF$X2, 3), "    p: ", round(pGOF$X2_p, 3), "\n",
    "K-S D: ", round(pGOF$KS_D, 2), "    p: ", round(pGOF$KS_p, 3)),
    xjust=0, text.col=4, bty="n", cex=0.7)
}


######################################################
# hard-coded variables:

# file paths:
outD = 'H:/2017_02_08_GRWL_Nature_Manuscript/GRWL/GRWL_by_hydroBASIN/'
pdfOut = 'H:/2017_02_08_GRWL_Nature_Manuscript/figs/fig2/fig2b_SAfits_NobsGt10000.pdf'
hydroBASINpath = 'H:/2017_02_08_GRWL_Nature_Manuscript/GRWL/basin_shapefiles/hydroBASINs/hybas_allMain.dbf'
gTabP = 'H:/2017_02_08_GRWL_Nature_Manuscript/figs/fig2/hydroBasin_SAfitStatistics_minW90_NobsGt10000.csv'
GRWLpath = "H:/2017_02_08_GRWL_Nature_Manuscript/GRWL/basin_shapefiles/hydroBASINs/hybas_allMain.dbf"
aridityPath = 'H:/2017_02_08_GRWL_Nature_Manuscript/climate/aridity/aridityByBasin.dbf'

# list basin-by-basin river width file paths:
fP = list.files(outD, 'csv', full.names=T)
fN = list.files(outD, 'csv')
fN = sub('.csv', '', fN)
hBASIN_code = sub("GRWL_", '', fN)


# user-defined constants:
# convert are units to reduce the size of numbers 
# to prevent overflow on MLE statistical tests:
reducer = 100

# mean pixel length:
mL = mean(c(30, 30*sqrt(2)))/reducer 
minL = 30/reducer
int = 100*mL  # histogram binning interval
wMin = 90 # minimum GRWL width
Amin = wMin*mL # minimum value to bin
minNobs = 1e3 # minimum number of GRWL observations in a basin for analysis
minElev = 0 # remove all rivers below this elevation (m) 
#fitMax = mL*2e3

# define minimum extrapolation bounds:
fOaMean = 0.321*mL # mean of the median widths of 1st order streams
fOaSD = 0.077*mL # sd width of 1st order streams
fOaMeans = c(fOaMean-fOaSD, fOaMean, fOaMean+fOaSD)

# create stats output table to store data:
gTabNames = c("hBASIN_code", 
              "basinA_km2", 
              "nObs", 
              "nChan_mean",
              "obSA_km2",
              "obSA_gt100_km2", 
              "pSA_km2",
              "pSA_rangeKm2",
              "pSA_prcnt",
              "pSA_prcnt_range",
              "fOw_uncertaintyKm2",
              "alpha_uncertaintyKm2",
              "pXm",
              "pAlpha", 
              "pAlpha_range", 
              "pX2_stat",
              "pX2_p",
              "pKS_D",
              "pKS_p")
gTab = data.frame(array(NA, c(1, length(gTabNames))))
names(gTab) = gTabNames

# global MLE fit:
gpFit = list(xm=Amin, 
             alpha=c(1.029474, 1.029690, 1.029905), 
             stdev=0.000215681)





######################################################
# plot width distributions by hBASIN basin:

# entire GRWL database:
# fP = 'H:/2017_02_08_GRWL_Nature_Manuscript/figs/fig1/analysis/allGRWL.csv'
#globalLandArea = 132773914

# read in hBASIN shapefile dbf:
hBASIN = foreign::read.dbf(sub('hybas_allMain', 'hybas_allMainCopy', hydroBASINpath))

pdf(pdfOut, width=6, height=6)
m = 1
for (i in 1:length(fP)){ # Amazon: i = 5687 # Miss: i=6829
  
  # read in widths:
  csv = read.csv(fP[i], header=T)
 
  # remove data with widths<100m, elev<1m, lakes, canals:
  keep = csv$width_m>wMin & 
    csv$elev_m>minElev & 
    csv$lakeFlag!=1 & 
    csv$lakeFlag!=3 
  
  # skip basin if there are less than the minimum number of observations:
  if (length(which(keep))<minNobs){ next }
  
  print(paste('i =',i, ' of ', length(fP),', fN = ', fN[i]))
  
  # calc width, distance, river area, sum river area, maximum area value: 
  w = csv$width_m[keep]/reducer
  d = distCalc(csv[keep,])
  A_raw = w*d
  A = A_raw[A_raw>Amin]
  sumA = sum(A)
  Amax = max(A)
  Alen = length(A)
  
  ####################################################
  # MLE fit: 
  pFit = pareto.mle(A) 
  # standard deviation of the MLE fit:
  #pFit$stdev = sqrt(pFit$alpha/Alen)
  
  # get uncertainty of fit by using an MLE optimization
  # algorithm to claculate the hessian (rate at which fit
  # falls off when moving away from optimum. Inverse of 
  # hessian is the variance-covariance matrix, from which 
  # standard deviation may be calculated.
  hessian = suppressWarnings(nlm(mleOptimizer, p=1, A, hessian=T))$hessian
  pFit$stdev = sqrt(diag(solve(hessian)))
  pFit$alpha = c(pFit$alpha-pFit$stdev, pFit$alpha, pFit$alpha+pFit$stdev)

  
  # calculate GOF statistics with X2 and KS test:
  breaks = seq(Amin, Amax+int, int)
  h = hist(A, breaks, plot=F)
  jA = jitter(A) # to prevent ties in the following GOF tests
  pGOF = GOF(gpFit, h, jA)

  
  ####################################################
  # total surface area estimate:
  
  # calculate total surface area of rivers and streams by extendending
  # area abundance to smaller rivers:
  pSA = extrapSA_calculator(pFit, fOaMeans, Amin, Amax, sumA)
  
  # combine uncertainty from MLE fit and extrapolation minimum
  pSA$totSD = c(
    pSA$fOWsd[2]-sqrt((pSA$fOWsd[2]-pSA$fOWsd[1])^2 + (pSA$pFitSD[2]-pSA$pFitSD[1])^2), 
    pSA$fOWsd[2], 
    pSA$fOWsd[2]+sqrt((pSA$fOWsd[2]-pSA$fOWsd[3])^2 + (pSA$pFitSD[2]-pSA$pFitSD[3])^2))
  
 
  # calculate the % of land surface occupied by rivers and streams: 
  pSA_prcnt = 100*pSA$totSD/hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]] #132773914
  obs_prcnt = 100*sumA*reducer*1e-6/hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]] #132773914
  
  print(paste("a =", round(pFit$alpha[2],3), "   ", round(100*obs_prcnt/pSA_prcnt[2], 2), 
              "% observed"))
  # Plot histogram and Pareto fits (Fig. 2b):
  histPlot(fOaMeans, Amax, pFit, gpFit, sumA, breaks, h,  Alen, int)
  
  
  ####################################################
  # add data to statsTable:
  gTab[m, ] = as.vector(c(fN[i], 
                          round(hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]]),
                          Alen, 
                          mean(csv$nchannels[csv$width_m>(Amin/mL) & 
                                               csv$elev_m>minElev & 
                                               csv$lakeFlag!=1 & 
                                               csv$lakeFlag!=3]),
                          round(sum(A_raw)*mL*reducer*1e-6), 
                          round(sumA*reducer*1e-6),
                          round(pSA$totSD[2]), 
                          paste0(pSA$totSD[1], '-', pSA$totSD[3]),
                          round(pSA_prcnt[2], 4), 
                          paste0(pSA_prcnt[1], '-', pSA_prcnt[3]),
                          paste0(pSA$fOWsd[1], '-', pSA$fOWsd[3]),
                          paste0(pSA$pFitSD[1], '-', pSA$pFitSD[3]),
                          round(pFit$xm), 
                          round(pFit$alpha[1], 3),
                          paste0(round(pFit$alpha[1], 3), '-', round(pFit$alpha[3], 3)),
                          round(pGOF$X2, 2), 
                          round(pGOF$X2_p, 3), 
                          round(pGOF$KS_D, 2), 
                          round(pGOF$KS_p, 3)
                        ))
  
  
  
  m = m + 1

}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)


# write out statistics table:
write.csv(gTab, gTabP, row.names=F)


x = strsplit(as.character(gTab$pSA_rangeKm2), '-')
lwr = sapply(x, '[[', 1)
upr = sapply(x, '[[', 2)
print(sum(as.numeric(gTab$pSA_km2)))
print(sum(as.numeric(lwr)))
print(sum(as.numeric(upr)))

print(sum(as.numeric(gTab$obSA_km2)))
print(paste(round(100*sum(as.numeric(gTab$obSA_km2))/sum(as.numeric(gTab$pSA_km2)), 1), 
            "% observed area"))
print(mean(as.numeric(gTab$pKS_D)))



####################################################
# Attach gtab attributes to hydroBASIN shapefile:

io = 0; if (io == 1){
  #print(gTab)
  
  gTab = read.csv(gTabP, header=T)
  if (length(grep("dbf", names(gTab)))>0){ gTab = gTab$dbf }
  
  # read in hBASIN dbf and attached gtab attributes to it:
  hBASIN = foreign::read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
  
  bindTab = data.frame(array(NA, c(nrow(hBASIN), ncol(gTab))))
  names(bindTab) = names(gTab)
  matchInd = match(gTab$hBASIN_code, hBASIN$MAIN_BAS)
  
  gMat = matrix(as.matrix(gTab), ncol = ncol(gTab), dimnames = NULL)
  for (i in 1:length(matchInd)){
    bindTab[matchInd[i], ] = gMat[i,]
  }
  
  # format table for arcmap shapefile:
  newhBASIN = as.data.frame(cbind(hBASIN, bindTab))
  newhBASIN[is.na(newhBASIN)] = 0
  newhBASIN = data.matrix(newhBASIN)
  
  write.dbf(newhBASIN, hydroBASINpath)
  
}




# sum extrap surface area
print(sum(as.numeric(gTab$pSA_km2)))

# read in hBASIN dbf and attached gtab attributes to it:
hBASIN = foreign::read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))

bindTab = data.frame(array(NA, c(nrow(hBASIN), ncol(gTab))))
names(bindTab) = names(gTab)
matchInd = match(gTab$hBASIN_code, hBASIN$MAIN_BAS)

gMat = matrix(as.matrix(gTab), ncol = ncol(gTab), dimnames = NULL)
for (i in 1:length(matchInd)){
  bindTab[matchInd[i], ] = gMat[i,]
}

# format table for arcmap shapefile:
newhBASIN = as.data.frame(cbind(hBASIN, bindTab))
newhBASIN[is.na(newhBASIN)] = 0
newhBASIN = data.matrix(newhBASIN)


# read in tables:
GRWL = as.data.frame(newhBASIN)
aridity = foreign::read.dbf(aridityPath)

# add FID column to GRWL table to match up data:
if (!'FID' %in% names(GRWL)){ FID = 1:nrow(GRWL)-1; GRWL = cbind(FID, GRWL) }

# only use basins with areas larger than 100k km2 to remove small basins
# that have a %SA value of 0 just because they are small:
k = GRWL$basinA_km2>1e5
length(k)
#k = order(GRWL$basinA_km2, decreasing=T)[1:20]

x1 = aridity$MEAN[match(GRWL$FID[k], aridity$FID_)]
x2 = GRWL$basinA_km2[k]
y = GRWL$pSA_prcnt[k]
w = GRWL$basinA_km2[k]

# fit a weighted multiple linear regression to aridity, basin size, and %SA data:
fit = lm(log(y)~log(x1)+log(x2), weights=w); print(summary(fit))


fitFun <- function(x1, x2, fit){
  return(exp(predict(fit, data.frame(x1=x1, x2=x2))))
}


# plots to pdf:
pdfOut = 'H:/temp/temp001.pdf'
pdf(pdfOut,  width=9, height=4)
layoutTab = rbind(c(1,1,2,2,3))
layout(layoutTab)

tab = data.frame(x1, x2, y, w)

with(tab, symbols(x1, y, circles=w, inches=0.3, bg=rgb(0,0,0,0.5), fg=NA, 
                  main="ED Figure 2a", 
                  xlab="Aridity Index", ylab = "Percent rivers and streams",
                  las=1))
xMin = 0; xMax = max(x1)
xSeq1= xMin+((0:50)^10/50^10)*(xMax-xMin)
xMin = 0; xMax = max(x2)
xSeq2= xMin+((0:50)^10/50^10)*(xMax-xMin)
lines(xSeq1, fitFun(x1=xSeq1, x2=xSeq2, fit), lwd=1.7)

text(xMin, max(y),  
     paste0("ln(%SA) = ", round(fit[[1]][[1]]*fit[[1]][[2]], 1), "ln(AI) + ",
            round(fit[[1]][[1]]*fit[[1]][[3]], 1), "ln(BA)"),
     pos=4)


with(tab, symbols(x2, y, circles=w, inches=0.3, bg=rgb(0,0,0,0.5), fg=NA, 
                  main="ED Figure 2b", 
                  xlab="Basin Size (km2)", ylab = "Percent rivers and streams",
                  las=1))
lines(xSeq2, fitFun(x1=xSeq1, x2=xSeq2, fit), lwd=1.7)
text(xMin, max(y),  
     paste0("ln(%SA) = ", round(fit[[1]][[1]]*fit[[1]][[2]], 1), "ln(AI) + ",
            round(fit[[1]][[1]]*fit[[1]][[3]], 1), "ln(BA)"),
     pos=4)


#display.brewer.all()
#cols = cols[GRWL$area_km2]
#col = brewer.pal(n=9,name="YlGnBu")
#plot(x, y, 
#     xlab="Aridity Index", ylab = "Percent rivers and streams",
#     main="ED Figure 2", log='xy',#axes=F, 
#     col=col, pch=16, cex=0.8, las=1)


# add legend:
dotSizeSeq = c(1e5, 2.5e5, 5e5, 1e6, 2e6, 6e6)
legendTab = data.frame(x=rep(1, 6), y=c(1:6), dotSizeSeq) 
suppressWarnings(with(legendTab, symbols(x, y, circles=dotSizeSeq, 
                                         inches=0.3,  bg=rgb(0,0,0,0.5), fg=NA, 
                                         main="Basin Area (BA)",
                                         axes=F, xlab='', ylab='')))


options(scipen=2)
text(rep(1, 4), c(1:4), paste(dotSizeSeq, "km2"), pos=1, offset=2.5)

dev.off() 
cmd = paste('open', pdfOut)
system(cmd)

summary(fit)


# use climate-%SA regression to fill in missing basins in hydroBASINs shapefile:

# set up model X and Y:
ai = aridity$MEAN[match(GRWL$FID, aridity$FID_)]
perSA = GRWL$pSA_prcnt
basinA = GRWL$area_km2
modSA = perSA

# for basins with no surface area estimates, use climate-SA model
# to estimate %SA: 
modSA[perSA==0] = fitFun(ai[perSA==0], GRWL$area_km2[perSA==0], fit)
modSA[is.na(modSA)] = 0

# set greenland icesheet where no rivers were observed to zero:
modSA[substr(GRWL$MAIN_BAS, 1,1) == "9" & perSA==0] = 0

range(modSA)

if ("modSA" %in% names(GRWL)){GRWL = GRWL[ ,names(GRWL)!="modSA"]}

summary(fit)
sum(GRWL$area_km2*modSA/100)
length(which(GRWL$pSA_km2 != 0))
mean(GRWL$pKS_D)
mean(GRWL$pAlpha[GRWL$pSA_km2 > 0])











