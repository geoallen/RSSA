
##############################################################################
# fig3_SAfit.R
# George Allen, March 14, 2017

# Description: coming soon


##############################################################################
# load packages
##############################################################################
if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"MASS" %in% rownames(installed.packages())){
  install.packages("MASS")}; require(MASS)


##############################################################################
# hard-coded variables
##############################################################################

# specify root paths:
wd = "H:/2018_03_16_GRWL_Manuscript/git/RSSA/"
outPath = paste0(wd, 'output/GRWL/GRWL_by_hydroBASIN/')
pdfOut =  paste0(wd, 'output/figs/fig3_SAfits_20basins.pdf')
gTabP = paste0(wd, 'output/figs/fig3_SAfits_20basins.csv')
hydroBASINpath = paste0(wd, 'input/basin_shapefiles/hydroBASINs/hybas_allMain.dbf')
valPath = paste0(wd, 'input/validation/Database_S1_GRWL_validation_data.csv')

# get file paths & codes: 
fP = list.files(outPath, 'csv', full.names=T)
fN = list.files(outPath, 'csv')
fN = sub('.csv', '', fN)
hBASIN_code = sub("GRWL_", '', fN)


# convert area units to reduce the size of numbers 
# to prevent overflow on MLE statistical tests:
reducer = 100
mL = mean(c(30, 30*sqrt(2)))/reducer 
minL = 30/reducer
int = 100*mL  # plotting histogram binning interval
wMin = 90
Amin = wMin*mL # minimum value to bin
minNobs = 2.5e5  # minimum number of GRWL observations in a basin for analysis
minElev = 0 # include rivers above this elevation (m)

# define minimum extrapolate bounds from Allen et al. Nature Communications
# DOI: s41467-018-02991-w
fOaMean = 0.321*mL # mean of the median widths of 1st order streams
fOaSD = 0.077*mL # sd width of 1st order streams
fOaMeans = c(fOaMean-fOaSD, fOaMean, fOaMean+fOaSD)

# define number of ensemble runs to estimate error:
nRun = 3

# create Stats tab to store data:
gTabNames = c("hBASIN_code", 
              "basinA_km2", 
              "nObs", 
              "nChan_mean",
              "Amax",
              "obSA_km2",
              "obSA_gt90m_km2", 
              # global
              "gSA_km2",
              "gSA_sd_km2",
              "gSA_pcnt",
              "gSA_sd_pcnt",
              "gXmin",
              "gAlpha", 
              "gAlpha_sd", 
              # basin
              "pSA_km2",
              "pSA_sd_km2",
              "pSA_pcnt",
              "pSA_sd_pcnt",
              "pXmin",
              "pAlpha", 
              "pAlpha_sd", 
              "pX2_stat",
              "pX2_p",
              "pKS_D",
              "pKS_p"
              )

globalLandArea = 132773914

gTab = data.frame(array(NA, c(1, length(gTabNames))))
names(gTab) = gTabNames

# global MLE fit:
gpFit = list(xm=Amin, 
             #alpha=1.01736, # minN = 65,000
             #stdev=0.0217053) # minN = 65,000
             #alpha=1.029690, # global fit
             #stdev=0.000215681) # global fit
             #stdev=0.1148863) # minN = 100000
             #alpha=0.9508378, # minN = 100,000
             #alpha=0.919, # 20 basin with most obs, no min elev
             #stdev=0.07176863) # 20 basins with most obs, no min elev
             alpha=0.90345, # 20 basin with most obs, min elev = 0
             stdev=0.06285026) # 20 basins with most obs, min elev = 0



##############################################################################
# functions
##############################################################################

# calculate quantized river surface area at each river observation:
distCalc <- function(csv){
  # calc along stream length at each width:
  nRow = nrow(csv)
  eDiff = abs(csv$utm_east[-1]-csv$utm_east[-nRow])
  nDiff = abs(csv$utm_north[-1]-csv$utm_north[-nRow])
  
  # set length values of large centerline jumps to 30m:
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


# specialized functions for the pareto fit:
# code based on Elvis's reply on stackexhange forum:
# http://stats.stackexchange.com/questions/78168/how-to-know-if-my-data-fits-pareto-distribution

# distribution, cdf, quantile and random functions for Pareto distributions:
dpareto <- function(x, xm, a) ifelse(x > 0, a*xm**a/(x**(a+1)), 0)
ppareto <- function(q, xm, a) ifelse(q > 0, 1 - (xm/q)**a, 0 )
qpareto <- function(p, xm, a) ifelse(p < 0 | p > 1, NaN, xm*(1-p)**(-1/a))
rpareto <- function(n, xm, a) qpareto(runif(n), xm, a)


# fit pareto distribution to data with MLE:
pareto.mle <- function(x, std = FALSE){
  
  # calculate alpha using calculus:
  xm = min(x)
  alpha = length(x)/(sum(log(x))-length(x)*log(xm))
  
  # std turns on calculating standard deviaton of MLE fit: 
  if (std == TRUE){
    # get uncertainty of fit by using an MLE optimization
    # algorithm to claculate the hessian (rate at which fit
    # falls off when moving away from optimum. Inverse of 
    # hessian is the variance-covariance matrix, from which 
    # standard deviation may be calculated.
    hessian = suppressWarnings(nlm(mleOptimizer, p=1, A, hessian=T))$hessian
    stdev = sqrt(diag(solve(hessian)))
  }else{
    stdev = NA
  }
  
  return( list(xm = xm, alpha = alpha, stdev=stdev))
  
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
GOF = function(pFit, h, jA){
  
  # Pearson's Chi-squared test for count data:
  # low X2 and high p-value signifies a good fit:
  f = dpareto(h$mids, pFit$xm, pFit$alpha)
  X2test = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=10)
  
  # Two sided One sample KS GOF test:
  # low D and high p-value signifies a good fit:
  KStest = ks.test(jA, "ppareto",  pFit$xm, pFit$alpha, alternative="two.sided")
  
  return(list(
    X2 = X2test$statistic,
    X2_p = X2test$p.value,
    KS_D = KStest$statistic[[1]],
    KS_p = KStest$p.value
  ))
}


extrapSA_calculator = function(pFit, fOaMeans, Amin, Amax, sumA){
  
  a = c(pFit$alpha-pFit$stdev, pFit$alpha, pFit$alpha+pFit$stdev)
  xm = fOaMeans[2]
  x1 = fOaMeans
  x2 = Amin
  x3 = Amax
  
  # calc. area under curve with different 1st order stream widths:
  obsIntegral = -(a[2]*x3*(xm/x3)^a[2])/(a[2]-1) + (a[2]*x2*(xm/x2)^a[2])/(a[2]-1)
  extrapIntegral = -(a[2]*x2*(xm/x2)^a[2])/(a[2]-1) + (a[2]*x1*(xm/x1)^a[2])/(a[2]-1)
  allIntegral = -(a[2]*x2*(xm/x2)^a[2])/(a[2]-1) + (a[2]*x1*(xm/x1)^a[2])/(a[2]-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  
  # add observed to estimated and convert to km2:
  fOWsd = (sumA*extrap2ObsRatio+sumA)*reducer*1e-6
    
  # calc. area under curve with different pareto MLE fits:
  obsIntegral = -(a*x3*(xm/x3)^a)/(a-1) + (a*x2*(xm/x2)^a)/(a-1)
  extrapIntegral = -(a*x2*(xm/x2)^a)/(a-1) + (a*x1[2]*(xm/x1[2])^a)/(a-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  
  # add observed to estimated and convert to km2:
  alphaSD = (sumA*extrap2ObsRatio+sumA)*reducer*1e-6

  return(list(
    fOWsd = fOWsd,
    alphaSD = alphaSD[!is.na(alphaSD)]
  ))
  
}


# combine uncertainty from MLE fit and extrapolation minimum
# they are multiplied together so...
mCarloUncertProp <- function(pFit, fOaMean, fOaSD, Amin, Amax, N, sumA){
  
  # calculate total surface area:
  a = pFit$alpha
  xm = fOaMean
  x1 = fOaMean
  x2 = Amin
  x3 = Amax
  
  # calc. area under curve with different 1st order stream widths:
  obsIntegral = -(a*x3*(xm/x3)^a)/(a-1) + (a*x2*(xm/x2)^a)/(a-1)
  extrapIntegral = -(a*x2*(xm/x2)^a)/(a-1) + (a*x1*(xm/x1)^a)/(a-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  # add observed to estimated area:
  meanSA = sumA*extrap2ObsRatio+sumA
  
  # estimate uncertainty using Monte Carlo error propagation:
  if (is.na(pFit$stdev)){
    a = pFit$alpha
  }else{
    a = rnorm(N, pFit$alpha, pFit$stdev)
  }
  
  xm = fOaMean
  x1 = rnorm(N, fOaMean, fOaSD)
  
  # calc. area under curve with different 1st order stream widths:
  obsIntegral = -(a*x3*(xm/x3)^a)/(a-1) + (a*x2*(xm/x2)^a)/(a-1)
  extrapIntegral = -(a*x2*(xm/x2)^a)/(a-1) + (a*x1*(xm/x1)^a)/(a-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  # add observed to estimated area:
  simSA = sumA*extrap2ObsRatio+sumA
  
  return(list(
    meanSimSA =  meanSA*reducer*1e-6,
    sdSimSA = sd(simSA, na.rm=T)*reducer*1e-6
  ))
  
}


# plot histogram and fit for each basin (Fig 2b):
histPlot = function(fOaMeans, Amax, pFit, gpFit, sumA, 
                    breaks, h, Alen, int, pMcarlo, gpMcarlo, 
                    obs_prcnt, pSA_prcnt, gpSA_prcnt, pGOF, m){
  
  # calculate total area of the bins: 
  pTotA = Alen*int

  # get pFit formated for the plotter (script could be improved):
  if (is.na(pFit$stdev)){ pFit$stdev = 0 }
  pFit$alpha = c(pFit$alpha-pFit$stdev, pFit$alpha, pFit$alpha+pFit$stdev)
  
  if (is.na(gpFit$stdev)){ gpFit$stdev = 0 }
  gpFit$alpha = c(gpFit$alpha-gpFit$stdev, gpFit$alpha, gpFit$alpha+gpFit$stdev)  
  
  # skip the first plot
  #if (m == 1){plot.new()}
  
  # plot empirical frequency histogram:
  xlim = c(Amin, 10725)
  ylim = c(1, 2e6)
  
  title = paste0('Basin #: ',fN[i])
  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab=paste("Area (m2)"), 
       ylab="N measurements",
       #main=title,
       type='n', 
       log='xy',
       xaxt='n',
       yaxt='n')
  if(m %in% c(16:20)){
  xAxis = axis(1, labels=F)
  axis(1, at=xAxis, labels=xAxis*reducer)
  }
  if(m %in% c(1,6,11,16)){
    yAxis = axis(2, labels=F)
    axis(2, at=yAxis, labels=yAxis)
  }
  mtext(paste(title), line = -1, cex=0.75)
  rect(breaks[-length(breaks)], 1, 
       breaks[-1], h$counts, 
       border=NA, col="gray")
  
 
  
  # basin-fit pareto curve:
  io=0; if(io==1){
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
          col=rgb(1,0.2,0,0.1), border=NA)
  }
  # Pareto fit over observed data:
  y1 = dpareto(Amin, pFit[[1]], pFit[[2]][2])*pTotA
  y2 = dpareto(Amax, pFit[[1]], pFit[[2]][2])*pTotA
  segments(Amin, y1, Amax, y2, col=4, lwd=1.2)
  
  io=0;if (io==1){
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
          col=rgb(0,0,1,0.2), border=NA)
  
  # Pareto MLE fit uncertainty: 
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][1])*pTotA
  y2 = dpareto(Amax, gpFit[[1]], gpFit[[2]][1])*pTotA
  y3 = dpareto(Amax, gpFit[[1]], gpFit[[2]][3])*pTotA
  y4 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][3])*pTotA
  polygon(x = c(fOaMeans[2], Amax, Amax, fOaMeans[2]),
          y = c(y1, y2, y3, y4),
          col=rgb(0,0,1,0.2), border=NA)
  
  # estimated area polygon:
  y1 = dpareto(fOaMeans[2], gpFit[[1]], gpFit[[2]][2])*pTotA
  y2 = dpareto(Amin, gpFit[[1]], gpFit[[2]][2])*pTotA
  polygon(x = c(fOaMeans[2], Amin, Amin, fOaMeans[2]),
          y = c(y1, y2, 1, 1),
          col=rgb(0,0.2,1,0.1), border=NA)
  
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
  }
  
  
  # add pareto statistics:
  io=1;if (io==1){
  legend("topright", #c(paste0("global fit: \n",
    #" SA: ", round(gpMcarlo$meanSimSA, 1), " +|- ", round(gpMcarlo$sdSimSA,1), " km2\n", 
    #" %SA: ", round(gpSA_prcnt, 2), " +|- ", round(gpSA_sd_prct, 2),"%\n", 
    #" alpha: ", round(gpFit$alpha[2], 3), " +|- ", round(gpFit$stdev, 3), "\n"),
    paste0("basin fit: \n",
    #" SA: ", round(pMcarlo$meanSimSA, 1), " +|- ", round(pMcarlo$sdSimSA,1), " km2\n", 
    #" %SA: ", round(pSA_prcnt, 2), " +|- ", round(pSA_sd_prct, 2),"%\n",
    " alpha: ", round(pFit$alpha[2], 3), " +|- ", round(pFit$stdev, 3), "\n",  
    " K-S D: ", round(pGOF$KS_D, 2), "    p: ", round(pGOF$KS_p, 3)),
    xjust=0, text.col=c(4,2), bty="n", cex=0.7)
  }
  
}


##############################################################################
# RivWidth-based error estimate
##############################################################################

# read in error spreadsheet
valTab = read.csv(valPath, header=T, sep = '\t')
keep = valTab$GAUGE_WIDTH_M >= 90
mod = lm(valTab$GRWL_WIDTH_M[keep] ~ valTab$GAUGE_WIDTH_M[keep])
#plot(valTab$GAUGE_WIDTH_M[keep], valTab$GRWL_WIDTH_M[keep])
#abline(mod, col=2)
dif = valTab$GRWL_WIDTH_M[keep] - valTab$GAUGE_WIDTH_M[keep]

# fit distribution with MLE:
nfit = fitdistr(dif, "normal")$estimate

lineInt = 1
lineSeq = seq( floor(min(dif)/lineInt)*lineInt, 
               ceiling(max(dif)/lineInt)*lineInt, 
               lineInt)
hBin = 5
hSeq = seq( floor(min(dif)/hBin)*hBin, 
            ceiling(max(dif)/hBin)*hBin, 
            hBin)

dn = dnorm(lineSeq, nfit[1], nfit[2])
fn = dn*length(dif)*hBin

# plot difference between GRWL and in situ:
plot(dif); abline(h=0, col=2)
# plot histograms of differences between GRWL and in situ
# to show a roughly normal distribution:
h = hist(dif, breaks=hSeq, 
     main = "Difference Histogram",
     xlab="GRWL Width - In Situ Width (m)")
#lines(lineSeq, fn, col=4)
abline(v=nfit[1], col=4, lwd=1.4)
abline(v=c(nfit[1]-nfit[2], nfit[1]+nfit[2]), col=4, lwd=1.4, lty=2)
print(paste0("mean = ", round(nfit[1]), ",  StDev = ", round(nfit[2])))




##############################################################################
# plot width distributions by hydroBASIN 
##############################################################################

# set up output pdf:
pdf(pdfOut, width=15, height=9)
layoutTab = par(mfrow=c(4,5))
par(oma=c(3,3,0,0), mar=c(0,0,0,0))
layout(layoutTab)


# read in hBASIN shapefile dbf:
hBASIN = foreign::read.dbf(sub('hybas_allMain', 'hybas_allMainCopy', hydroBASINpath))

m = 1
i = 1
j = 1

for (i in 1:length(fP)){ # Amazon: i = 5687 # Miss: i=6829
  
  for (j in 1:nRun){
  
    # read in and filter GRWL data: ####
    if (j == 1){
      # entire GRWL database:
      # fP = paste0(wd, 'output/figs/allGRWL.pdf')
      csv_raw = read.csv(fP[i], header=T)
      
      # remove data with widths<100m, elev<1m, lakes, canals:
      keep = csv_raw$width_m>wMin & 
        csv_raw$elev_m>minElev & 
        csv_raw$lakeFlag!=1 & 
        csv_raw$lakeFlag!=3 
    } 
    
    csv = csv_raw
    
    
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
    
    # calculate the average number of channels: 
    nChan = mean(csv$nchannels[keep])
    
    
    # MLE fit: ####
    
    # set std=T to calc. std dev of fit using MLE optimization,
    # which is slow and produces very small std dev vals:
    pFit = pareto.mle(x=A, std=F) 
  
    # calculate GOF statistics with X2 and KS test:
    breaks = seq(Amin, Amax+int, int)
    h = hist(A, breaks, plot=F)
    jA = jitter(A) # to prevent ties in the following GOF tests
    pGOF = suppressWarnings(GOF(pFit, h, jA))
  
    
    # total surface area estimate: ####
    
    # calculate total surface area of rivers and streams by extendending
    # area abundance to smaller rivers:
    #pSA = extrapSA_calculator(pFit, fOaMeans, Amin, Amax, sumA)
    #gpSA = extrapSA_calculator(gpFit, fOaMeans, Amin, Amax, sumA)
    
    
    # calculate total surface area of rivers and streams by extendending
    # area abundance to smaller rivers, and calculate confidence intervals
    # using monteCarlo simulations:
    pMcarlo = mCarloUncertProp(pFit, fOaMean, fOaSD, Amin, Amax, N=2e4, sumA)
    gpMcarlo = mCarloUncertProp(gpFit, fOaMean, fOaSD, Amin, Amax, N=2e4, sumA)
    
    # calculate the % of land surface occupied by rivers and streams: 
    basinArea = hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]] #132773914
    obs_prcnt = 100*sumA*reducer*1e-6/basinArea
    pSA_prcnt = 100*pMcarlo$meanSimSA/basinArea
    pSA_sd_prct = 100*pMcarlo$sdSimSA/basinArea
    gpSA_prcnt = 100*gpMcarlo$meanSimSA/basinArea
    gpSA_sd_prct = 100*gpMcarlo$sdSimSA/basinArea
    
    # Plot histogram and Pareto fits (Fig. 2b):
    histPlot(fOaMeans, Amax, pFit, gpFit, sumA, 
             breaks, h, Alen, int, pMcarlo, gpMcarlo, 
             obs_prcnt, pSA_prcnt, gpSA_prcnt, pGOF, m)
    
    
    ####################################################
    # add data to statsTable:
    gTab[m, ] = as.vector(c(fN[i], 
                            round(basinArea),
                            length(which(keep)), 
                            nChan,
                            round(Amax),
                            round(sum(reducer*mL*1e-6*
                                        csv$width_m[csv$lakeFlag!=1 & csv$lakeFlag!=3])), 
                            round(sumA*reducer*1e-6),
                            # global fits:
                            round(gpMcarlo$meanSimSA), 
                            round(gpMcarlo$sdSimSA), 
                            round(gpSA_prcnt, 4), 
                            round(gpSA_sd_prct, 4), 
                            round(gpFit$xm*reducer, 2), 
                            round(gpFit$alpha, 3),
                            round(gpFit$stdev, 3),
                            # basin fits:
                            round(pMcarlo$meanSimSA), 
                            round(pMcarlo$sdSimSA), 
                            round(pSA_prcnt, 4), 
                            round(pSA_sd_prct, 4), 
                            round(pFit$xm*reducer, 2), 
                            round(pFit$alpha, 3),
                            round(pFit$stdev, 3),
                            round(pGOF$X2, 2), 
                            round(pGOF$X2_p, 3), 
                            round(pGOF$KS_D, 2), 
                            round(pGOF$KS_p, 3)
                          ))
    
    
    m = m + 1
  
  } # end Monte Carlo simulation
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)






trainingBasinRA = sum(as.numeric(gTab$gSA_km2[order(as.numeric(gTab$nObs), decreasing=T)][1:20]))
trainingBasinRA = sum(as.numeric(gTab$pSA_km2[order(as.numeric(gTab$nObs), decreasing=T)][1:20]))

print(paste("training basins RA:", round(100*trainingBasinRA/globalLandArea,2), "%"))

print(paste("observed RA in extrapolated basins:", 
            sum(as.numeric(gTab$obSA_km2)), "km2"))

totKm2 = c(
  sum(as.numeric(gTab$pSA_km2))-sum(as.numeric(gTab$pSA_sd_km2)),
  sum(as.numeric(gTab$pSA_km2)),
  sum(as.numeric(gTab$pSA_km2))+sum(as.numeric(gTab$pSA_sd_km2)))
print("extrapolated RA:")
print(totKm2)
print(100*totKm2/globalLandArea)


paste(round(100*sum(as.numeric(gTab$obSA_gt90m_km2))/sum(as.numeric(gTab$gSA_km2)), 1), "% observed area used in extrap")
paste(round(100*sum(as.numeric(gTab$obSA_km2))/sum(as.numeric(gTab$gSA_km2)), 1), "% observed area")


# write out statistics table:
#write.csv(gTab, gTabP, row.names=F)


## elev>5, alphaN=37; minN=5000 --> 600095, -464410, +856176
## elev>1, alphaN=38; minN=100,000 --> 356679.6, -465322, +678286.9




# Attach gtab attributes to hydroBASIN shapefile:
##############################################################################
io = 0; if (io == 1){
  #print(gTab)
  
 # gTab = read.csv(gTabP, header=T)
  if (length(grep("dbf", names(gTab)))>0){ gTab = gTab$dbf }
  
  # read in hBASIN dbf and attached gtab attributes to it:
  hBASIN = read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
  if (length(grep("dbf", names(hBASIN)))>0){ hBASIN = hBASIN$dbf }
  
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







# Import real fits from large river networks:
####################################################
io = 0; if (io == 1){
  #print(gTab)
  
 # gTab = read.csv(gTabP, header=T)
  if (length(grep("dbf", names(gTab)))>0){ gTab = gTab$dbf }
  if (length(grep("Amax", names(gTab)))>0){ gTab = gTab[, -which(names(gTab)=="Amax")] }
  
  bPath = 'H:/2017_02_08_GRWL_Nature_Manuscript/figs/fig2/hydroBasin_SAfitStatistics_minW90_NobsGt100000.csv'
  bgTab = read.csv(bPath, header=T)
  if (length(grep("dbf", names(bgTab)))>0){ bgTab = gTab$dbf }
  
  j = match(bgTab$hBASIN_code, gTab$hBASIN_code)
  gTab = gTab[-j,]
  gTab = rbind(gTab, bgTab)
  
  # read in hBASIN dbf and attached gtab attributes to it:
  hBASIN = read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
  if (length(grep("dbf", names(hBASIN)))>0){ hBASIN = hBASIN$dbf }
  
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
