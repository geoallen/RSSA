##############################################################################
# RSSA_calc.R
##############################################################################
# George H. Allen

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
tabDirPath = paste0(wd, 'output/tabs/')
nGRWLperBasinOutPath = paste0(tabDirPath, "nGRWLperBasinTab/nGRWLperBas.csv")
mnTabP = paste0(wd, 'output/figs/fig3_SAfits_20basins_MCmn.csv')
sdTabP = paste0(wd, 'output/figs/fig3_SAfits_20basins_MCsd.csv')
hydroBASINpath = paste0(wd, 'input/basin_shapefiles/hydroBASINs/hybas_allMain.dbf')
valPath = paste0(wd, 'input/validation/Database_S1_GRWL_validation_data.csv')


# get file paths & codes: 
fP = list.files(outPath, 'csv', full.names=T)
fN = list.files(outPath, 'csv')
fN = sub('.csv', '', fN)
hBASIN_code = sub("GRWL_", '', fN)
ensembleTabPath = sub(outPath, paste0(tabDirPath, 'ensembleOuputTabs/'), fP)
hTabPath = sub(outPath, paste0(tabDirPath, 'ensembleHistTabs/'), fP)

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

# number of Monte Carlo simulation ensemble runs to estimate error:
nRun = 20

# create table to store outputs of each Monte Carlo enseble run:
ensembleTabNames = c("hBASIN_code", 
              "basinA_km2", 
              "nObs",
              "nObsA",
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
ensembleTab = data.frame(array(NA, c(nRun, length(ensembleTabNames)))) 
names(ensembleTab) = ensembleTabNames

# create tabs to store ensemble mean and stdev output data:
mnTab = data.frame(array(NA, c(1, length(ensembleTabNames))))
names(mnTab) = ensembleTabNames
sdTab = mnTab

# global MLE mean fit: 
gpFit = list(xm=Amin, 
             alpha=0.90345, # 20 basin with most obs, min elev = 0
             stdev=0.06285026) # 20 basins with most obs, min elev = 0

# total surface area of Earth's non glaciated land surface:
globalLandArea = 132773914

##############################################################################
# functions
##############################################################################

# list N GRWL measurements in each hBasin. This list is used to efficiently
# sort through the classes of basins used in RSSA extrapolation:
nGRWLperBasListGenerator <- function(fP, tabDirPath, nGRWLperBasinOutPath){
  # generate list of GRWL measurements in each hBasin:
  NfP = length(fP)
  fP_ind = 1:NfP
  nGRWLperBas = rep(NA, NfP)
  for (i in 1:NfP){ 
    print(paste("reading in", i, "of", NfP))
    # read in GRWL by basin and remove data with widths<100m, elev<1m, etc:
    csv_raw = read.csv(fP[i], header=T)
    keep = csv_raw$width_m>wMin & 
      csv_raw$elev_m>minElev & 
      csv_raw$lakeFlag!=1 & 
      csv_raw$lakeFlag!=3 
    nGRWLperBas[i] = length(which(keep))
  }
  # create table and sort rows by n GRWL obs:
  nGRWLperBasTab = data.frame(fP_ind, fN, fP, nGRWLperBas)
  nGRWLperBasTab = nGRWLperBasTab[order(nGRWLperBasTab$nGRWLperBas, decreasing=T),]
  # write table out:
  if (!dir.exists(paste0(tabDirPath, "nGRWLperBasinTab/"))){
    dir.create(paste0(tabDirPath, "nGRWLperBasinTab/"))
  }
  write.csv(nGRWLperBasTab, nGRWLperBasinOutPath, row.names=F)
}

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
       breaks[-1], h, 
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


# round columns of mnTab:
tabRounder <- function(tab){
  
  tab$hBASIN_code = tab$hBASIN_code
  tab$basinA_km2 = tab$basinA_km2
  tab$nObs = tab$nObs
  tab$nObsA = round(Alen)
  tab$nChan_mean = tab$nChan_mean
  tab$Amax = round(tab$Amax, 1)
  tab$obSA_km2 = round(tab$obSA_km2, 1)
  tab$obSA_gt90m_km2 = round(tab$obSA_gt90m_km2, 1)
  # global fits:
  tab$gSA_km2 = round(tab$gSA_km2, 1)
  tab$gSA_sd_km2 = round(tab$gSA_sd_km2, 1)
  tab$gSA_pcnt = round(tab$gSA_pcnt, 4)
  tab$gSA_sd_pcnt = round(tab$gSA_sd_pcnt, 4)
  tab$gXmin = round(tab$gXmin, 1)
  tab$gAlpha = round(tab$gAlpha, 3)
  tab$gAlpha_sd = round(tab$gAlpha_sd, 3)
  # basin fits:
  tab$pSA_km2 = round(tab$pSA_km2, 1)
  tab$pSA_sd_km2 = round(tab$pSA_sd_km2, 1)
  tab$pSA_pcnt = round(tab$pSA_pcnt, 4)
  tab$pSA_sd_pcnt = round(tab$pSA_sd_pcnt, 4)
  tab$pXmin = round(tab$pXmin, 2)
  tab$pAlpha = round(tab$pAlpha, 3)
  tab$pAlpha_sd = round(tab$pAlpha_sd, 3)
  tab$pX2_stat = round(tab$pX2_stat, 2)
  tab$pX2_p = round(tab$pX2_p, 3)
  tab$pKS_D = round(tab$pKS_D, 2)
  tab$pKS_p = round(tab$pKS_p, 3)

  return(tab)
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
lines(lineSeq, fn, col=4)
abline(v=nfit[1], col=4, lwd=1.4)
abline(v=c(nfit[1]-nfit[2], nfit[1]+nfit[2]), col=4, lwd=1.4, lty=2)
print(paste0("GRWL error Monte Carlo  mean = ", round(nfit[1]), 
             ",  StDev = ", round(nfit[2])))



##############################################################################
# plot RSSA distributions in each hydroBASIN 
##############################################################################

# if list file does not exist, create a list N GRWL obs. in each hBasin: 
if (!file.exists(nGRWLperBasinOutPath)){
  nGRWLperBasListGenerator(fP, tabDirPath, nGRWLperBasinOutPath)
}
# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)

# read in hBASIN shapefile dbf:
hBASIN = foreign::read.dbf(sub('hybas_allMain', 'hybas_allMainCopy', hydroBASINpath))


# Class A basins:

# Class A basins contain >250,000 river measurements:
classA_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>250000])
classA_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>250000]
nClassA = length(classA_fP)
print(paste("Class A Basins:", paste(classA_fN, collapse=" ")))

# in calss A basins, We use a Monte Carlo simulation to characterize 
# RSSA error due to the uncertainty of GRWL measurements. For each ensemble 
# run, we fit a Pareto distribution of RSSA measurements using maximum
# likelihood estimation (MLE). To fit the distribution, we set the 
# Pareto scale parameter, xm, equal to the lowermost observed RSSA value  
# and solve for the Pareto shape parameter. This lowermost observed RSSA, 
# Omin=3259 m2, is the product of the minimum observed river width threshold 
# (90 m) and the mean distance between river centerline pixels in 30-m 
# resolution Landsat imagery (36.21 m). 
# Quantify the distribution of these statistical fits

# get start time:
old <- Sys.time() 

# set up output pdf:
pdf(pdfOut, width=15, height=9)
par(mfrow=c(4,5))
par(oma=c(3,3,0,0), mar=c(0,0,0,0))


m = 1
for (i in 1:nClassA){ 
  
  # read in and filter GRWL data: ####
  csv_raw = read.csv(classA_fP[i], header=T)
  
  # remove data with widths<100m, elev<1m, lakes, canals:
  keep = csv_raw$width_m>wMin & 
    csv_raw$elev_m>minElev & 
    csv_raw$lakeFlag!=1 & 
    csv_raw$lakeFlag!=3 
  
  N = length(which(keep))
  
  print(paste('i =',i, ' of ', nClassA,', fN = ', classA_fN[i]))
  
  # calculate the average number of channels: 
  nChan = mean(csv_raw$nchannels[keep])
  # calculate distance between each adjacent GRWL centerline cell:
  d = distCalc(csv_raw[keep,])
  
  
  # Monte Carlo error propogation: ####
  for (j in 1:nRun){
    
    print(paste("Run:", j))
    # reset csv to original table:
    csv = csv_raw[keep, ]
    
    # generate monte carlo simulation pertubations:
    w_perturb = rnorm(N, nfit[1], nfit[2])
    
    # calc width, distance, river area, sum river area, maximum area value: 
    w_raw = csv$width_m + w_perturb
    w = w_raw/reducer
    A_raw = w*d
    A = A_raw[A_raw>Amin]
    sumA = sum(A)
    Amax = max(A)
    Alen = length(A)
    

    
    # MLE fit: ####
    
    # set std=T to calc. std dev of fit using MLE optimization,
    # which is slow and produces negligible std dev vals:
    pFit = pareto.mle(x=A, std=F) 
  
    # calculate GOF statistics with X2 and KS test:
    breaks = seq(Amin, Amax+int, int)
    h = hist(A, breaks, plot=F)
    jA = jitter(A) # to prevent ties in the following GOF tests
    pGOF = suppressWarnings(GOF(pFit, h, jA))
  
    # add histogram information to table to be in plot:
    if (j == 1){
      hTab = as.data.frame(array(NA, c(nRun, length(h$mids))))
    }
    hCounts = c(h$counts, rep(0, 100))
    hTab[j,] = hCounts[1:ncol(hTab)]
    
    
    # RSSA estimate with added uncertainty: ####
    
    # calculate total surface area of rivers and streams by extending the
    # RSSA abundance to smaller rivers, and calculate confidence intervals
    # using a Monte Carlo simulation:
    pMcarlo = mCarloUncertProp(pFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    # global 
    gpMcarlo = mCarloUncertProp(gpFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    
    # calculate the % of land surface occupied by rivers and streams: 
    basinArea = hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]] #132773914
    obs_prcnt = 100*sumA*reducer*1e-6/basinArea
    pSA_prcnt = 100*pMcarlo$meanSimSA/basinArea
    pSA_sd_prct = 100*pMcarlo$sdSimSA/basinArea
    gpSA_prcnt = 100*gpMcarlo$meanSimSA/basinArea
    gpSA_sd_prct = 100*gpMcarlo$sdSimSA/basinArea
    
    # fill in table with outputs of ensemble run: ####
    ensembleTab[j,] =
      as.vector(c(
                fN[i], 
                basinArea,
                length(which(keep)), 
                Alen,
                nChan,
                Amax,
                sum(reducer*mL*1e-6*
                      csv$width_m[csv$lakeFlag!=1 & csv$lakeFlag!=3]), 
                sumA*reducer*1e-6,
                # global fits:
                gpMcarlo$meanSimSA, 
                gpMcarlo$sdSimSA, 
                gpSA_prcnt, 
                gpSA_sd_prct, 
                gpFit$xm*reducer, 
                gpFit$alpha,
                gpFit$stdev,
                # basin fits:
                pMcarlo$meanSimSA, 
                pMcarlo$sdSimSA, 
                pSA_prcnt, 
                pSA_sd_prct, 
                pFit$xm*reducer, 
                pFit$alpha,
                pFit$stdev,
                pGOF$X2, 
                pGOF$X2_p, 
                pGOF$KS_D, 
                pGOF$KS_p
                ))
    
  
  } # end Monte Carlo simulation
  
  # write out ensemble histogram table:
  write.csv(hTab, hTabPath[i], row.names=F)
  # write out ensemble output table:
  write.csv(ensembleTab, ensembleTabPath[i], row.names=F)
  
  
  # take column means and add them to the mnTab output table:
  mat = matrix(as.numeric(unlist(ensembleTab)), nrow=nRun)
  mnTab[m, ] = colMeans(mat, na.rm=T)
  sdTab[m, ] = apply(mat, 2, sd)
  
  # Plot histogram and Pareto fits (Fig. 2b):
  histPlot(fOaMeans, 
           mnTab$Amax[m], 
           pFit=list(xm=mnTab$pXmin[m]/reducer, 
                     alpha=mnTab$pAlpha[m],
                     stdev=sdTab$pAlpha[m]),
           gpFit, 
           mnTab$obSA_gt90m_km2[m]*reducer^2, #sdTab$obSA_gt90m_km2[m] 
           breaks, 
           round(as.vector(colMeans(hTab, na.rm=T))),
           mnTab$nObsA[m], #sdTab$nObsA[m]
           int, 
           pMcarlo = list(meanSimSA=mnTab$pSA_km2[m],
                          sdSimSA=sdTab$pSA_km2[m]),
           gpMcarlo, 
           100*mnTab$obSA_gt90m_km2[m]/basinArea, # 100*sdTab$obSA_gt90m_km2[m]/basinArea  
           mnTab$pSA_pcnt[m], #sdTab$pSA_pcn[m]t
           mnTab$gSA_pcnt[m], #sdTab$gSA_pcnt[m], 
           pGOF = list(X2=mnTab$pX2_stat[m], 
                       X2_p=mnTab$pX2_p[m], 
                       KS_D=mnTab$pKS_D[m], 
                       KS_p=mnTab$pKS_p[m]), 
           # sdTab$pX2_stat[m], sdTab$pX2_p[m], sdTab$pKS_D[m], sdTab$pKS_p[m], 
           m)
  
  m = m + 1
}

# round columns of mnTab:
mnTab = tabRounder(mnTab)
sdTab = tabRounder(sdTab)

# write out mean and stdev output tables:
#write.csv(mnTab, mnTabP, row.names=F)
#write.csv(sdTab, sdTabP, row.names=F)

#mnTab = read.csv(mnTabP, header=T)
#sdTab = read.csv(sdTabP, header=T)

dev.off()
cmd = paste('open', pdfOut)
system(cmd)

new = Sys.time() - old
print(new) 

# add up the total observed & extrapolated river surface area in Class A basins:
classA_obs_gt90m = sum(mnTab$obSA_gt90m_km2)
classA_obs = sum(mnTab$obSA_km2)
classA_RSSA = sum(mnTab$pSA_km2)
classA_BasinA = sum(mnTab$basinA_km2)
print(paste("Class A basins observed %RSSA > 90m:", round(100*classA_obs_gt90m/classA_BasinA,2), "%"))
print(paste("Class A basins observed %RSSA:", round(100*classA_obs/classA_BasinA,2), "%"))
print(paste("Class A basins extrapolated %RSSA:", round(100*classA_RSSA/classA_BasinA,2), "%"))

# add uncertainty:
totKm2 = c(
  sum(as.numeric(mnTab$pSA_km2))-sum(as.numeric(sdTab$pSA_km2)),
  sum(as.numeric(mnTab$pSA_km2)),
  sum(as.numeric(mnTab$pSA_km2))+sum(as.numeric(sdTab$pSA_km2)))

print("extrapolated RSSA:")
print(totKm2)
print(100*totKm2/globalLandArea)



# plot Fig. S6: Pareto fits in each Class A basin


# set up output pdf:
pdf(pdfOut, width=15, height=9)
par(mfrow=c(4,5))
par(oma=c(3,3,0,0), mar=c(0,0,0,0))

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
       breaks[-1], h, 
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

dev.off()
cmd = paste('open', pdfOut)
system(cmd)


## elev>5, alphaN=37; minN=5000 --> 600095, -464410, +856176
## elev>1, alphaN=38; minN=100,000 --> 356679.6, -465322, +678286.9




# Attach mnTab attributes to hydroBASIN shapefile:
##############################################################################
io = 0; if (io == 1){
  #print(mnTab)
  
 # mnTab = read.csv(mnTabP, header=T)
  if (length(grep("dbf", names(mnTab)))>0){ mnTab = mnTab$dbf }
  
  # read in hBASIN dbf and attached mnTab attributes to it:
  hBASIN = read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
  if (length(grep("dbf", names(hBASIN)))>0){ hBASIN = hBASIN$dbf }
  
  bindTab = data.frame(array(NA, c(nrow(hBASIN), ncol(mnTab))))
  names(bindTab) = names(mnTab)
  matchInd = match(mnTab$hBASIN_code, hBASIN$MAIN_BAS)
  
  gMat = matrix(as.matrix(mnTab), ncol = ncol(mnTab), dimnames = NULL)
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
  #print(mnTab)
  
 # mnTab = read.csv(mnTabP, header=T)
  if (length(grep("dbf", names(mnTab)))>0){ mnTab = mnTab$dbf }
  if (length(grep("Amax", names(mnTab)))>0){ mnTab = mnTab[, -which(names(mnTab)=="Amax")] }
  
  bPath = 'H:/2017_02_08_GRWL_Nature_Manuscript/figs/fig2/hydroBasin_SAfitStatistics_minW90_NobsGt100000.csv'
  bgTab = read.csv(bPath, header=T)
  if (length(grep("dbf", names(bgTab)))>0){ bgTab = mnTab$dbf }
  
  j = match(bgTab$hBASIN_code, mnTab$hBASIN_code)
  mnTab = mnTab[-j,]
  mnTab = rbind(mnTab, bgTab)
  
  # read in hBASIN dbf and attached mnTab attributes to it:
  hBASIN = read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
  if (length(grep("dbf", names(hBASIN)))>0){ hBASIN = hBASIN$dbf }
  
  bindTab = data.frame(array(NA, c(nrow(hBASIN), ncol(mnTab))))
  names(bindTab) = names(mnTab)
  matchInd = match(mnTab$hBASIN_code, hBASIN$MAIN_BAS)
  
  gMat = matrix(as.matrix(mnTab), ncol = ncol(mnTab), dimnames = NULL)
  for (i in 1:length(matchInd)){
    bindTab[matchInd[i], ] = gMat[i,]
  }
  
  # format table for arcmap shapefile:
  newhBASIN = as.data.frame(cbind(hBASIN, bindTab))
  newhBASIN[is.na(newhBASIN)] = 0
  newhBASIN = data.matrix(newhBASIN)
  
  write.dbf(newhBASIN, hydroBASINpath)
  
}

















#### CLASS B BASIN CODE:

# get start time:
old <- Sys.time() 

# set up output pdf:
pdf(pdfOut, width=15, height=9)
par(mfrow=c(4,5))
par(oma=c(3,3,0,0), mar=c(0,0,0,0))


m = 1
for (i in 1:length(fP)){ 
  
  # read in and filter GRWL data: ####
  csv_raw = read.csv(fP[i], header=T)
  
  # remove data with widths<100m, elev<1m, lakes, canals:
  keep = csv_raw$width_m>wMin & 
    csv_raw$elev_m>minElev & 
    csv_raw$lakeFlag!=1 & 
    csv_raw$lakeFlag!=3 
  
  N = length(which(keep))
  
  # skip basin if there are less than the minimum number of observations:
  if (N < minNobs){ next }
  
  print(paste('i =',i, ' of ', length(fP),', fN = ', fN[i]))
  
  # calculate the average number of channels: 
  nChan = mean(csv_raw$nchannels[keep])
  # calculate distance between each adjacent GRWL centerline cell:
  d = distCalc(csv_raw[keep,])
  
  
  # Monte Carlo error propogation simulation ####
  for (j in 1:nRun){
    
    print(paste("Run:", j))
    # reset csv to original table:
    csv = csv_raw[keep, ]
    
    # generate monte carlo simulation pertubations:
    w_perturb = rnorm(N, nfit[1], nfit[2])
    
    # calc width, distance, river area, sum river area, maximum area value: 
    w_raw = csv$width_m + w_perturb
    w = w_raw/reducer
    A_raw = w*d
    A = A_raw[A_raw>Amin]
    sumA = sum(A)
    Amax = max(A)
    Alen = length(A)
    
    
    
    # MLE fit: ####
    
    # set std=T to calc. std dev of fit using MLE optimization,
    # which is slow and produces negligible std dev vals:
    pFit = pareto.mle(x=A, std=F) 
    
    # calculate GOF statistics with X2 and KS test:
    breaks = seq(Amin, Amax+int, int)
    h = hist(A, breaks, plot=F)
    jA = jitter(A) # to prevent ties in the following GOF tests
    pGOF = suppressWarnings(GOF(pFit, h, jA))
    
    # add histogram information to table to be in plot:
    if (j == 1){
      hTab = as.data.frame(array(NA, c(nRun, length(h$mids))))
    }
    hCounts = c(h$counts, rep(0, 100))
    hTab[j,] = hCounts[1:ncol(hTab)]
    
    
    # RSSA estimate with added uncertainty: ####
    
    # calculate total surface area of rivers and streams by extending the
    # RSSA abundance to smaller rivers, and calculate confidence intervals
    # using a Monte Carlo simulation:
    pMcarlo = mCarloUncertProp(pFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    # global 
    gpMcarlo = mCarloUncertProp(gpFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    
    # calculate the % of basin occupied by rivers and streams: 
    basinArea = hBASIN$area_km2[hBASIN$MAIN_BAS == fN[i]]
    obs_prcnt = 100*sumA*reducer*1e-6/basinArea
    pSA_prcnt = 100*pMcarlo$meanSimSA/basinArea
    pSA_sd_prct = 100*pMcarlo$sdSimSA/basinArea
    gpSA_prcnt = 100*gpMcarlo$meanSimSA/basinArea
    gpSA_sd_prct = 100*gpMcarlo$sdSimSA/basinArea
    
    # fill in table with outputs of ensemble run: ####
    ensembleTab[j,] =
      as.vector(c(
        fN[i], 
        basinArea,
        length(which(keep)), 
        Alen,
        nChan,
        Amax,
        sum(reducer*mL*1e-6*
              csv$width_m[csv$lakeFlag!=1 & csv$lakeFlag!=3]), 
        sumA*reducer*1e-6,
        # global fits:
        gpMcarlo$meanSimSA, 
        gpMcarlo$sdSimSA, 
        gpSA_prcnt, 
        gpSA_sd_prct, 
        gpFit$xm*reducer, 
        gpFit$alpha,
        gpFit$stdev,
        # basin fits:
        pMcarlo$meanSimSA, 
        pMcarlo$sdSimSA, 
        pSA_prcnt, 
        pSA_sd_prct, 
        pFit$xm*reducer, 
        pFit$alpha,
        pFit$stdev,
        pGOF$X2, 
        pGOF$X2_p, 
        pGOF$KS_D, 
        pGOF$KS_p
      ))
    
    
  } # end Monte Carlo simulation
  
  # write out ensemble histogram table:
  write.csv(hTab, hTabPath[i], row.names=F)
  # write out ensemble output table:
  write.csv(ensembleTab, ensembleTabPath[i], row.names=F)
  
  
  # take column means and add them to the mnTab output table:
  mat = matrix(as.numeric(unlist(ensembleTab)), nrow=nRun)
  mnTab[m, ] = colMeans(mat, na.rm=T)
  sdTab[m, ] = apply(mat, 2, sd)
  
  # Plot histogram and Pareto fits (Fig. 2b):
  histPlot(fOaMeans, 
           mnTab$Amax[m], 
           pFit=list(xm=mnTab$pXmin[m]/reducer, 
                     alpha=mnTab$pAlpha[m],
                     stdev=sdTab$pAlpha[m]),
           gpFit, 
           mnTab$obSA_gt90m_km2[m]*reducer^2, #sdTab$obSA_gt90m_km2[m] 
           breaks, 
           round(as.vector(colMeans(hTab, na.rm=T))),
           mnTab$nObsA[m], #sdTab$nObsA[m]
           int, 
           pMcarlo = list(meanSimSA=mnTab$pSA_km2[m],
                          sdSimSA=sdTab$pSA_km2[m]),
           gpMcarlo, 
           100*mnTab$obSA_gt90m_km2[m]/basinArea, # 100*sdTab$obSA_gt90m_km2[m]/basinArea  
           mnTab$pSA_pcnt[m], #sdTab$pSA_pcn[m]t
           mnTab$gSA_pcnt[m], #sdTab$gSA_pcnt[m], 
           pGOF = list(X2=mnTab$pX2_stat[m], 
                       X2_p=mnTab$pX2_p[m], 
                       KS_D=mnTab$pKS_D[m], 
                       KS_p=mnTab$pKS_p[m]), 
           # sdTab$pX2_stat[m], sdTab$pX2_p[m], sdTab$pKS_D[m], sdTab$pKS_p[m], 
           m)
  
  m = m + 1
}

# round columns of mnTab:
mnTab = tabRounder(mnTab)
sdTab = tabRounder(sdTab)

# write out mean and stdev output tables:
write.csv(mnTab, mnTabP, row.names=F)
write.csv(sdTab, sdTabP, row.names=F)

#mnTab = read.csv(mnTabP, header=T)
#sdTab = read.csv(sdTabP, header=T)

dev.off()
cmd = paste('open', pdfOut)
system(cmd)


new = Sys.time() - old
print(new) 








