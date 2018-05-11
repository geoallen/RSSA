##############################################################################
# RSSA_calc.R
##############################################################################
# George H. Allen
#
# Description: 
# This code calculate River and Stream Surface Area (RSSA) for all non-
# glaciated regions globally using the river width and length measurements
# from the GRWL database. 
# 
# Classify each HydroBASIN (hBASIN) into one of three categories:
# Class A: basins that contain >250,000 measurements (N=20)
# Class B: basins that contain 10,000 to 250,000 measurements
# Class C: basins that contain ≤10,000 measurements 
# In Class A basins, we estimated the total RSSA by extending a fitted Pareto 
# frequency distribution down to the median first-order wetted stream width of
# 32±8 cm (22) (Fig. 3C). In Class B basins (N=294), which contain insufficient
# GRWL data to exhibit a well-developed fractal RSSA distribution, we use the 
# average Pareto shape parameter established in Class A basins (fig. S2) to 
# extend the RSSA distribution to first order streams (Fig. 3D). Class C 
# basins contain very little GRWL data, so we developed an empirical power-law 
# relationship between climate aridity (23), basin area (21), and percent 
# basin occupied by RSSA.
# 
# Inputs: 
# 1. GRWL data Subsetted by HydroBASIN in ArcGIS (used Spatial Join) from 
# Allen and Pavelsky (2018)
# 2. Table listing Class A basin codes and their names
# 3. HydroBASINS: Lehner et al. (2013) https://doi.org/10.1002/hyp.9740
# 4. In situ river width validation data from Allen et al. (2018)
# 5. Aridity Index from Zomer et al. (2008): 
# https://doi.org/10.1016/j.agee.2008.01.014
#
# Outputs:
# tables and figures presented in Allen et al. (2018)
#
#
##############################################################################
# Load R packages
##############################################################################
if (!"foreign" %in% rownames(installed.packages())){
  install.packages("foreign")}; require(foreign)
if (!"MASS" %in% rownames(installed.packages())){
  install.packages("MASS")}; require(MASS)


##############################################################################
# Hard-coded variables
##############################################################################

# specify root paths:
wd = "/Users/geoallen/Documents/research/2018_01_08_Science_submission_GRWL/git/RSSA/"

# Input file paths:
# GRWL data Subsetted by HydroBASIN in ArcGIS (used Spatial Join) from 
# Allen and Pavelsky (2018):
GRWLpath = paste0(wd, 'input/GRWL/GRWL_by_hydroBASIN/')
# HydroBASINS: Lehner et al. (2013) https://doi.org/10.1002/hyp.9740
hydroBASINpath = paste0(wd, 'input/basin_shapefiles/hydroBASINs/hybas_allMain.dbf')
# In situ river width validation data from Allen et al. (2018)
valPath = paste0(wd, 'input/validation/Database_S1_GRWL_validation_data.csv')
# Table listing Class A basin codes and their names
classA_bNamesPath = paste0(wd, 'input/misc/classA_bNames.csv')
# Aridity Index from Zomer et al. (2008): 
# https://doi.org/10.1016/j.agee.2008.01.014
aridityPath = paste0(wd, 'input/aridity/aridityByBasin.dbf')

# Output file paths:
tabDirPath = paste0(wd, 'output/tabs/')
figDirPath = paste0(wd, 'output/figs/')
nGRWLperBasinOutPath = paste0(tabDirPath, "nGRWLperBasinTab/nGRWLperBas.csv")

# get file paths & codes: 
fP = list.files(GRWLpath, 'csv', full.names=T)
fN = list.files(GRWLpath, 'csv')
fN = sub('.csv', '', fN)
hBASIN_code = sub("GRWL_", '', fN)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs'), fP)

# define constants: 
reducer = 100 # convert area units to prevent overflow on max. liklihood. est.
minL = 30/reducer # pixel resoluton
mL = mean(c(30, 30*sqrt(2)))/reducer # mean length between GRWL centerline pixels
int = 100*mL  # plotting histogram binning interval
wMin = 90 # minimum river width to threshold
Amin = wMin*mL # minimum river area value to threshold
minNobs = 2.5e5  # minimum number of GRWL observations in a basin for analysis
minElev = 0 # include rivers above this elevation (m)
hMax = 10725 # maximum limit to Class A histogram graphs
figBreaks = seq(Amin, hMax+int, int)*reducer
nRun = 500 # number of Monte Carlo simulation ensemble runs
# set lower extrapolation bounds from Allen et al. (2018): 
# https://doi.org/10.1038/s41467-018-02991-w : 
fOaMean = 0.321*mL # mean of the median widths of 1st order streams
fOaSD = 0.077*mL # sd width of 1st order streams
fOaMeans = c(fOaMean-fOaSD, fOaMean, fOaMean+fOaSD)
globalLandArea = 132773914 # area of Earth's non-glaciated land surface
# Class A MLE mean fit: 
aFit = list(
  xm=Amin*reducer,
  alpha=0.9034167, # class A basins with most obs, min elev = 0
  stdev=0.06285026 # class A basins with most obs, min elev = 0
) 

# specify plotting parameters for Class C basin regression:
yRange = c(0,3)
xRange = c(0, 2.2e4)
NxSeq = 100
xSeq1 = seq(xRange[1], xRange[2], length.out=NxSeq)
xSeq2 = seq(min(x2), max(x2), length.out=NxSeq)
confInt = 0.68

# create large table to store outputs of each Monte Carlo enseble run:
ensembleTabNames = c("hBASIN_code", 
              "BA_km2", 
              "nObs",
              "bClass",
              "nChan_mn",
              "Amax",
              "obRSSA_km2",
              "obRSSAgt90_km2", 
              # global
              "aRSSA_km2",
              "aMCRSSA_km2",
              "aRSSA_pc",
              "gMCRSSA_pc",
              "gXmin",
              "gMCalpha", 
              "gMCfOw", 
              # basin
              "bRSSA_km2",
              "bMCRSSA_km2",
              "bRSSA_pc",
              "pMCRSSA_pc",
              "pXmin",
              "pMCalpha", 
              "pMCfOw", 
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


##############################################################################
# Functions
##############################################################################

# list N GRWL measurements in each hBasin. This list is used to efficiently
# sort through the classes of basins used in RSSA extrapolation:
nGRWLperBasListGenerator <- function(fP, fN, tabDirPath, wMin, minElev, nGRWLperBasinOutPath){
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
GOF = function(bFit, h, jA){
  
  # Pearson's Chi-squared test for count data:
  # low X2 and high p-value signifies a good fit:
  f = dpareto(h$mids, bFit$xm, bFit$alpha)
  X2test = chisq.test(h$counts, p=f, rescale.p=T, simulate.p.value=T, B=10)
  
  # Two sided One sample KS GOF test:
  # low D and high p-value signifies a good fit:
  KStest = ks.test(jA, "ppareto",  bFit$xm, bFit$alpha, alternative="two.sided")
  
  return(list(
    X2 = X2test$statistic,
    X2_p = X2test$p.value,
    KS_D = KStest$statistic[[1]],
    KS_p = KStest$p.value
  ))
}

extrapSA_calculator = function(bFit, fOaMeans, Amin, Amax, sumA){
  
  a = c(bFit$alpha-bFit$stdev, bFit$alpha, bFit$alpha+bFit$stdev)
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

# include uncertainty of extrapolation minimum and, for Class B basins, 
# uncertainty of Pareto fit to produce a Monte Carlo simulated RSSA estimate:
RSSAextrapolater <- function(bFit, fOaMean, fOaSD, Amin, Amax, N, sumA){
  
  # calculate total surface area:
  a = bFit$alpha
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
  
  # If Class B basin, use uncertainty of Pareto fits of Class A
  # basins to estimate RSSA using a Monte Carlo simulation. This isn't
  # needed in Class A basins becaues the Pareto fit is fitted directly
  # on the width data:
  if (is.na(bFit$stdev)){
    a = bFit$alpha
  }else{
    a = rnorm(N, bFit$alpha, bFit$stdev)
  }
  
  # Estimate RSSA uncertainty associated with uncertain lower extrap. bound
  # using a Monte Carlo simulation:
  xm = fOaMean
  x1 = rnorm(N, fOaMean, fOaSD)
  
  obsIntegral = -(a*x3*(xm/x3)^a)/(a-1) + (a*x2*(xm/x2)^a)/(a-1)
  extrapIntegral = -(a*x2*(xm/x2)^a)/(a-1) + (a*x1*(xm/x1)^a)/(a-1)
  extrap2ObsRatio = extrapIntegral/obsIntegral
  # add observed to estimated area:
  simSA = mean(sumA*extrap2ObsRatio+sumA)
  
  return(list(
    meanRSSA =  meanSA*reducer*1e-6,
    MCRSSA = mean(simSA, na.rm=T)*reducer*1e-6,
    MCAlpha = a,
    MCfOw = x1
  ))
}

# round columns of mnTab:
tabRounder <- function(tab){
  
  tab$hBASIN_code = tab$hBASIN_code
  tab$BA_km2 = tab$BA_km2
  tab$nObs = tab$nObs
  tab$bClass = tab$bClass
  tab$nChan_mn = tab$nChan_mn
  tab$Amax = round(tab$Amax, 1)
  tab$obRSSA_km2 = round(tab$obRSSA_km2, 1)
  tab$obRSSAgt90_km2 = round(tab$obRSSAgt90_km2, 1)
  # Class B RSSA:
  tab$aRSSA_km2 = round(tab$aRSSA_km2, 1)
  tab$aMCRSSA_km2 = round(tab$aMCRSSA_km2, 1)
  tab$aRSSA_pc = round(tab$aRSSA_pc, 4)
  tab$gMCRSSA_pc = round(tab$gMCRSSA_pc, 4)
  tab$gXmin = round(tab$gXmin, 1)
  tab$gMCalpha = round(tab$gMCalpha, 3)
  tab$gMCfOw = round(tab$gMCfOw, 3)
  # Class A RSSA:
  tab$bRSSA_km2 = round(tab$bRSSA_km2, 1)
  tab$bMCRSSA_km2 = round(tab$bMCRSSA_km2, 1)
  tab$bRSSA_pc = round(tab$bRSSA_pc, 4)
  tab$pMCRSSA_pc = round(tab$pMCRSSA_pc, 4)
  tab$pXmin = round(tab$pXmin, 2)
  tab$pMCalpha = round(tab$pMCalpha, 3)
  tab$pMCfOw = round(tab$pMCfOw, 3)
  tab$pX2_stat = round(tab$pX2_stat, 2)
  tab$pX2_p = round(tab$pX2_p, 3)
  tab$pKS_D = round(tab$pKS_D, 2)
  tab$pKS_p = round(tab$pKS_p, 3)

  return(tab)
}

# Large function that fits Pareto function to GRWL data. 
# Runs in two modes:
# First is non Monte Carlo simulation that simply estimates mean parameters
# using non pertubated GRWL width data.
# Second is Monte Carlo simulation that estimates uncertainty. 
classAB_MC <- function(nEnsRun, bClass){
  
  # if calculating mean Pareto alpha, just do one run and don't perturb data:
  nRun = nEnsRun
  
  # determine basin class:
  if (bClass==1){ classAB_fN=classA_fN }
  if (bClass==2){ classAB_fN=classB_fN }
  
  # Monte Carlo error propogation:
  for (j in 1:nRun){
    
    print(paste("Run:", j, "of", nRun))
    
    # reset perturbed GRWL table to original GRWL data: 
    csv = csv_raw
    
    # generate monte carlo simulation width pertubations:
    if (nRun > 1){
      w_perturb = rnorm(N_raw, nfit[1], nfit[2])
      csv$width_m = csv$width_m + w_perturb
    }
    
    # remove GRWL data with widths<90m, elev<0m, lakes, canals:
    keep = csv$width_m>wMin & 
      csv$elev_m>minElev & 
      csv$lakeFlag!=1 & 
      csv$lakeFlag!=3 
    csv = csv[keep,]
    
    # calc river distance, width, surface area, sum river area, max discrete area value: 
    # calculate distance between each adjacent GRWL centerline cell:
    N  = nrow(csv)
    d = distCalc(csv)
    w = csv$width_m/reducer
    A_raw = w*d
    A = A_raw[A_raw>Amin]
    sumA = sum(A)
    Amax = max(A)
    Alen = length(A)
    # calculate the average number of channels: 
    nChan = mean(csv$nchannels)
    
    # Maximum Liklihood Estimation (MLE) fit:
    
    # set std=T to calc. std dev of fit using MLE optimization,
    # which is slow and produces negligible std dev vals:
    bFit = pareto.mle(x=A, std=F) 
    
    # calculate GOF statistics with X2 and KS test:
    h = hist(A[A<max(figBreaks/reducer)], figBreaks/reducer, plot=F)
    jA = jitter(A) # to prevent ties in the following GOF tests
    bGOF = suppressWarnings(GOF(bFit, h, jA))
    
    # add histogram information to table to be in plot:
    if (j == 1){
      hTab = as.data.frame(array(0, c(nRun, length(figBreaks))))
      names(hTab) = round(figBreaks)
    }
    hTab[j,1:length(h$counts)] = h$counts
    
    # RSSA estimate with uncertainty: 
    # calculate total surface area of rivers and streams by extending 
    # RSSA abundance to smaller rivers, and calc confidence intervals
    # using a pareto fit (as in Class B basins):
    bRSSAextrap = RSSAextrapolater(bFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    # using mean of Class A basin fits: 
    aRSSAextrap = RSSAextrapolater(aFit, fOaMean, fOaSD, Amin, Amax, N=1, sumA)
    
    
    # calculate the % of land surface occupied by rivers and streams: 
    # obs_pc = 100*sumA*reducer*1e-6/BA[i]
    bRSSA_pc = 100*bRSSAextrap$meanRSSA/BA[i]
    bRSSA_MC_pc = 100*bRSSAextrap$MCRSSA/BA[i]
    aRSSA_pc = 100*aRSSAextrap$meanRSSA/BA[i]
    aRSSA_MC_pc = 100*aRSSAextrap$MCRSSA/BA[i]
    
    # if not a MC run, remove all rows but the first one:
    if (nRun == 1){ ensembleTab = ensembleTab[1,] }
    
    # fill in table with outputs of ensemble run:
    ensembleTab[j,] =
      as.vector(c(
        classAB_fN[i], 
        BA[i],
        length(which(keep)), 
        bClass,
        nChan,
        Amax,
        sum(reducer*mL*1e-6*csv$width_m[csv$lakeFlag!=1 & csv$lakeFlag!=3]), 
        sumA*reducer*1e-6,
        # Class B basins:
        aRSSAextrap$meanRSSA, 
        aRSSAextrap$MCRSSA, 
        aRSSA_pc, 
        aRSSA_MC_pc, 
        aFit$xm*reducer, 
        aRSSAextrap$MCAlpha,
        aRSSAextrap$MCfOw,
        # Class A basins:
        bRSSAextrap$meanRSSA, 
        bRSSAextrap$MCRSSA, 
        bRSSA_pc, 
        bRSSA_MC_pc, 
        bFit$xm*reducer, 
        bRSSAextrap$MCAlpha,
        bRSSAextrap$MCfOw,
        bGOF$X2, 
        bGOF$X2_p, 
        bGOF$KS_D, 
        bGOF$KS_p
      ))
  } # end Monte Carlo simulation

  return(list(
    hTab=hTab,
    ensembleTab = ensembleTab
    ))
  
}

# weighted multiple power-law regression for Class A basins:
fitFun <- function(x1, x2, fit, confInt=NA){
  # if a confidence interval is specified, return a table with mean fit
  # with confidence intervals, otherwise only return mean fit vector:
  if (!is.na(confInt)){
    return(exp(predict(fit, data.frame(x1=x1, x2=x2), interval="confidence", level=confInt)))
  }else{
    return(as.vector(exp(predict(fit, data.frame(x1=x1, x2=x2)))))
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
#abline(0, 1, lty=2); abline(mod, col=2)
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
h = hist(dif, hSeq, 
     main = "Difference Histogram",
     xlab="GRWL Width - In Situ Width (m)")
lines(lineSeq, fn, col=4)
abline(v=nfit[1], col=4, lwd=1.4)
abline(v=c(nfit[1]-nfit[2], nfit[1]+nfit[2]), col=4, lwd=1.4, lty=2)
print(paste0("GRWL error Monte Carlo  mean = ", round(nfit[1]), 
             ",  StDev = ", round(nfit[2])))


##############################################################################
# Class A: Calculate RSSA (aRSSA)
##############################################################################
# In Class A basins, We use a Monte Carlo simulation to characterize 
# RSSA error due to the uncertainty of GRWL measurements. For each ensemble 
# run, we fit a Pareto distribution of RSSA measurements using maximum
# likelihood estimation (MLE). To fit the distribution, we set the 
# Pareto scale parameter, xm, equal to the lowermost observed RSSA value  
# and solve for the Pareto shape parameter. This lowermost observed RSSA, 
# Omin=3259 m2, is the product of the minimum observed river width threshold 
# (90 m) and the mean distance between river centerline pixels in 30-m 
# resolution Landsat imagery (36.21 m). 
# Quantify the distribution of these statistical fits

# read in hBASIN shapefile dbf:
hBASIN = foreign::read.dbf(sub('hybas_allMain', 'hybas_allMainCopy', hydroBASINpath))

# if list file does not exist, create a list N GRWL obs. in each hBasin: 
if (!file.exists(nGRWLperBasinOutPath)){
  nGRWLperBasListGenerator(fP, fN, tabDirPath, wMin, minElev, nGRWLperBasinOutPath)
}
# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")


# Class A basins contain >250,000 river measurements:
bClass = "1"
classA_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>250000])
classA_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>250000]
#mnTabP = paste0(wd, 'output/figs/figS6_SAfits_classAB_MCmn.csv')
#sdTabP = paste0(wd, 'output/figs/figS6_SAfits_classAB_MCsd.csv')
hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classA_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classA_fP)
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsMnTab.csv')

nClassA = length(classA_fP)
mnHtab = as.data.frame(array(0, c(nClassA, length(figBreaks))))
names(mnHtab) = round(figBreaks)
mnEnsTab = data.frame(array(NA, c(nClassA, length(ensembleTabNames)))) 
names(mnEnsTab) = ensembleTabNames


print(paste("Class A Basins:", paste(classA_fN, collapse=" ")))

# get basin area in km2:
BA = hBASIN$area_km2[match(classA_fN, hBASIN$MAIN_BAS)] 

# get start time:
old <- Sys.time() 

# calculate RSSA in the basin:
for (i in 1:nClassA){
  
  print(paste('i =',i, ' of ', nClassA,', fN = ', classA_fN[i]))
  
  # read in original GRWL data by basin:
  csv_raw = read.csv(classA_fP[i], header=T)
  N_raw = nrow(csv_raw)
  
  # use original unperturbed GRWL data to get mean RSSA estimate:
  print("Getting raw unperturbed parameters...")
  mnList = classAB_MC(nEnsRun=1, bClass=1)
  
  # concatenate mnList into one large table for all Class A basins:
  mnHtab[i,] = mnList$hTab
  mnEnsTab[i,] = mnList$ensembleTab
  
  # Monte Carlo error propogation to estimate uncertainty:
  print("Running Monte Carlo Simulation to get parameter uncertainty..")
  MClist = classAB_MC(nEnsRun=nRun, bClass=1)
  
  # write out ensemble histogram table:
  write.csv(MClist$hTab, hTabPath[i], row.names=F)
  # write out ensemble output table:
  write.csv(MClist$ensembleTab, ensembleTabPath[i], row.names=F)
  
  # # plots:
  # hist(as.numeric(MClist$ensembleTab$bMCRSSA_km2))
  # abline(v=MClist$ensembleTab$bRSSA_km2)
  # # Sensitivity scatter for first order width threshold:
  # plot(MClist$ensembleTab$pMCfOw, MClist$ensembleTab$bMCRSSA_km2)
  # # Sensitivity scatter for pareto alpha (slope param):
  # plot(MClist$ensembleTab$pMCalpha, MClist$ensembleTab$bMCRSSA_km2)
  # # Sensitivity scatter for observed river surface area:
  # plot(MClist$ensembleTab$obRSSAgt90_km2, MClist$ensembleTab$bMCRSSA_km2)
  
} # end basin RSSA calculation

new = Sys.time() - old
print(new) 

# write out mnHtab and mnEnsTab:
names(mnHtab) = round(figBreaks)
write.csv(mnHtab, mnHtabPath, row.names=F)
write.csv(mnEnsTab, mnEnsTabPath, row.names=F)

# round columns and write out mean and stdev output tables: 
# mnTabRound = tabRounder(mnTab)
# sdTabRound = tabRounder(sdTab)
# write.csv(mnTabRound, mnTabP, row.names=F)
# write.csv(sdTabRound, sdTabP, row.names=F)


# # add up the total observed & extrapolated river surface area in Class A basins:
# classA_obs_gt90m = sum(mnTab$obRSSAgt90_km2)
# classA_obs = sum(mnTab$obRSSA_km2)
# classA_RSSA = sum(mnTab$bRSSA_km2)
# classA_BasinA = sum(mnTab$BA_km2)
# print(paste("Class A basins observed %RSSA widths >90m:", round(100*classA_obs_gt90m/classA_BasinA,2), "%"))
# print(paste("Class A basins observed %RSSA all widths:", round(100*classA_obs/classA_BasinA,2), "%"))
# print(paste("Class A basins extrapolated %RSSA:", round(100*classA_RSSA/classA_BasinA,2), "%"))
# 
# # add uncertainty:
# RSSA_km2 = round(c(sum(as.numeric(mnTab$bRSSA_km2)), 
#                    + sum(as.numeric(sdTab$bRSSA_km2)), 
#                    sum(as.numeric(sdTab$bRSSA_km2))))
#   
#   
# print("extrapolated Class A RSSA:")
# print(RSSA_km2)
# print(100*RSSA_km2/sum(mnTab$BA_km2))

##############################################################################
# Class A: Plot Fig. S6: Pareto fits 
##############################################################################

# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")
classA_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>250000])
classA_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>250000]

# get paths of ensemble tables: 
hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classA_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classA_fP)

# read in non-perturbed mnHistTab and mnEnsTab used to establish mean alpha:
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsMnTab.csv')
mnHtab = read.csv(mnHtabPath, header=T)
mnEnsTab = read.csv(mnEnsTabPath, header=T)

# read in table with Class A names:
classA_bNames = read.csv(classA_bNamesPath, header=T)
basinOrder = match(classA_bNames$bID, classA_fN)

# define plotting bounds:
xlim = range(figBreaks)
ylim = c(1, 3e6)

# set up output pdf:
pdfOut = paste0(wd, 'output/figs/figS6_ClassA_fits.pdf')
pdf(pdfOut, width=6, height=4)
par(mfrow=c(4,5))
par(oma=c(3,3.5,0,0.5), mar=c(0,0,0,0))

# for each class A basin, plot RSSA histogram:
for (h in 1:nClassA){
  
  # plot in order of map in Fig. S6:
  i = basinOrder[h]
  
  # read in unperturbed mean table to get mean alpha:
  mnTab = mnEnsTab[i,]
  pTotA = mnTab$nObs*int*reducer
  
  # read in table containing stastistical results of each Monte Carlo simulation
  # run and take stdev of ensembles:
  ensembleTab = read.csv(ensembleTabPath[i], header=T)
  sdTab = as.data.frame(t(apply(ensembleTab, 2, sd)))
  
  # read in binned histogram data and take mean and stdev of ensembles:
  dTab = read.csv(hTabPath[i], header=T)
  dSD = (apply(dTab, 2, sd))
  dMn = mnHtab[i,]
  dMn[dMn < ylim[1] & dMn!=0] = ylim[1]
  
  # get histogram breaks:
  lB = c(figBreaks[-length(figBreaks)])
  rB = figBreaks[-1]
  mid = colMeans(rbind(lB, rB))
  lD = dMn - dSD
  lD[lD < ylim[1]] = ylim[1]
  uD = dMn + dSD
  uD[uD < ylim[1]] = ylim[1]
  zC = dMn>0 # & lD>0 & uD>0

  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab="", 
       ylab="",
       bty='n',
       type='n', 
       log='xy',
       xaxt='n',
       yaxt='n',
       las=T)
  box(lwd=0.5) 
  if(h %in% c(16:20)){
    xAxis = c(1e4, 1e5, 1e6)
    axis(1, at=xAxis, labels=formatC(xAxis, digits=0, format="e"), 
         tcl=-0.5, lwd=0.5, cex.axis=0.7)
  }
  if(h %in% c(1,6,11,16)){
    yAxis = c(1e0, 1e2, 1e4, 1e6)
    axis(2, at=yAxis, labels=formatC(yAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  
  polygon(x=rbind(lB[zC], rB[zC], rB[zC], lB[zC], NA),
          y=rbind(dMn[zC], dMn[zC], 
                  ylim[1], ylim[1], NA),
          bty="n", lwd=0.7, border=NA, col=rgb(0.7, 0.7, 0.7, 1))
  
  # add standard deviation (1-alpha) bar to each bin:
  segments(x0=mid[zC], x1=mid[zC], y0=lD[zC], y1=uD[zC], 
           col=rgb(1, 0, 0, 0.5), lwd=0.4)

  # Pareto fit over observed data:
  # # FIX ME: rerun MC ensemble then remove reducer multiplier below:
  y1 = dpareto(mnTab$pXmin, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  segments(mnTab$pXmin, y1, mnTab$Amax*reducer, y2, col=4, lwd=1)
  
  # add legend:
  legend("topright",
         c(paste0(substr(paste0("0", h), nchar(h), nchar(h)+1), ': ', classA_bNames$bName[h]),
           paste0("a = ", round(mnTab$pMCalpha, 2), "±", formatC(round(sdTab$pMCalpha,4), format="g"))),
           xjust=0, text.col=c(1,4), bty="n", cex=0.7)
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)


##############################################################################
# Class A: Plot Fig. 3C: Pareto extrapolations
##############################################################################

# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")
classA_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>250000])
classA_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>250000]

# read in ensemble plots:
hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classA_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classA_fP)

# read in non-perturbed mnHistTab and mnEnsTab used to establish mean alpha:
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsMnTab.csv')
mnHtab = read.csv(mnHtabPath, header=T)
mnEnsTab = read.csv(mnEnsTabPath, header=T)

# read in table with Class A names:
classA_bNames = read.csv(classA_bNamesPath, header=T)
basinOrder = match(classA_bNames$bID, classA_fN)

# define plotting bounds:
xlim = c(fOaMeans[1]*reducer, max(figBreaks))
ylim = c(1, 1e12)



# set up output pdf:
pdfOut = paste0(wd, 'output/figs/figS3b_ClassA_extraps.pdf')
pdf(pdfOut, width=6, height=4)
par(mfrow=c(4,5))
par(oma=c(3,3.5,0,0.5), mar=c(0,0,0,0))

# for each class A basin, plot RSSA histogram:
for (h in 1:nClassA){
  
  # plot in order of map in Fig. S6:
  i = basinOrder[h]
  
  # read in unperturbed mean table to get mean alpha:
  mnTab = mnEnsTab[i,]
  pTotA = mnTab$nObs*int*reducer
  
  # read in table containing stastistical results of each Monte Carlo simulation
  # run and take stdev of ensembles:
  ensembleTab = read.csv(ensembleTabPath[i], header=T)
  sdTab = as.data.frame(t(apply(ensembleTab, 2, sd)))
  
  # read in binned histogram data and take mean and stdev of ensembles:
  dTab = read.csv(hTabPath[i], header=T)
  dSD = (apply(dTab, 2, sd))
  dMn = as.numeric(mnHtab[i,])
  dMn[dMn < ylim[1] & dMn!=0] = ylim[1]
  
  # get histogram breaks:
  lB = c(figBreaks[-length(figBreaks)])
  rB = figBreaks[-1]
  mid = colMeans(rbind(lB, rB))
  lD = dMn - dSD
  lD[lD < ylim[1]] = ylim[1]
  
  uD = dMn + dSD
  uD[uD < ylim[1]] = ylim[1]
  zC = dMn>0 # & lD>0 & uD>0
  
  # plot:
  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab="", 
       ylab="",
       bty='n',
       type='n', 
       log='xy',
       xaxt='n',
       yaxt='n',
       las=T); box(lwd=0.5) 
  if(h %in% c(16:20)){
    xAxis = c(1e1, 1e3, 1e5)
    axis(1, at=xAxis, labels=formatC(xAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  if(h %in% c(1,6,11,16)){
    yAxis = c(1e0, 1e6, 1e12)
    axis(2, at=yAxis, labels=formatC(yAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  # add histogram & std (1 sigma) segments to each bin:
  polygon(x=rbind(lB[zC], rB[zC], rB[zC], lB[zC], NA),
          y=rbind(dMn[zC], dMn[zC], 
                  ylim[1], ylim[1], NA),
          bty="n", lwd=0.7, border=NA, col=rgb(0.7, 0.7, 0.7, 1))
  segments(x0=mid[zC], x1=mid[zC], y0=lD[zC], y1=uD[zC], 
           col=rgb(0,0,0,0.5), lwd=0.4)
  
  # add Class A extrapolation polygon:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(mnTab$pXmin, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  polygon(x = c(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pXmin, fOaMeans[2]*reducer),
          y = c(y1, y2, ylim[1], ylim[1]),
          col=rgb(0.85,0.85,0.85,1), border=NA)
  
  # add Class A mean fit:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  segments(fOaMeans[2]*reducer, y1, mnTab$Amax*reducer, y2, col=1, lwd=0.8)
  
  # add first order width uncertainty segments: 
  y1 = dpareto(fOaMeans[1]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(fOaMeans[3]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  segments(fOaMeans[1]*reducer, ylim[1], fOaMeans[1]*reducer, y1, col=1, lwd=0.8, lty=3)
  segments(fOaMeans[3]*reducer, ylim[1], fOaMeans[3]*reducer, y2, col=1, lwd=0.8, lty=3)
  # error arrows at top:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(fOaMeans[3]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  arrows(fOaMeans[2]*reducer, y1, fOaMeans[3]*reducer, y2, 0.02, 90, col=1, lwd=0.8)
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(fOaMeans[1]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  arrows(fOaMeans[2]*reducer, y1, fOaMeans[1]*reducer, y2, 0.02, 90, col=1, lwd=0.8)
  # error arrows at bottom:
  y2 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  arrows(fOaMeans[2]*reducer, ylim[1], fOaMeans[1]*reducer, ylim[1], 0.02, 90, col=1, lwd=0.8)
  y2 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  arrows(fOaMeans[2]*reducer, ylim[1], fOaMeans[3]*reducer, ylim[1], 0.02, 90, col=1, lwd=0.8)
  
  # add labels:
  text(x=mid[1], 
       y=exp(mean(log(c(dMn[1], ylim[1])))/2),
       expression(bold(italic('Observed'))), pos=4, cex=0.6) 
  text(x=exp(mean(log(c(xlim[1], mnTab$pXmin)))),
       y=exp(mean(log(c(dMn[1], ylim[1])))/2),
       expression(bold(italic('Estimated'))), cex=0.6) 
  
  # add legend:
  legend("topright",
         c(paste0(substr(paste0("0", h), nchar(h), nchar(h)+1), ': ', classA_bNames$bName[h]),
           paste0("%SA: ", round(mnTab$pMCRSSA_pc, 2), "±", round(sdTab$pMCRSSA_pc, 2),"%"),
           paste0("a = ", round(mnTab$pMCalpha, 2), "±", formatC(round(sdTab$pMCalpha,4), format="g"))),
         xjust=0, text.col=c(1,1,4), bty="n", cex=0.7)
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)



##############################################################################
# Class B: Calculate RSSA
##############################################################################
# In Class B basins, Class B basins contain an intermediate amount of GRWL 
# data, specifically between 10,000 and 250,000 river measurements. The same 
# approach applied in Class A basins is also applied in Class B basins, except 
# that rather than fitting the Pareto shape parameter to each individual basin, 
# we use the mean and standard deviation of the Pareto shape parameter that was 
# fit in Class A basins (fig. S6). In these basins, we calculate 
# the combined uncertainty of the Pareto fit and the extrapolation minimum by 
# Monte-Carlo uncertainty propagation on the RSSA definite integral (Nruns=500)
# Quantify the distribution of these statistical fits


# read in hBASIN shapefile dbf:
hBASIN = foreign::read.dbf(sub('hybas_allMain', 'hybas_allMainCopy', hydroBASINpath))

# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")

# if previously reran analysis in Class A basins, 
# recalculate mean of raw fits and stdev of ensemble fits for Class A basins:
classA_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>250000])
classA_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>250000]
nClassA = length(classA_fP)

# read in unperturbed class A mnHistTab to establish mean alpha:
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsMnTab.csv')
mnEnsTab = read.csv(mnEnsTabPath, header=T)

# charactarize uncertainty of alpha by finding spread between ensemble runs:
for (i in 1:length(classA_fP)){
  if (i == 1){
    classAensembleTab = read.csv(ensembleTabPath[i], header=T)
  }else{
    classAensembleTab = rbind(classAensembleTab, read.csv(ensembleTabPath[i], header=T))
  }
}
sdTab = as.data.frame(t(apply(classAensembleTab, 2, sd)))

# create list of statistical parameters derived from class A basins:
aFit = list(
  xm = Amin*reducer,
  alpha = mean(mnEnsTab$pMCalpha),
  stdev = sdTab$pMCalpha
)
print(aFit)

# Class B basins contain between 10k and 250k river measurements:
bClass = "2"
classBboo = nGRWLperBasTab$nGRWLperBas>10000 & nGRWLperBasTab$nGRWLperBas<=250000
classB_fP = as.character(nGRWLperBasTab$fP[classBboo])
classB_fN = nGRWLperBasTab$fN[classBboo]
nClassB = length(classB_fP)
print(paste("Class B Basins:", paste(classB_fN, collapse=" ")))
#mnTabP = paste0(wd, 'output/figs/figSX_SAfits_classB_MCmn.csv')
sdTabP = paste0(wd, 'output/figs/figSX_SAfits_classB_MCsd.csv')
hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classB_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classB_fP)
# unperturbed histogram and parameter table:
mnHtab = as.data.frame(array(0, c(nClassB, length(figBreaks))))
names(mnHtab) = round(figBreaks)
mnEnsTab = data.frame(array(NA, c(nClassB, length(ensembleTabNames)))) 
names(mnEnsTab) = ensembleTabNames
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsMnTab.csv')

# get basin area in km2:
BA = hBASIN$area_km2[match(classB_fN, hBASIN$MAIN_BAS)] 


# get start time:
old <- Sys.time() 

# calculate RSSA in each basin:
for (i in 1:nClassB){
  
  print(paste('i =',i, ' of ', nClassB,', fN = ', classB_fN[i]))
  
  # read in GRWL data by basin:
  csv_raw = read.csv(classB_fP[i], header=T)
  N_raw = nrow(csv_raw)
  
  # use original unperturbed GRWL data to get mean RSSA estimate:
  print("Getting raw unperturbed parameters...")
  mnList = classAB_MC(nEnsRun=1, bClass=2)
  
  # concatenate mnList into one large table for all Class A basins:
  print("Running Monte Carlo Simulation to get parameter uncertainty..")
  mnHtab[i,] = mnList$hTab
  mnEnsTab[i,] = mnList$ensembleTab
  
  # Monte Carlo error propogation to estimate uncertainty:
  MClist = classAB_MC(nEnsRun=nRun, bClass=2)
  
  # write out ensemble histogram table:
  write.csv(MClist$hTab, hTabPath[i], row.names=F)
  # write out ensemble output table:
  write.csv(MClist$ensembleTab, ensembleTabPath[i], row.names=F)

}

new = Sys.time() - old
print(new) 

# write out mnHtab and mnEnsTab:
names(mnHtab) = round(figBreaks)
write.csv(mnHtab, mnHtabPath, row.names=F)
write.csv(mnEnsTab, mnEnsTabPath, row.names=F)


##############################################################################
# Class B: Plot Fig. SX: Pareto fits
##############################################################################

# read in non-perturbed mnHistTab and mnEnsTab used to establish mean alpha:
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsMnTab.csv')
mnHtab = read.csv(mnHtabPath, header=T)
mnEnsTab = read.csv(mnEnsTabPath, header=T)

hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classB_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classB_fP)

xlim = range(figBreaks)
ylim = c(1, 3e6)


# set up output pdf:
pdfOut = paste0(wd, 'output/figs/figSX_ClassB_fits.pdf')
pdf(pdfOut, width=6, height=4)
par(mfrow=c(4,5))
par(oma=c(3,3.5,0,0.5), mar=c(0,0,0,0))

# for each class A basin, plot RSSA histogram:
for (i in 1:nClassB){
  # read in unperturbed mean table to get mean alpha:
  mnTab = mnEnsTab[i,]
  pTotA = mnTab$nObs*int*reducer
  
  # read in table containing stastistical results of each Monte Carlo simulation
  # run and take stdev of ensembles:
  ensembleTab = read.csv(ensembleTabPath[i], header=T)
  sdTab = as.data.frame(t(apply(ensembleTab, 2, sd)))
  
  # read in binned histogram data and take mean and stdev of ensembles:
  dTab = read.csv(hTabPath[i], header=T)
  dSD = (apply(dTab, 2, sd))
  dMn = mnHtab[i,]
  dMn[dMn < ylim[1] & dMn!=0] = ylim[1]
  
  # get histogram breaks:
  lB = c(figBreaks[-length(figBreaks)])
  rB = figBreaks[-1]
  mid = colMeans(rbind(lB, rB))
  lD = dMn - dSD
  lD[lD < ylim[1]] = ylim[1]
  uD = dMn + dSD
  uD[uD < ylim[1]] = ylim[1]
  zC = dMn>0 # & lD>0 & uD>0

  # plot:
  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab="", 
       ylab="",
       bty='n',
       type='n', 
       log='xy',
       xaxt='n',
       yaxt='n',
       las=T)
  box(lwd=0.5) 
  xAxisInd = rep(c(1:5)+15, 50) + rep(20, 50)*rep(c(0:9), each=5)
  if(i %in% xAxisInd){
    xAxis = c(1e4, 1e5, 1e6)
    axis(1, at=xAxis, labels=formatC(xAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  yAxisInd = rep(c(1,6,11,16), 50) + rep(20, 50)*rep(c(0:49), each=4)
  if(i %in% yAxisInd){
    yAxis = c(1e0, 1e2, 1e4, 1e6)
    axis(2, at=yAxis, labels=formatC(yAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  
  polygon(x=rbind(lB[zC], rB[zC], rB[zC], lB[zC], NA),
          y=rbind(dMn[zC], dMn[zC], 
                  ylim[1], ylim[1], NA),
          bty="n", lwd=0.7, border=NA, col=rgb(0.7, 0.7, 0.7, 1))
  
  # add standard deviation (1-alpha) bar to each bin:
  segments(x0=mid[zC], x1=mid[zC], y0=lD[zC], y1=uD[zC], col=1, lwd=0.4)
  
  # add Class A mean fit:
  y1 = dpareto(mnTab$pXmin, mnTab$pXmin, mnTab$gMCalpha)*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, mnTab$gMCalpha)*pTotA
  segments(mnTab$pXmin, y1, mnTab$Amax*reducer, y2, col=2, lwd=1)
  
  # add Class B MLE fit:
  y1 = dpareto(mnTab$pXmin, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  segments(mnTab$pXmin, y1, mnTab$Amax*reducer, y2, col=4, lwd=1)


  # add legend:
  legend("topright",
         c(paste0(substr(paste0("0", i), nchar(i), nchar(i)+1), ': ', classB_fN[i]),
           paste0("a = ", round(mnTab$pMCalpha, 2), "±", formatC(round(sdTab$pMCalpha,4), format="g")),
           paste0("a = ", round(aFit$alpha, 2), "±", formatC(round(aFit$stdev,2), format="g"))),
         xjust=0, text.col=c(1,4,2), bty="n", cex=0.7)
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)




##############################################################################
# Class B: Plot Fig. 3D: Pareto extrapolations
##############################################################################

# read in non-perturbed mnHistTab and mnEnsTab used to establish mean alpha:
mnHtabPath = paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsHtab.csv')
mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsMnTab.csv')
mnHtab = read.csv(mnHtabPath, header=T)
mnEnsTab = read.csv(mnEnsTabPath, header=T)

xlim = range(fOaMeans[1]*reducer, figBreaks)
ylim = c(1, 1e12)

hTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleHistTabs/'), classB_fP)
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classB_fP)


# set up output pdf:
pdfOut = paste0(wd, 'output/figs/figS3b_ClassB_extraps.pdf')
pdf(pdfOut, width=6, height=4)
par(mfrow=c(4,5))
par(oma=c(3,3.5,0,0.5), mar=c(0,0,0,0))

# for each class A basin, plot RSSA histogram:
for (i in 1:nClassB){
  # read in unperturbed mean table to get mean alpha:
  mnTab = mnEnsTab[i,]
  pTotA = mnTab$nObs*int*reducer
  
  # read in table containing stastistical results of each Monte Carlo simulation
  # run and take stdev of ensembles:
  ensembleTab = read.csv(ensembleTabPath[i], header=T)
  sdTab = as.data.frame(t(apply(ensembleTab, 2, sd)))
  
  # read in binned histogram data and take mean and stdev of ensembles:
  dTab = read.csv(hTabPath[i], header=T)
  dSD = (apply(dTab, 2, sd))
  dMn = as.numeric(mnHtab[i,])
  dMn[dMn < ylim[1] & dMn!=0] = ylim[1]
  
  # get histogram breaks:
  lB = c(figBreaks[-length(figBreaks)])
  rB = figBreaks[-1]
  mid = colMeans(rbind(lB, rB))
  lD = dMn - dSD
  lD[lD < ylim[1]] = ylim[1]
  uD = dMn + dSD
  uD[uD < ylim[1]] = ylim[1]
  zC = dMn>0 # & lD>0 & uD>0
  # alternative: instead of mean and std, use 1st, 2nd, & 3rd quartile:
  # zC = quarts[1,]>0 & quarts[2,]>0 & quarts[3,]>0
  # arrows(mid[zC], quarts[1,zC], mid[zC], quarts[3,zC],
  #        code=3, length=0.018, angle=90, lwd=0.5)
  
  # plot RSSA histogram:
  # # quartile histogram:
  # plotLimInd = which(figBreaks<=xlim[2])
  # quarts = apply(dTab, 2, quantile, probs=c(0.25, 0.5, 0.75))[,plotLimInd[-length(figBreaks)]]
  plot(NA, 
       xlim=xlim, 
       ylim=ylim,
       xlab="", 
       ylab="",
       bty='n',
       type='n', 
       log='xy',
       xaxt='n',
       yaxt='n',
       las=T); box(lwd=0.5) 
  xAxisInd = rep(c(1:5)+15, 100) + rep(20, 100)*rep(c(0:9), each=5)
  if(i %in% xAxisInd){
    xAxis = c(1e1, 1e3, 1e5)
    axis(1, at=xAxis, labels=formatC(xAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  yAxisInd = rep(c(1,6,11,16), 100) + rep(20, 100)*rep(c(0:99), each=4)
  if(i %in% yAxisInd){
    yAxis = c(1e0, 1e6, 1e12)
    axis(2, at=yAxis, labels=formatC(yAxis, digits=0, format="e"), 
         tcl=-0.5, las=1, lwd=0.5, cex.axis=0.7)
  }
  # add histogram:
  polygon(x=rbind(lB[zC], rB[zC], rB[zC], lB[zC], NA),
          y=rbind(dMn[zC], dMn[zC], 
                  ylim[1], ylim[1], NA),
          bty="n", lwd=0.7, border=NA, col=rgb(0.7, 0.7, 0.7, 1))
  # add standard deviation (1-alpha) bar to each bin:
  segments(x0=mid[zC], x1=mid[zC], y0=lD[zC], y1=uD[zC], col=1, lwd=0.4)
  
  # add Class A extrapolation polygon:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  y2 = dpareto(mnTab$pXmin, mnTab$pXmin, aFit[[2]])*pTotA
  polygon(x = c(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pXmin, fOaMeans[2]*reducer),
          y = c(y1, y2, ylim[1], ylim[1]),
          col=rgb(0.85,0.85,0.85,1), border=NA)
  
  # add Class A fit uncertainty segments: 
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]]+aFit[[3]])*pTotA
  y2 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]]-aFit[[3]])*pTotA
  y3 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, aFit[[2]]+aFit[[3]])*pTotA
  y4 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, aFit[[2]]-aFit[[3]])*pTotA
  segments(fOaMeans[2]*reducer, y1, mnTab$Amax*reducer, y3, col=1, lwd=0.8, lty=3)
  segments(fOaMeans[2]*reducer, y2, mnTab$Amax*reducer, y4, col=1, lwd=0.8, lty=3)
  # add Class A mean fit:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  segments(fOaMeans[2]*reducer, y1, mnTab$Amax*reducer, y2, col=1, lwd=0.8, lty=1)
  
  # add first order width uncertainty segments: 
  y1 = dpareto(fOaMeans[1]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  y2 = dpareto(fOaMeans[3]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  segments(fOaMeans[1]*reducer, ylim[1], fOaMeans[1]*reducer, y1, col=1, lwd=0.8, lty=3)
  segments(fOaMeans[3]*reducer, ylim[1], fOaMeans[3]*reducer, y2, col=1, lwd=0.8, lty=3)
  # error arrows at top:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  y2 = dpareto(fOaMeans[3]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  arrows(fOaMeans[2]*reducer, y1, fOaMeans[3]*reducer, y2, 0.02, 90, col=1, lwd=0.8)
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  y2 = dpareto(fOaMeans[1]*reducer, mnTab$pXmin, aFit[[2]])*pTotA
  arrows(fOaMeans[2]*reducer, y1, fOaMeans[1]*reducer, y2, 0.02, 90, col=1, lwd=0.8)
  # error arrows at bottom:
  y2 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]]-aFit[[3]])*pTotA
  arrows(fOaMeans[2]*reducer, ylim[1], fOaMeans[1]*reducer, ylim[1], 0.02, 90, col=1, lwd=0.8)
  y2 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, aFit[[2]]+aFit[[3]])*pTotA
  arrows(fOaMeans[2]*reducer, ylim[1], fOaMeans[3]*reducer, ylim[1], 0.02, 90, col=1, lwd=0.8)
  
  # add labels:
  text(x=mid[1], 
       y=exp(mean(log(c(dMn[1], ylim[1])))/2),
       expression(bold(italic('Observed'))), pos=4, cex=0.6) 
  text(x=exp(mean(log(c(xlim[1], mnTab$pXmin)))),
       y=exp(mean(log(c(dMn[1], ylim[1])))/2),
       expression(bold(italic('Estimated'))), cex=0.6) 
  
  # add Class B MLE fit:
  y1 = dpareto(fOaMeans[2]*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  y2 = dpareto(mnTab$Amax*reducer, mnTab$pXmin, mnTab$pMCalpha)*pTotA
  segments(fOaMeans[2]*reducer, y1, mnTab$Amax*reducer, y2, col=4, lwd=0.5)
  
  # add legend:
  legend("topright",
         c(paste0(substr(paste0("0", i), nchar(i), nchar(i)+1), ': ', classB_fN[i]),
           paste0("%SA: ", round(mnTab$gMCRSSA_pc, 2), "±", round(sdTab$gMCRSSA_pc, 2),"%"),
           paste0("a = ", round(aFit[[2]], 2), "±", formatC(round(aFit[[3]],4), format="g")),
           paste0("a = ", round(mnTab$pMCalpha, 2), "±", formatC(round(sdTab$pMCalpha,4), format="g"))),
         xjust=0, text.col=c(1,1,1,4), bty="n", cex=0.7)
}

dev.off()
cmd = paste('open', pdfOut)
system(cmd)





##############################################################################
# Attach Class A & B attributes to hydroBASIN shapefile
##############################################################################
# Allows for the examination of RSSA spatial distribution in Class A&B Basins:

# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")

# run through all Class A and Class B esemble output files and generate
# a basin mean and std table:
classAB_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>10000])
classAB_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>10000]
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classAB_fP)

# read in non-perturbed mnEnsTab used to establish mean alpha for 
# Class A and Class B basins: 
mnEnsTabPathAB =  paste0(tabDirPath, 'nonEnsRunTabs/classAB_nonEnsMnTab.csv')
io = 0; if (io == 1){
  mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classA_nonEnsMnTab.csv')
  mnTabA = read.csv(mnEnsTabPath, header=T)
  mnEnsTabPath =  paste0(tabDirPath, 'nonEnsRunTabs/classB_nonEnsMnTab.csv')
  mnTabB = read.csv(mnEnsTabPath, header=T)
  mnTab = rbind(mnTabA, mnTabB)
  write.csv(mnTab, mnEnsTabPathAB, row.names=F)
}else{
  mnTab = read.csv(mnEnsTabPathAB, header=T)
}

# run this short calculation if Class A & B basins have been regenerated:
# read in table containing stastistical results of each Monte Carlo simulation
# run and take stdev of ensembles:
sdTabP = paste0(figDirPath, 'figSX_SAfits_classAB_MCsd.csv')

io = 0; if (io == 1){
  for (i in 1:length(classAB_fN)){
    print(i)
    enTab = read.csv(ensembleTabPath[i], header=T)
    if (i == 1){
      sdTab = as.data.frame(t(apply(enTab, 2, sd, na.rm=T)))
    }else{
      sdTab = rbind(sdTab, as.data.frame(t(apply(enTab, 2, sd, na.rm=T))))
    }
  }
  write.csv(sdTab, sdTabP, row.names=F)
}else{
  sdTab = read.csv(sdTabP, header=T)
}

# add column that contains the RSSA valuess used in manuscript:
classAboo = mnTab$bClass == 1
classBboo = mnTab$bClass == 2
mnTab$RSSA_km2[classAboo] = mnTab$bMCRSSA_km2[classAboo]
mnTab$RSSA_km2[classBboo] = mnTab$aMCRSSA_km2[classBboo]
sdTab$RSSA_km2[classAboo] = sdTab$bMCRSSA_km2[classAboo]
sdTab$RSSA_km2[classBboo] = sdTab$aMCRSSA_km2[classBboo]
# RSSA percentage basin:
mnTab$RSSA_pc = 100*mnTab$RSSA_km2/mnTab$BA_km2
sdTab$RSSA_pc = 100*sdTab$RSSA_km2/mnTab$BA_km2
# update headers of sdTab to distinguish them from mnTab headers:
names(sdTab) = paste0("sd_",names(mnTab))


# read in hBASIN dbf and attached mnTab & sdTab attributes to it:
hBASIN = foreign::read.dbf(sub('hybas_allMain.dbf', 'hybas_allMainCopy.dbf', hydroBASINpath))
bindTab = data.frame(array(NA, c(nrow(hBASIN), ncol(mnTab)+ncol(sdTab))))
names(bindTab) = c(names(mnTab), names(sdTab))
matchInd = match(mnTab$hBASIN_code, hBASIN$MAIN_BAS)
gMat = matrix(as.matrix(cbind(mnTab,sdTab)), ncol = ncol(mnTab)+ncol(sdTab), dimnames = NULL)
for (i in 1:length(matchInd)){
  bindTab[matchInd[i], ] = gMat[i,]
}

# format table for arcmap shapefile:
newhBASIN = as.data.frame(cbind(hBASIN, bindTab))
newhBASIN[is.na(newhBASIN)] = 0
newhBASIN = data.matrix(newhBASIN)
write.dbf(newhBASIN, hydroBASINpath)



##############################################################################
# Class C: Calculate RSSA 
##############################################################################
# Class C basins contain <10,000 GRWL width measurements ≥90 m, tend to be 
# small and/or dry basins. For these basins, we develop a relationship between 
# basin percent RSSA (%RSSA), basin area (BA), and aridity index (AI). 
# We use a least squares multiple linear regression of log-transformed data 
# weighted by basin area to interpolate RSSA in (Fig. 3E). Larger basins have 
# a larger percent %RSSA because they contain higher-order rivers. Uncertainty 
# in these basins is based on the 1σ confidence intervals of the multiple 
# linear regression. Use climate-%RSSA regression to fill in missing basins 
# in hydroBASINs shapefile:

# create MC ensemble concatenated table output path:
conTabP = paste0(wd, 'output/figs/figSX_SAfits_classAB_MCconTab.csv')

# read in and process hydroBasin table:
hBASIN = foreign::read.dbf(hydroBASINpath)
if ("fit" %in% names(hBASIN)){hBASIN = hBASIN[ ,(1:30)]}
# add FID column to hBASIN table to match up data:
if (!'FID' %in% names(hBASIN)){ 
  FID = 1:nrow(hBASIN)-1; hBASIN = cbind(FID, hBASIN) 
}

# read in sorted table of N GRWL obs in each hBasin: 
nGRWLperBasTab = read.csv(nGRWLperBasinOutPath, header=T)
nGRWLperBasTab$fP = paste0(wd, "input/GRWL/GRWL_by_hydroBASIN/", nGRWLperBasTab$fN, ".csv")

# run through all Class A and Class B ensemble output files and generate
# a basin mean and std table:
classAB_fP = as.character(nGRWLperBasTab$fP[nGRWLperBasTab$nGRWLperBas>10000])
classAB_fN = nGRWLperBasTab$fN[nGRWLperBasTab$nGRWLperBas>10000]
ensembleTabPath = sub(GRWLpath, paste0(tabDirPath, 'ensembleOutputTabs/'), classAB_fP)

# label Class C basins:
classABboo = hBASIN$bClass == 1 | hBASIN$bClass == 2
classCboo = !classABboo
hBASIN$bClass[classCboo] = 3

# To produce aridity table, downloaded Zomer et al., 2007 and 
# in ArcMap --> ArcToolbox --> Zonal Statistics as Table --> Export as DBF
# read in aridity table:
aridity = foreign::read.dbf(aridityPath)

# define regression inputs:
AI = aridity$MEAN[match(hBASIN$FID, aridity$FID_)] 
BA = hBASIN$area_km

# run weight multiple power law regression and get confidence interval 
# which will be included in uncertainty calc:
x1 = AI[classABboo]
x2 = BA[classABboo]
y = hBASIN$RSSA_pc[classABboo]
weight = hBASIN$area_km[classABboo]
fit = lm(log(y)~log(x1)+log(x2), weights=weight)
mnFitTab = as.data.frame(fitFun(AI, BA, fit, confInt = 0.68))

# create unperturbed mean & conf interval table:
mnFitL = fitFun(x1=xSeq1, x2=xSeq2, fit, confInt=confInt)

# # Aridity index vs River Area plot:
# # create plotting table:
# tab = data.frame(x1, x2, y, weight)
# # rescale dot radius for plot:
# dSize = 2*sqrt(weight/pi)
# dotSize = dSize - min(dSize)+100
# with(tab, symbols(x1, y,
#                   circles=dotSize, inches=0.12,
#                   bg=rgb(0,0,0,.3), fg=NA,
#                   xlim=xRange, ylim=yRange,
#                   main="", xlab="Aridity Index (AI)", ylab = "%RSSA", cex.lab=0.7,
#                   las=1, bty='n')); box(lwd=0.5)
# # plot regressions:

# lines(xSeq1, mnFitL[,1], lwd=1)
# lines(xSeq1, mnFitL[,2], lwd=1, lty=3) # lower confidence
# lines(xSeq1, mnFitL[,3], lwd=1, lty=3) # upper confidence



# get uncertainty estimate of regression using Monte Carlo simulation: 

# read in each ensemble tab and concatenate data into one large table:
# if new ClassAB ensembles were updated, concatenate ensembles into big tab:
io = 0; if (io == 1){
  print("Generating MC concat table...")
  for (i in 1:length(classAB_fN)){
    enTab = read.csv(ensembleTabPath[i], header=T)
    if (i == 1){conTab = enTab
    }else{ conTab = rbind(conTab, enTab) }
  }
  write.csv(conTab, conTabP, row.names=F)
}else{
  print("Reading in MC concat table...")
  conTab = read.csv(conTabP, header=T)
}

# add a column to conTab that contains the RSSA values used in manuscript:
classAboo = conTab$bClass == 1
conTab$RSSA_km2[classAboo] = conTab$bMCRSSA_km2[classAboo]
conTab$RSSA_km2[!classAboo] = conTab$aMCRSSA_km2[!classAboo]
conTab$RSSA_pc = 100*conTab$RSSA_km2/conTab$BA_km2


# Continue Monte Carlo simulation using the outputs of class A and B basins
# as the inputs to the  RSSA-BA-aridity relationship, run with nRun ensembles

# create index sequence:
iSeq = seq(1, nRun*length(classAB_fN), nRun)

# match ensemble dataset to hBasin and aridity dataset:
FIDmatchAB = match(conTab$hBASIN_code[iSeq], hBASIN$MAIN_BAS)
x1 = aridity$MEAN[match(hBASIN$FID[FIDmatchAB], aridity$FID_)]
x2 = hBASIN$area_km2[FIDmatchAB]
weight = hBASIN$area_km2[FIDmatchAB]

# for each ensemble run, develope RSSA-BA-aridity relationship:
modRSSAtabP = paste0(tabDirPath, "input_classC_MC.csv")
modRSSAtab = as.data.frame(array(NA, c(length(AI), nRun)))
mnClassBC = sdClassBC = rep(NA, length(nRun))
# for each ensembel run, develope a table for mean regression plot:
fitLmean = as.data.frame(array(NA, c(length(xSeq1), nRun)))


# # PLOT
# iSeq = seq(1, nRun*length(classAB_fN), nRun)
# y = conTab$RSSA_pc[iSeq]
# tab = data.frame(x1, x2, y, weight)
# yRange = c(0,3)
# xRange = c(0, 2.2e4)
# # Aridity index vs River Area plot:
# # rescale dot radius for plot:
# dSize = 2*sqrt(weight/pi)
# dotSize = dSize - min(dSize)+100
# 
# with(tab, symbols(x1, y,
#                   circles=dotSize, inches=0.12,
#                   bg=rgb(1,0,0,1), fg=NA,
#                   xlim=xRange, ylim=yRange,
#                   main="", xlab="Aridity Index (AI)", ylab = "%RSSA", cex.lab=0.7,
#                   las=1, bty='n')); box(lwd=0.5)
# fit = lm(log(y)~log(x1)+log(x2), weights=weight)
# fitFun1 = fitFun(x1=xSeq1, x2=xSeq2, fit, confInt=confInt)
# lines(xSeq1, fitFun1[,1], lwd=0.3, col=rgb(1,0,0,1))
# lines(xSeq1, fitFun1[,2], lwd=1, lty=3, col=rgb(1,0,0,1)) # lower confidence
# lines(xSeq1, fitFun1[,3], lwd=1, lty=3, col=rgb(1,0,0,1)) # upper confidence



# slower to read in this huge table than just generate it:
print('generating regression table...')  
# fitList = list()
for (i in 1:nRun){
  # print(i)
  iSeq = seq(i, nRun*length(classAB_fN), nRun)
  y = conTab$RSSA_pc[iSeq]
  
  # fit a weighted multiple linear regression to RSSA, aridity, BA: 
  fit = lm(log(y)~log(x1)+log(x2), weights=weight)
  # and add fit parameters to a list:
  # fitList[[i]] = fit
  
  # get a RSSA estimate for every basin (rows) for nRuns (cols):
  fitTab = fitFun(AI, BA, fit, confInt=NA)
  modRSSAtab[,i] = fitTab
  
  # using the above fit, generate three tables (mean regression,
  # upr, & lwr confidence intervals), for each ensemble:
  fitLtab = fitFun(xSeq1, xSeq2, fit, confInt=confInt)
  fitLmean[,i] = fitLtab[,1]
  
  # par(new=T)
  # with(tab, symbols(x1, y,
  #                   circles=dotSize, inches=0.12,
  #                   bg=rgb(0,0,1,0.1), fg=NA,
  #                   xlim=xRange, ylim=yRange,
  #                   main="", xlab="", ylab = "%RSSA", cex.lab=0.7,
  #                   las=1, bty='n', xaxt='n', yaxt='n')); box(lwd=0.5)
  # lines(xSeq1, fitLtab[,1], lwd=0.3, col=rgb(0,0,1,1))
  # lines(xSeq1, fitLtab[,2], lwd=0.3, lty=3, col=rgb(0,0,1,1)) # lower confidence
  # lines(xSeq1, fitLtab[,3], lwd=0.3, lty=3, col=rgb(0,0,1,1)) # upper confidence
  
}


# calculate total RSSA uncertainty for Class C basins:
# take mean of each unperturbed basins:
unPertModRSSAmn = mnFitTab[,1]
unPertModRSSAlwr = mnFitTab[,2]
unPertModRSSAupr = mnFitTab[,3]
unPertMnConfInt = rowMeans(cbind(unPertModRSSAmn-unPertModRSSAlwr, 
                               unPertModRSSAupr-unPertModRSSAmn))
# take standard deviation of each pertubed basin ensemble mean:
mnModRSSAuncert = t(as.data.frame(t(apply(modRSSAtab, 1, sd, na.rm=T))))
# calculate combined uncertainty:
sdModRSSA = unPertMnConfInt + mnModRSSAuncert
mnModRSSA = mnFitTab[,1]


# calculate total uncertainty for regression lines:
# take mean of each unperturbed basins:
unPertFitLmn = mnFitL[,1]
# take standard deviation of each pertubed basin ensemble mean:
mnFitLuncert = t(as.data.frame(t(apply(fitLmean, 1, sd, na.rm=T))))
# calculate combined uncertainty:
fitMean = c(mnFitL[,1])
fitLwr = c(mnFitL[,2] - mnFitLuncert)
fitUpr = c(mnFitL[,3] + mnFitLuncert)

# write out regression lines:
fitLinePath = paste0(figDirPath, "fig3E_regressionXY.csv")
write.csv(cbind(xSeq1, fitMean, fitLwr, fitUpr), fitLinePath, row.names=F)



# fill in hBASIN tabe:

# Fill in areas and uncertainty from area extrapolations and 
# combine error from area extrapolation and climate interpolation:
RSSA_pc = hBASIN$RSSA_pc
sd_RSSA_pc = hBASIN$sd_RSSA_pc
RSSA_km2 = hBASIN$RSSA_km2
sd_RSSA_km = hBASIN$sd_RSSA_km

RSSA_pc[classCboo] = mnModRSSA[classCboo]
sd_RSSA_pc[classCboo] = sdModRSSA[classCboo]

# set greenland icesheet to zero:
greenlandBoo = substr(hBASIN$MAIN_BAS, 1,1) == "9"
RSSA_pc[greenlandBoo] = 0
sd_RSSA_pc[greenlandBoo] = 0

# set any NA elements to zero since these are typically in small basins:
RSSA_pc[is.na(RSSA_pc)] = 0
sd_RSSA_pc[is.na(sd_RSSA_pc)] = 0

# calculate RSSA in terms of km2:
RSSA_km2[classCboo] = RSSA_pc[classCboo]*hBASIN$area_km2[classCboo]/100
sd_RSSA_km[classCboo] = sd_RSSA_pc[classCboo]*hBASIN$area_km2[classCboo]/100

# print out some important numbers and statistics:
print(range(RSSA_pc)) # range of basin %RSSA values
print(sum(RSSA_km2)) # global RSSA
print(sum(sd_RSSA_km)) # global RSSA error
print(100*sum(RSSA_km2)/globalLandArea) # global %RSSA
print(100*sum(sd_RSSA_km)/globalLandArea) # global RSSA error


# add final RSSA and uncertainty to dbf:
hBASIN$cRSSA_km2 = RSSA_km2
hBASIN$sd_cRSSA_km = sd_RSSA_km
hBASIN$cRSSA_pc = RSSA_pc
hBASIN$sd_cRSSA_pc = sd_RSSA_pc

# format table for arcmap shapefile and write out dbf:
hBASIN[is.na(hBASIN)] = 0
hBASIN = data.matrix(hBASIN)
write.dbf(hBASIN, hydroBASINpath)






# Compare results from previous studies:
raymond_mean = 536000
raymond_lwr = 399000
raymond_upr = 673000
downing_lwr = 585000
downing_upr = 662000
hBASIN_mean = sum(RSSA_km2)
hBASIN_lwr = sum(RSSA_km2) - sum(sd_RSSA_km)
hBASIN_upr = sum(RSSA_km2) + sum(sd_RSSA_km)

areaCompareTab = cbind(downing_lwr, downing_upr, raymond_mean, hBASIN_mean)
par(mar=c(8,5,1,1))
barplot(areaCompareTab, 
        ylim=c(0, max(c(areaCompareTab, hBASIN_upr))),
        ylab="Global River Area (km)",
        names.arg = c("Downing et al. \n lower",
                      "Downing et al. \n upper",
                      "Raymond et al.",
                      "This study"),
        cex.names=1, las=3, col="light gray")
arrows(3.1,raymond_mean,3.1,raymond_lwr, 0.1, 90)
arrows(3.1,raymond_mean,3.1,raymond_upr, 0.1, 90)
arrows(4.3,hBASIN_mean,4.3,hBASIN_lwr, 0.1, 90)
arrows(4.3,hBASIN_mean,4.3,hBASIN_upr, 0.1, 90)

perDif = c(round(100*(hBASIN_mean-downing_lwr)/hBASIN_mean), 
           round(100*(hBASIN_mean-downing_upr)/hBASIN_mean), 
           round(100*(hBASIN_mean-raymond_mean)/hBASIN_mean),
           round(100*(hBASIN_mean-raymond_lwr)/hBASIN_mean),
           round(100*(hBASIN_mean-raymond_upr)/hBASIN_mean))

text(c(1, 2, 3, 3, 3), as.numeric(cbind(downing_lwr, downing_upr, 
                                        raymond_mean, raymond_lwr, raymond_upr)), 
     paste(-perDif, "%"), pos=3)











##############################################################################
# Class C: Plot Fig. 3E Climate-RSSA multiple regresion
##############################################################################
# plot RSSA, BA and aridity (Fig. 3E):

# read in aridity table:
aridity = foreign::read.dbf(aridityPath)

# read in and process hydroBasin table:
hBASIN = foreign::read.dbf(hydroBASINpath)
# add FID column to hBASIN table to match up data:
if (!'FID' %in% names(hBASIN)){ 
  FID = 1:nrow(hBASIN)-1; hBASIN = cbind(FID, hBASIN) 
}

# determine which basins are not Class C:
classABboo = hBASIN$bClass != 3


# fit a weighted multiple linear regression to RSSA, aridity, BA:
x1 = aridity$MEAN[match(hBASIN$FID[classABboo], aridity$FID_)] 
x2 = hBASIN$area_km2[classABboo]
y = 100*hBASIN$RSSA_pc[classABboo]
sdY = 100*hBASIN$sd_RSSA_p[classABboo]
weight = hBASIN$area_km2[classABboo]

fit = lm(log(y)~log(x1)+log(x2), weights=weight); print(summary(fit))

# Plot multiple regression (Fig. 3E)
pdfOut = paste0(wd, 'output/figs/fig3E_ClassC_regression.pdf')
pdf(pdfOut,  width=3, height=2.4)
layoutTab = rbind(c(1,1,1,1,2))
layout(layoutTab)
par(mar=c(5.1,4.1,1,1))

# create plotting table:
tab = data.frame(x1, x2, y, weight)
yRange = c(0,5)
xRange = c(0, 2.2e4)
# Aridity index vs River Area plot:
# rescale dot radius for plot:
dSize = 2*sqrt(weight/pi)
dotSize = dSize - min(dSize)+100
with(tab, symbols(x1, y, 
                  circles=dotSize, inches=0.12, 
                  bg=rgb(0,0,0,.3), fg=NA, 
                  xlim=xRange, ylim=yRange,
                  main="", xlab="Aridity Index (AI)", ylab = "%RSSA", cex.lab=0.7,
                  las=1, bty='n', xaxt='n', yaxt='n')); box(lwd=0.5) 
segments(x1, y-sdY, x1, y+sdY,
         col=rgb(0,0,0,0.2))

xTx = seq(xRange[1], xRange[2], length.out=5)
xAxis = seq(xRange[1], xRange[2], length.out=3)
axis(1, at=xTx, labels=NA, lwd=0.5, cex.axis=0.7)
axis(1, at=xAxis, labels=formatC(xAxis, digits=0, format="e"), 
     lwd=0.5, cex.axis=0.7)
yTx =  seq(yRange[1], yRange[2], length.out=7)
yAxis =  round(seq(yRange[1], yRange[2], length.out=4))
axis(2, at=yTx, labels=NA, las=1, lwd=0.5, cex.axis=0.7)
axis(2, at=yAxis, labels=formatC(yAxis, digits=0, format="f"), las=1, 
     lwd=0.5, cex.axis=0.7)
# plot regressions:
xSeq1= seq(xRange[1], xRange[2], length.out=500) #xMin+((0:50)^10/50^10)*(xMax-xMin)
fitFun1 = fitFun(x1=xSeq1, x2=xSeq2, fit)
lines(xSeq1, fitFun1[,1], lwd=1)
lines(xSeq1, fitFun1[,2], lwd=1, lty=3) # lower confidence
lines(xSeq1, fitFun1[,3], lwd=1, lty=3) # upper confidence
# add regression equation:
text(xMin, yMax,  
     paste0("%RSSA = e^", round(round(fit[[1]][[1]], 1)),
            " * AI^", round(fit[[1]][[2]],2), 
            " * BA^", round(fit[[1]][[3]],2)), 
     cex=0.7, pos=4)

# add legend:
par(mar=c(4.1,1,2.1,1))
par(mai=c(0.4,0,0.2,0))

legendAreas = c(1e5, 5e5, 2e6, 6e6)
dSizeLeg = 2*sqrt(legendAreas/pi)
dotSizeLeg = dSizeLeg - min(dSizeLeg)+100
legendTab = data.frame(x=rep(1, length(legendAreas)), 
                       y=c(1:length(legendAreas)), legendAreas) 
suppressWarnings(with(legendTab, symbols(x, rev(y), circles=dotSizeLeg, 
                                         inches=0.12,  bg=rgb(0,0,0,0.5), fg=NA,
                                         ylim = c(0, length(legendAreas)+1),
                                         main="", axes=F, xlab='', ylab='')))
title("Basin\nArea\n(BA)", line=-2, cex.main=1)
options(scipen=2)
text(rep(1, length(legendAreas)), 
     c(1:length(legendAreas)) + 0.5, 
     rev(paste(labels=formatC(legendAreas, digits=0, format="e"), "km2")), 
     pos=1, offset=150*dotSizeLeg/max(dotSizeLeg)-3.5,
     cex=0.7)

dev.off() 
cmd = paste('open', pdfOut)
system(cmd)

summary(fit)

