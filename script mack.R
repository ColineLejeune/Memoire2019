rm(list =ls())

#Chargement des librairies
library(readxl)
library(ChainLadder)
library(readr)
library (lubridate)
library(purrr)
library(MASS)
library(actuar)

#Chargement des data
Mack <- read_excel("...", col_names = FALSE)
data = Mack
data=as.matrix(data)
fct_triangle_payments = as.triangle(data)
fct_triangle_payments

#plot(fct_triangle_payments, lattice = TRUE)

# Ajustement du triangle pour que les facteurs de dev soint > ou = a 100

# fct_triangle_payments[1,12] = fct_triangle_payments[1,11]
# fct_triangle_payments[1,14] = fct_triangle_payments[1,13]
# fct_triangle_payments[1,13] = fct_triangle_payments[1,12]
# fct_triangle_payments[2,13] = fct_triangle_payments[2,12]
# fct_triangle_payments[1,14] = fct_triangle_payments[1,13]

####### Mack - Chain Ladder #######
mack_payments = MackChainLadder(Triangle = fct_triangle_payments, est.sigma = "Mack")

#facteurs de developpement
mack_payments$f 

#Triangle complet
mack_payments$FullTriangle

#Erreur de processus et de parametre
mack_payments$Total.ProcessRisk
mack_payments$Total.ParameterRisk

# Calcul des reserves
CDR(mack_payments)

####### Mack - Bootstrap #######

set.seed(0)
Bootstrap_payments = BootChainLadder (fct_triangle_payments, R=10000, process.distr="gamma") 
#Bootstrap_payments = BootChainLadder (fct_triangle_payments, R=10000) 
mean_BE =   (sum(Bootstrap_payments$IBNR.Totals))/10000
# plot de l'empirical cumulative distribution function
plot(ecdf(Bootstrap_payments$IBNR.Totals))

# plot sans les valeurs negatives
#plot(ecdf(Bootstrap_payments$IBNR.Totals[Bootstrap_payments$IBNR.Totals>0]))

#Il faut maintenant choisir la distribution empirique a choisir pour modeliser cette ecdf

## fit log-normal distribution
fit_payments <- fitdistr(Bootstrap_payments$IBNR.Totals[Bootstrap_payments$IBNR.Totals>0], "lognormal")
fit_payments
curve(plnorm(x,fit_payments$estimate["meanlog"], fit_payments$estimate["sdlog"]),col="red", add=TRUE)  

## fit normal distribution
fit_payments_N <- fitdistr(Bootstrap_payments$IBNR.Totals, "normal")
fit_payments_N
curve(pnorm(x,fit_payments_N$estimate["mean"], fit_payments_N$estimate["sd"]),col="blue", add=TRUE)  

## fit expo distribution
fit_payments_ex <- fitdistr(Bootstrap_payments$IBNR.Totals[Bootstrap_payments$IBNR.Totals>0], "exponential")
fit_payments_ex
curve(pexp(x,fit_payments_ex$estimate["rate"]),col="green", add=TRUE)

## fit weibull distribution
fit_payments_w <- fitdistr(Bootstrap_payments$IBNR.Totals[Bootstrap_payments$IBNR.Totals>0], "weibull")
fit_payments_w
curve(pweibull(x,fit_payments_w$estimate["shape"],fit_payments_w$estimate["scale"]),col="pink", add=TRUE)  

## fit Geom distribution
fit_payments_geo <- fitdistr(Bootstrap_payments$IBNR.Totals, "geometric")
fit_payments_geo
curve(pgeom(x,fit_payments_geo$estimate["prob"]),col="orange", add=TRUE)

#check avec la densite
plot(density(Bootstrap_payments$IBNR.Totals)$y, type = "l", col = "black")
plot(density(Bootstrap_payments$IBNR.Totals[Bootstrap_payments$IBNR.Totals>0])$y, type = "l", col = "black")

##lognormale
stock = rlnorm(1000000,fit_payments$estimate["meanlog"], fit_payments$estimate["sdlog"])
lines(density(stock)$y, col = "red")

##normal
stock_N = rnorm(1000000,fit_payments_N$estimate["mean"], fit_payments_N$estimate["sd"])
lines(density(stock_N)$y, col = "blue")

##expo
stock_ex = rexp(1000000,fit_payments_ex$estimate["rate"])
lines(density(stock_ex)$y, col = "green")

## weibull
stock_w = rweibull(1000000,fit_payments_w$estimate["shape"],fit_payments_w$estimate["scale"])
lines(density(stock_w)$y, col = "pink")

## geometric
stock_g = rgeom(1000000,fit_payments_geo$estimate["prob"])
lines(density(stock_g)$y, col = "orange")

#### en combinant l'analyse ecdf et densite, la log normale semble la plus appropriee

####### Risk Adjustment #######

#aggregate dist : Compute the aggregate claim amount cumulative distribution function of a portfolio over a period
# On simule 1000000 log normale de parametre  meanlog: 17.528144721 et sdlog: 0.545867631   

set.seed(0)
lognorm_payments_model <- expression(data = rlnorm(fit_payments$estimate["meanlog"] , fit_payments$estimate["sdlog"]))
Fstw <- aggregateDist("simulation", model.sev=lognorm_payments_model, nb.simul = 1000000)

#On calcule la VaR a 60% 
VAR_60= actuar::VaR(Fstw, conf.level=0.6)
RA_VAR60 = VAR_60 -mean_BE
taux_RA1 = RA_VAR60/mean_BE

#VaR a 65%
VAR_65 = VaR(Fstw, conf.level=0.65)
RA_VAR_65 = VAR_65 -mean_BE
taux_RA2 = RA_VAR_65/mean_BE

#TVaR
CTE = TVaR(Fstw, conf.level = 0.4)
RA_TVaR = CTE - mean_BE
taux_RA3 = RA_TVaR/mean_BE

#CoC
RA_COC = (VaR(Fstw, conf.level=0.9))*0.06
taux_RA4 = RA_COC/mean_BE

#Expected Shortfall
library(VaRES)
eslognorm(0.7, mu=fit_payments$estimate["meanlog"], sigma=fit_payments$estimate["sdlog"])

####### Hypothesis verification to apply Mack model #######

dev_factor = function(triangle)
{
  store_dev_f = c()
  for(i in (1:(nrow(triangle)-1)))
  {
    store_dev_f[i] = sum(triangle[1:(nrow(triangle)-i), i+1])/sum(triangle[1:(nrow(triangle)-i), i])
  }
  return (store_dev_f)
}


plotCij_Cik = function(fct_triangle_payments, j)
{
  plot(fct_triangle_payments[1,j-1], fct_triangle_payments[1,j], 
       ylab = paste("C i", j) , xlab =  paste("C i", j-1), main = paste("C i", j, "vs", "C i", j-1))
  
  for(i in 2:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    points(fct_triangle_payments[i,j-1], fct_triangle_payments[i,j])
  }
  abline(0, mack_payments$f[j-1])
}

plotWeightedRes_fk1 = function(fct_triangle_payments, j)
{
  
  dev_factor = dev_factor(fct_triangle_payments)
  maximus = 0
  minimus = Inf
  for(i in 1:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    maximus = max((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/sqrt(fct_triangle_payments[i,j-1]), maximus)
    minimus = min((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/sqrt(fct_triangle_payments[i,j-1]), minimus)
  }
  
  plot(fct_triangle_payments[1,j-1], (fct_triangle_payments[1,j] - fct_triangle_payments[1,j-1] * dev_factor[j-1])/sqrt(fct_triangle_payments[1,j-1]), 
       ylab = "weighted residuals", main = "Residual Plots for fk1", xlab =  paste("C i", j-1), ylim = c(minimus, maximus),
       xlim = c(min(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))]), 
                max(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))])))
  
  for(i in 2:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    points(fct_triangle_payments[i,j-1], (fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/sqrt(fct_triangle_payments[i,j-1]))
  }
}

plotWeightedRes_fk0 = function(fct_triangle_payments, j)
{
  #Implementation of equation (12) paper Mack
  dev_factor = rep(0, dim(fct_triangle_payments)[2])
  for(i in 2:dim(fct_triangle_payments)[2])
  {
    rectangle = fct_triangle_payments[1:(dim(fct_triangle_payments)[2]-i+1), (i-1):i]
    if(is.null(dim(rectangle)))
    {
      dev_factor[i-1] = (rectangle[2]*rectangle[1])/(rectangle[1]^2)
    }
    else{
      dev_factor[i-1] = sum(rectangle[,2]*rectangle[,1]) / sum(rectangle[,1]^2)
    }
  }
  
  maximus = 0
  minimus = Inf
  for(i in 1:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    maximus = max((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1]), maximus)
    minimus = min((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1]), minimus)
  }
  
  plot(fct_triangle_payments[1,j-1], (fct_triangle_payments[1,j] - fct_triangle_payments[1,j-1] * dev_factor[j-1]), 
       ylab = "weighted residuals", main = "Residual Plots for fk0"  ,xlab =  paste("C i", j-1), ylim = c(minimus, maximus), 
       xlim = c(min(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))]), 
                max(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))])))
  
  for(i in 2:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    points(fct_triangle_payments[i,j-1], (fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1]))
  }
}

plotWeightedRes_fk2 = function(fct_triangle_payments, j)
{
  #Implementation of equation (12) paper Mack
  dev_factor = rep(0, dim(fct_triangle_payments)[2])
  for(i in 2:dim(fct_triangle_payments)[2])
  {
    rectangle = fct_triangle_payments[1:(dim(fct_triangle_payments)[2]-i+1), (i-1):i]
    if(is.null(dim(rectangle)))
    {
      dev_factor[i-1] = (rectangle[2]/rectangle[1])
    }
    else{
      dev_factor[i-1] = (1/dim(rectangle)[1]) * sum(rectangle[,2]/rectangle[,1])
    }
  }
  
  maximus = 0
  minimus = Inf
  for(i in 1:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    maximus = max((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/(fct_triangle_payments[i,j-1]), maximus)
    minimus = min((fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/(fct_triangle_payments[i,j-1]), minimus)
  }
  
  plot(fct_triangle_payments[1,j-1], (fct_triangle_payments[1,j] - fct_triangle_payments[1,j-1] * dev_factor[j-1])/(fct_triangle_payments[1,j-1]), 
       ylab = "weighted residuals", main = "Residual Plots for fk2", xlab =  paste("C i", j-1), ylim = c(minimus, maximus),
       xlim = c(min(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))]), 
                max(c(fct_triangle_payments[,j-1])[!is.na(c(fct_triangle_payments[,j-1]))])))
  
  for(i in 2:length((fct_triangle_payments[,j])[!is.na(fct_triangle_payments[,j])]))
  {
    points(fct_triangle_payments[i,j-1], (fct_triangle_payments[i,j] - fct_triangle_payments[i,j-1] * dev_factor[j-1])/(fct_triangle_payments[i,j-1]))
  }
}

par(mfrow = c(2,2))
for(j in 2:14)
{
  plotCij_Cik(fct_triangle_payments = fct_triangle_payments, j = j)
  plotWeightedRes_fk1(fct_triangle_payments = fct_triangle_payments, j = j)
  plotWeightedRes_fk0(fct_triangle_payments = fct_triangle_payments, j = j)
  plotWeightedRes_fk2(fct_triangle_payments = fct_triangle_payments, j = j)
}


# Hypothesis 4) Independence

factorial = function(d)
{
  if(d == 1)
  {
    return(d)
  }
  else if(d == 0)
  {
    return(1)
  }
  else
  {
    return(d*factorial(d-1))
  }
}
combin = function(n,k)
{
  return(factorial(n)/(factorial(k) * (factorial(n-k))))
}

OrderingTriangle = function(fct_triangle_payments)
{
  factoredTriangle = fct_triangle_payments
  
  for(i in 2:dim(fct_triangle_payments)[2])
  {
    factoredTriangle[,i] = c((fct_triangle_payments[,i])/(fct_triangle_payments[,i-1]))
    
  }
  factoredTriangle = factoredTriangle[,-1]
  #factoredTriangle = factoredTriangle[-14,]
  factoredTriangle = factoredTriangle[-(nrow(fct_triangle_payments)),]
  
  orderedTriangle = factoredTriangle
  
  for(i in 1:dim(factoredTriangle)[2])
  {
    orderedTriangle[,i] = c(order(factoredTriangle[,i][!is.na(factoredTriangle[,i])]), rep(NA, dim(factoredTriangle)[2]-length(factoredTriangle[,i][!is.na(factoredTriangle[,i])])))
  }
  # return(orderedTriangle)
  stringedTriangle = orderedTriangle
  for(i in 1:dim(factoredTriangle)[2])
  {
    med = median(orderedTriangle[,i][!is.na(orderedTriangle[,i])])
    for(j in 1:length((factoredTriangle[,i])[!is.na(factoredTriangle[,i])]))
    {
      if(orderedTriangle[j,i] < med){
        stringedTriangle[j,i] = "S"   
      }
      else if(orderedTriangle[j,i] == med )
      {
        stringedTriangle[j,i] = "*"  
      }
      else if(orderedTriangle[j,i] > med)
      {
        stringedTriangle[j,i] = "L"
      }
    }
  }
  
  diagVectorS = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorL = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorZ = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorN = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorM = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorE = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  diagVectorV = rep(0, length(2:(nrow(fct_triangle_payments)-1)))
  
  Lengths = c(2:dim(factoredTriangle)[2])
  for(j in Lengths)
  {
    for(i in 1:j)
    {
      if(stringedTriangle[i,j+1-i] == "S"){
        diagVectorS[j-1] = diagVectorS[j-1] + 1
      }
      else if(stringedTriangle[i,j+1-i] == "L")
      {
        diagVectorL[j-1] = diagVectorL[j-1] + 1  
      }
      
    }
    diagVectorZ[j-1] = min(diagVectorS[j-1],diagVectorL[j-1])
    diagVectorN[j-1] = diagVectorS[j-1] + diagVectorL[j-1]
    diagVectorM[j-1] = floor((diagVectorN[j-1]-1)/2)
    diagVectorE[j-1] = (diagVectorN[j-1] /2) - combin(diagVectorN[j-1]-1, diagVectorM[j-1]) * (diagVectorN[j-1] / (2^diagVectorN[j-1]))
    diagVectorV[j-1] = ((diagVectorN[j-1] * (diagVectorN[j-1] - 1))/4) - combin(diagVectorN[j-1]-1, diagVectorM[j-1]) * ((diagVectorN[j-1] * (diagVectorN[j-1] - 1))/(2^diagVectorN[j-1])) + diagVectorE[j-1] - diagVectorE[j-1]^2
  }
  
  bindings = cbind(diagVectorS, diagVectorL, diagVectorZ, diagVectorN, diagVectorM, diagVectorE, diagVectorV)
  bindings = data.frame(bindings)
  return(list("bindings" = bindings, "order" = stringedTriangle))
}

keep = OrderingTriangle(fct_triangle_payments)
stat_test_Z = sum(keep$bindings$diagVectorZ)
IC_up = sum(keep$bindings$diagVectorE) + qnorm(0.975) * sqrt(sum(keep$bindings$diagVectorV))
IC_down = sum(keep$bindings$diagVectorE) - qnorm(0.975) * sqrt(sum(keep$bindings$diagVectorV))
IC = c(IC_down, IC_up) # and stat_test_Z is in the IC and so the independence is validated
IC
stat_test_Z
keep
