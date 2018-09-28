library(INLA)
library(TMB)
library(fields) #use image.plot()
library(dplyr) #use image.plot()
library(animation) #use saveGIF()

#Comile and load c++ code-------
compile("./SPDExAR1_varying_stations/SPDExAR1_varying_stations.cpp")
dyn.load(dynlib("./SPDExAR1_varying_stations/SPDExAR1_varying_stations"))
#--------------------------------
data(PRborder)
k <- 12
data(PRprec)
coords <- as.matrix(PRprec[sample(1:nrow(PRprec)), 1:2])
source('./SPDExAR1_varying_stations/spde-book-functions.R')
params <- c(variance = 1, kappa = 1)
prdomain <- inla.nonconvex.hull(as.matrix(PRprec[, 1:2]),
                                convex = -0.03, concave = -0.05,
                                resolution = c(100, 100))
prmesh1 <- inla.mesh.2d(boundary = prdomain,
                        max.edge = c(0.7, 0.7), cutoff = 0.35,
                        offset = c(-0.05, -0.05))

set.seed(1)
x.k <- book.rspde(coords, range = sqrt(8) / params[2],
                  sigma = sqrt(params[1]), n = k, mesh = prmesh1,
                  return.attributes = TRUE)
dim(x.k)

c100 <- rainbow(101)
par(mfrow=c(4,3), mar=c(0,0,0,0))
for (j in 1:k)
  plot(coords, col=c100[round(100*(x[,j]-min(x[,j]))/diff(range(x[,j])))],
       axes=FALSE, asp=1, pch=19, cex=0.5)

rho <- 0.7
x <- x.k
for (j in 2:k)
  x[, j] <- rho * x[, j - 1] + sqrt(1 - rho^2) * x.k[, j]

n <- nrow(coords)
set.seed(2)
# ccov <- factor(sample(LETTERS[1:3], n * k, replace = TRUE))
# table(ccov)
# beta <- -1:1
x2 <- rnorm(n * k, 10, 1)
sd.y <- 0.1
#y <- beta[unclass(ccov)] + x + rnorm(n * k, 0, sd.y)
y <- x + 2*x2 + rnorm(n * k, 0, sd.y)

#tapply(y, ccov, mean)

isel <- sample(1:(n * k), n * k / 2)
dat <- data.frame(y = as.vector(y),
                  x2 = x2,
                  #w = ccov,
                  time = rep(1:k, each = n),
                  xcoo = rep(coords[, 1], k),
                  ycoo = rep(coords[, 2], k))[isel, ]

dat <- dat[order(dat$time),]
rownames(dat) <- 1:nrow(dat)

#Read data------------
borders <-read.table("./SPDExAR1_varying_stations/Piemonte_borders.csv",header=TRUE,sep=",")
Piemonte_data <-read.table("./SPDExAR1_varying_stations/Piemonte_data_byday.csv",header=TRUE,sep=",")
# #Convert dates to POSIX
Piemonte_data$Date <- as.Date(Piemonte_data$Date, format="%d/%m/%y")
# #Order the Piemonte data by date
Piemonte_data <- Piemonte_data[order(Piemonte_data$Date),]
#Piemonte_data <- Piemonte_data[-c(1:2000),]
# #Tweaking to create dataset with non-identical station locations per day------------

tweak_index <- which(Piemonte_data$Station.ID %in% c(2,10,15,21))
Piemonte_data[tweak_index,'UTMX'] <- Piemonte_data[tweak_index,'UTMX']+runif(length(tweak_index),-10,10)
Piemonte_data[tweak_index,'UTMY'] <- Piemonte_data[tweak_index,'UTMY']+runif(length(tweak_index),-10,10)

#Additionally remove some stations from some days (arbitrarily)
drop_stations_index <- sample.int(nrow(Piemonte_data), 1000)

Piemonte_data <- Piemonte_data[-drop_stations_index,]
rownames(Piemonte_data) <- 1:nrow(Piemonte_data)
#End of dataset tweaking------------

#coordinates <-read.table("./SPDExAR1_varying_stations/coordinates.csv",header=TRUE,sep=",")
#rownames(coordinates) = coordinates[,"Station.ID"]

#Coordinates are no longer same for all days and stations 
coordinates <- Piemonte_data[, colnames(Piemonte_data) %in% c("UTMX", "UTMY")]
plot(borders, cex=0.1)
points(coordinates)

#--------------------------------

#Define some variables used------
#n_stations <- length(coordinates$Station.ID)
n_stations <- as.data.frame(Piemonte_data %>% count(Date))
n_stations <- n_stations[order(n_stations$Date),]
n_stations <- n_stations[,2] # Vector with number of stations sampled per day

n_data <- length(Piemonte_data$Station.ID)
#n_days <- as.integer(n_data/n_stations)
n_days <- length(unique(Piemonte_data$Date))

# #dat
#n_stations <- length(coordinates$Station.ID)
n_stations2 <- as.data.frame(dat %>% count(time))
n_stations2 <- n_stations2[order(n_stations2$time),]
n_stations2 <- n_stations2[,2] # Vector with number of stations sampled per day

n_data2 <- length(dat$y)
#n_days <- as.integer(n_data/n_stations)
n_days2 <- length(unique(dat$time))

#--------------------------------

##Standardize covariates and calculate log of PM10--------
mean_covariates = apply(Piemonte_data[,3:10],2,mean)
sd_covariates = apply(Piemonte_data[,3:10],2,sd)
Piemonte_data[,3:10] =
  scale(Piemonte_data[,3:10],
        mean_covariates, sd_covariates)
Piemonte_data$logPM10 = log(Piemonte_data$PM10)
Piemonte_data$logPM10[which(is.na(Piemonte_data$logPM10))]=-999
#---------------------------------------------------------

#Calcualtes the days from first observation---------------
# Piemonte_data$day = rep(0,n_data)
# for(i in 1:n_days){
#   for(j in 1:n_stations){
#     Piemonte_data$day[(i-1)*n_stations + j] = i-1
#   }
# }

#This now can be done in a different way
Piemonte_data$day <- as.numeric(difftime(Piemonte_data$Date, min(Piemonte_data$Date), units="days"))
dat$day <- dat$time-1
#---------------------------------------------------------

#Group the days into intervalls used as time steps in the AR1 process--------
dtLength = 21#Length of time intervalls
dtLength2 = 1#Length of time intervalls
Piemonte_data$dt = floor(Piemonte_data$day/dtLength)
dat$dt = dat$day
#lengthDt = sum(Piemonte_data$dt==0)#Number of observations within each time intervall
#lengthDt <- as.data.frame(Piemonte_data %>% group_by(dt) %>% count(Station.ID))
#lengthDt <- lengthDt[,3]
#lengthDt2 <- as.data.frame(dat %>% group_by(dt) %>% count(Station.ID))
#lengthDt2 <- lengthDt2[,3]

maxDt = max(Piemonte_data$dt)+1
maxDt2 = max(as.numeric(dat$dt))+1
#---------------------------------------------------------

#Defines the GMRF and represenation of the sparce precision matrix------
# locStations = cbind(coordinates$UTMX, coordinates$UTMY)

locStations = cbind(coordinates$UTMX, coordinates$UTMY)


mesh =
  inla.mesh.2d(loc=locStations,
               loc.domain=borders,
               offset=c(10, 140),
               max.edge=c(40, 80),
               min.angle=c(26, 21),
               #cutoff=0
               cutoff=10
  )

spde = inla.spde2.matern(mesh=mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")]
A = inla.spde.make.A(mesh,locStations)
#aLoc = rep((1:dim(locStations)[1]),dtLength)-1


# 
locStations2 = coords

mesh2 =
  prmesh1

spde2 = inla.spde2.matern(mesh=mesh2, alpha=2)
spdeMatrices2 = spde2$param.inla[c("M0","M1","M2")]
A2 = inla.spde.make.A(mesh2,locStations2)
#aLoc2 = rep((1:dim(locStations2)[1]),dtLength2)-1

#---------------------------------------------------------

#Plot the triangulation-----------------------------------
plot(mesh)
lines(borders, lwd=3)
points(locStations, pch=20, cex=1, col=2)
# 
plot(mesh2)
lines(prdomain, lwd=3)
points(locStations2, pch=20, cex=1, col=2)
#---------------------------------------------------------

#Define the design matrix for the fixed effects-------
X <- model.matrix( ~ 1 + Piemonte_data$A + Piemonte_data$UTMX + Piemonte_data$UTMY +
                     Piemonte_data$WS + Piemonte_data$TEMP + Piemonte_data$HMIX +
                     Piemonte_data$PREC + Piemonte_data$EMI , data = Piemonte_data)

X2 <- model.matrix( ~ 1 + dat$x2 , data = dat)
#-----------------------------------------------------

#Define data and parameter object--------------------
data <- list(logPM10 = Piemonte_data$logPM10,
             #n_data = n_data,
             #lengthDt = lengthDt,
             maxDt = maxDt,
             A = A,
             #aLoc = aLoc,
             X = as.matrix(X),
             spdeMatrices = spdeMatrices,
             time_index = as.integer(Piemonte_data$dt)
)

parameters <- list(beta      = c(3,rep(0,8)),
                   log_tau   = 4,
                   log_kappa = -4,
                   rhoTan = 0,
                   x  = array(0,dim = c(mesh$n,maxDt)),
                   logSigmaE = -1)

data2 <- list(logPM10 = dat$y,
             #n_data = n_data2,
             #lengthDt = lengthDt2,
             maxDt = maxDt2,
             A = A2,
             #aLoc = aLoc2,
             X = as.matrix(X2),
             spdeMatrices2 = spdeMatrices2,
             time_index2 = as.integer(dat$dt)
)

parameters2 <- list(beta      = c(2),
                   log_tau   = 4,
                   log_kappa = -4,
                   rhoTan = 0,
                   x  = array(0,dim = c(mesh2$n,maxDt2)),
                   logSigmaE = -1)

#-----------------------------------------------------

#Fit model with SEPARABLE formulation-----------------
startTime <- Sys.time()
#data$flag = 1
data2$flag = 1

obj <- TMB::MakeADFun(data2,parameters2,random = c("x"),DLL = "SPDExAR1_varying_stations")
obj <- normalize(obj, flag="flag")
opt<-stats::nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
rep<- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)

# startTime <- Sys.time()
# data$flag = 1
# obj <- TMB::MakeADFun(data,parameters,random = c("x"),DLL = "SPDExAR1_varying_stations")
# obj <- normalize(obj, flag="flag")
# opt<-stats::nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=1000, iter.max=1000))
# rep<- sdreport(obj)
# endTime <- Sys.time()
# timeUsed = endTime - startTime
# print(timeUsed)
#-----------------------------------------------------


#Extract the range---------------
rangeIndex = which(row.names(summary(rep,"report"))=="range")
range = summary(rep,"report")[rangeIndex,]
#--------------------------------

#Construct .git file with spatio-termporal illustration------
# ani.options(convert = "C:/Program Files/ImageMagick/convert.exe")
# 
# system("where convert", intern = TRUE)
# system('"C:\\ImageMagick-7.0.8-Q16\\convert.exe"')
ani.options('autobrowse'= FALSE)
saveGIF({
  for(dtPlot in 1:(maxDt-1))
  {
    proj = inla.mesh.projector(mesh)
    gammaindeks = which(names(rep$par.random)=="gamma")
    gammaindeks = which(names(rep$par.random)=="xVec")
    XYUTM = mesh$loc[,1:2]
    indexOfGammaToPlot = ((dtPlot)*mesh$n+1):((dtPlot+1)*mesh$n)
    gamma = rep$par.random[indexOfGammaToPlot]/exp(rep$par.fixed[which(names(rep$par.fixed)=="log_tau")])
    
    latentFieldML = gamma
    image.plot(proj$x,proj$y, inla.mesh.project(proj, latentFieldML),col =  colorRampPalette(c("white","yellow", "red"))(12),
               xlab = 'Easting', ylab = 'Northing',
               zlim = c(-2,2),
               main = paste("MAP estimate of the spatio-temporal field day ", (dtPlot-1)*dtLength +1 ," to ",dtPlot*dtLength,sep = ""),
               cex.lab = 1.1,cex.axis = 1.1, cex.main=1, cex.sub= 1.1)
    contour(proj$x, proj$y,inla.mesh.project(proj, latentFieldML),nlevels = 6 ,add = T,labcex  = 1,cex = 1)
    lines(borders, lwd=3)
  }
}, movie.name = "./SPDExAR1_varying_stations/SPDExAR1_varying_stations.gif")
#-----------------------------------------------------------
