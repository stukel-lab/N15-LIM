#Note that first you have to set the directory: Session > Set Working Directory
rm(list=ls())  # Housekeeping

## Load libraries 
library(R.matlab)  # Used to read in/out matlab binary files
library(openxlsx)  # Used to read xlsx files
library(limSolve)  # The original LIM Package
library(MASS)      # Matrix library

source('xsampleN15.r')  # The N15 xsample algorithm 
source('ExternalFunctions.r')  # The place for ancillary functions

###############################
#### 1. Read in input data ####
###############################

Model <- 'NEMURO.Offshore'
datafile <- paste('N15InverseModelRW.',Model,'.RW.mat',sep="")
data <- readMat(datafile, fixNames=TRUE)
h <- data[['h']]
b <- data[['b']]
ba <- data[['ba']]
be <- data[['be']]
Ae <- data[['Ae']]
Aa <- data[['Aa']]
G <- data[['G']]
sdba <- data[['sdba']]
InputCol <- data[['InputCol']]
A <- data[['A']]
Inputs <- data[['Inputs']]

###############################
#### 2. Set run parameters ####
###############################

KeeperRatio <- 10000
IterLength <- 10000000
OutputLength <- IterLength / KeeperRatio
BurninLength <- IterLength / 100

## Set based on which forward model
if (Model == "DIAZO.Coastal"){
  jmpLength <- 0.02
}else if (Model == "DIAZO.Mesohaline") {
  jmpLength <- 0.005
}else if (Model == "NEMURO.Coastal") {
  jmpLength <- 0.001
}else if (Model == "NEMURO.Offshore") {
  jmpLength <- 0.0001
}
jmpLength15 <- 0.02

## Cut to isolate non-15N equations (11 -> 15 here)
Aa2 <- Aa[11:14,]
ba2 <- ba[11:14]
sdba2 <- sdba[11:14]
d15N0 <- c(0, 0, 0, 0, 0, 0)

## Read N15 values from forward model
upNO315N <- Inputs[8, InputCol]
Mes15N <- Inputs[9, InputCol]
Det15N <- Inputs[7, InputCol]
Doc15N <- Inputs[22, InputCol]
del15Nknown <- c(upNO315N, Mes15N, Doc15N, Det15N)

## L2MN Solution (without N15)
# Start
test <- lsei(A = Aa2, B = ba2, E = Ae, F = be, G = G, H = h, type=2) # No N15 here
X <- test[['X']]
solNorm <- test[['solutionNorm']]
lseisol <- as.matrix(X)
rm(test, X)  # Housekeeping


#Central Value Solution
center <- xranges(E = Ae, F = be, G = G, H = h,
                  ispos = FALSE, tol = 1e-8, central = TRUE, full=FALSE)
centralval <- center[,3]

##############################################################
#### 3a. Iterative burn-in in order to establish matricies ####
##############################################################
## First burnin
test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
                 H = h, sdB = sdba2*10, iter = BurninLength/2,
				 type="mirror", jmp=runif(1, 0, 1)*jmpLength/5,
				 x0 = centralval, fulloutput='TRUE')
				 
Burninmat <- test2[['X']]
Startpt <- Burninmat[BurninLength/2,]
rm(Burninmat)  # Housekeeping

## Second burnin
test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
                 H = h, sdB = sdba2, iter = BurninLength/2, type="mirror",
				 jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, x0 = Startpt,
				 fulloutput='TRUE')
				 
Burninmat <- test2[['X']]
Startpt <- Burninmat[BurninLength/2,]
rm(Burninmat)  # Housekeeping

############################################
#### 3b. Run the full model without N15 ####
############################################

test2 <- xsample(A = Aa2, B = ba2, E = Ae, F = be, G = G,
                 H = h, sdB = sdba2, iter = IterLength, outputlength = OutputLength,
				 type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5, x0 = Startpt,
				 fulloutput='TRUE')

MCMCmatplain <- test2[['X']]


##############################
#### 4. N15 model section ####
##############################

# Use full model solution to initialize the N15 matricies 
MCMCmatmean <- colMeans(MCMCmatplain)
RN2 <- 0.0036765;
Eps_Remin <- -1
Eps_Eg <- -2
Eps_NH4up <- -10
Alpha_NH4up <- exp(Eps_NH4up/1000)
Alpha_Eg <- exp(Eps_Eg/1000)
Alpha_Remin <- exp(Eps_Remin/1000)
sdbafactor <- (RN2*Alpha_Remin - RN2)/10

sdba[1] <- sum(MCMCmatmean[c(1,3,10)]) * sdbafactor
sdba[2] <- sum(MCMCmatmean[c(4,11,18,21,25,32)]) * sdbafactor
sdba[3] <- sum(MCMCmatmean[c(2,3,4,5,6,7,8)]) * sdbafactor
sdba[4] <- sum(MCMCmatmean[c(9,10,11,12,13,14,15)]) * sdbafactor
sdba[5] <- sum(MCMCmatmean[c(12,16,17,18,19,20,28)]) * sdbafactor
sdba[6] <- sum(MCMCmatmean[c(5,13,16,21,22,23,24,29)]) * sdbafactor
sdba[7] <- sum(MCMCmatmean[c(6,17,22,25,26,27,30,34)]) * sdbafactor
sdba[8] <- sum(MCMCmatmean[c(7,14,19,23,26,28,29,30,31,33)]) * sdbafactor
sdba[9] <- sum(MCMCmatmean[c(8,15,20,24,27,31,32)]) * sdbafactor
sdba[10] <- sum(MCMCmatmean[c(1,2,9,33,34)]) * sdbafactor


## Burnin - N15 (iterative process)
test2 <- xsampleN15(A = Aa, B = ba, E = Ae, F = be, G = G,
                 H = h, sdB = sdba*10, iter = BurninLength/2,
				 type="mirror", jmp=runif(1, 0, 1)*jmpLength/5,
				 jmp2=jmpLength15/5, x0 = centralval, del15N1=d15N0,
				 del15Nknown=del15Nknown, fulloutput='TRUE')
				 
Burninmat <- test2[['X']]
del15Ntemp <- test2[['del15Ntrack']]
d15N0 <- del15Ntemp[BurninLength/2,]
Startpt <- Burninmat[BurninLength/2,]
rm(del15Ntemp, Burninmat)  # Housekeeping


## Second Burnin - N15
test2 <- xsampleN15(A = Aa, B = ba, E = Ae, F = be, G = G,
                    H = h, sdB = sdba, iter = BurninLength/2,
					type="mirror", jmp=runif(1, 0, 1)*jmpLength/5,
					jmp2=jmpLength15/5, x0 = Startpt, del15N1=d15N0,
					del15Nknown=del15Nknown, fulloutput='TRUE')
					
Burninmat <- test2[['X']]
del15Ntemp <- test2[['del15Ntrack']]
d15N0 <- del15Ntemp[BurninLength/2,]
Startpt <- Burninmat[BurninLength/2,]
rm(del15Ntemp, Burninmat)

################################
#### 5. Full model Run #########
################################

test2 <- xsampleN15(A = Aa, B = ba, E = Ae, F = be, G = G,
                 H = h, sdB = sdba, iter = IterLength, outputlength = OutputLength,
				 type="mirror", jmp=jmpLength+runif(1, 0, 1)*jmpLength/5,
				 jmp2=jmpLength15, x0 = Startpt, del15N1=d15N0, del15Nknown=del15Nknown,
				 fulloutput='TRUE')
				 

MCMCmatN15 <- test2[['X']]
del15N <- test2[['del15Ntrack']]
randomnumber <- test2[['randomnumber']]
acceptedratio <- test2[['acceptedratio']]

##############################
#### 6. Save the output ######
##############################

fileout.mat <- paste('N15InverseModel.',Model,'.ROutputs.mat',sep="")
fileout.r <- paste('N15InverseModel.',Model,'.ROutputs.rdata',sep="")

save(file=fileout.r,
 list(MCMCmatN15 = MCMCmatN15,
 MCMCmatplain = MCMCmatplain,
 del15N = del15N,
 Aa = Aa,
 ba = ba,
 Ae = Ae,
 be = be,
 G = G,
 h = h,
 sdba = sdba,
 lseisol = lseisol,
 Inputs = Inputs,
 InputCol = InputCol,
 jmpLength = jmpLength,
 jmpLength15 = jmpLength15,
 acceptedratio = acceptedratio,
 Startpt = Startpt,
 centralval = centralval,
 KeeperRatio = KeeperRatio,
 Aa2 = Aa2,
 ba2 = ba2,
 sdba2 = sdba2)
 )
 
writeMat(fileout.mat,
 MCMCmatN15 = MCMCmatN15,
 MCMCmatplain = MCMCmatplain,
 del15N = del15N,
 Aa = Aa,
 ba = ba,
 Ae = Ae,
 be = be,
 G = G,
 h = h,
 sdba = sdba,
 lseisol = lseisol,
 Inputs = Inputs,
 InputCol = InputCol,
 jmpLength = jmpLength,
 jmpLength15 = jmpLength15,
 acceptedratio = acceptedratio,
 Startpt = Startpt,
 centralval = centralval,
 KeeperRatio = KeeperRatio,
 Aa2 = Aa2,
 ba2 = ba2,
 sdba2 = sdba2
 )
