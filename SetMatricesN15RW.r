
# The code communicates with an open excel spreadsheet containing the A,b,G,
# and h matrices for Ax=b & Gx>h.  The new version also allows the option of 
# creating an additional set of matrices Cx~d.  (In this case, the matrices
# A, b, C, d will be named Ae, be, Aa, ba, respectively).  Note that this must 
# be an excel 1997-2003 spreadsheet (not a .xlsx).
library(openxlsx)

model = 'DIAZO.Mesohaline'

# first, get the name of spreadsheet with data
model.file = 'N15InverseModelRW'
sheet.name = 'N15'
rskip = 5
cskip = 4


## Helper function to convert data.frames and others into matrixes.
to.matrix = function(x) {
    m = matrix(as.numeric(unlist(x)), nrow=nrow(x))
    m[is.na(m)] = 0  ## NB: may be dangerous if row and col numbers of inaccurate (no signs of failure)
    m
}

## Read in the excel file (model)
sheet = read.xlsx(paste0(model.file,'.xlsx'), sheet = sheet.name)

## Now start reading in the data, starting with sizes
datsize = as.numeric(sheet[1,1:4])
neeq = datsize[1]   #Number of exact equalities
naeq = datsize[2]   #Number of approximate equalities
ngt0 = datsize[3]   #Number of inequalities
nvar = datsize[4]   #Number of variables (flows)

ngt = ngt0 + nvar   # total no. ineqs, including >0
neq = neeq + naeq      # total number of equalities


## Determine spreadsheet bounds 
ccl = cskip + 1	# upper left number to be read
crl = rskip + 1
crr = rskip + nvar + 1
ccr = cskip + neq + ngt0 + 2
crr2 = rskip + nvar + 2

## Master sheet
M0 = sheet[c(crl:crr), c(ccl:ccr)]

## now, sort out
## do some minimal prcocessing
weight = M0[1:nvar, 1]
M = cbind(M0[, 2:(neq+ngt0+1)], diag(nvar+1))


#### Set 'A' matrix

# A.wt = diag(1 / weight)
# A.raw = to.matrix(M[1:nvar, 1:neq])
# A = A.raw %*% A.wt    ## NB: These steps may have to be varified before use
A = to.matrix( M[1:nvar, 1:neq] )

Ae = A[1:neeq, ] # exact
Aa = A[(neeq+1):(neeq+naeq),] # approx

#### Setup 'G' matrix

# G.raw = to.matrix(M[1:nvar, (neq+1):(neq+ngt)]) %*% A.wt
# G = G.raw    ## NB: These steps may have to be varified before use
G = M[1:nvar, (neq+1):(neq+ngt)]
G = to.matrix(G)

Inputs = read.xlsx('Inputs.xlsx', 'Sheet1')
Inputs = to.matrix(Inputs[4:26, 10:16])  ## Determined from Rows and Columns of the spreadsheet (i.e. N4:T26)

## Each Model is a different col.
if (model == 'DIAZO.Coastal') {
    col = 1
}
if (model == 'DIAZO.Mesohaline') {
    col = 2
}
if (model == 'NEMURO.Coastal') {
    col = 3
}
if (model == 'NEMURO.Offshore') {
    col = 4
}
if (model == 'Base') {
    col = 5
}

#### Setup 'b' vectors ####

b = to.matrix(M[nvar+1, 1:neq])

b.end = length(b)
b[b.end] = as.numeric(Inputs[6,col])      # SinkingFlux
b[b.end - 1] = as.numeric(Inputs[2,col])  # New Production (Nitrate uptake - Nitrification)
b[b.end - 2] = as.numeric(Inputs[5,col])  # Mesozoo Grazing
b[b.end - 3] = as.numeric(Inputs[16,col]) + as.numeric(Inputs[17,col]) # NPP

## Set ba & be from b
be = b[1:neeq]
ba = b[(neeq+1):(neeq+naeq)]


#### Setup 'h' vector ####

h = to.matrix(M[nvar+1, (neq+1):(neq+ngt)])

temp = as.numeric(Inputs[14, col])

Weight = 7.5;
h[3] = -1.7 * Weight^-0.25 * exp(0.0693 * (temp-20)) * as.numeric(Inputs[10, col])   #Protistan Max Respiration

Weight = 3800000
h[5] = -14 * Weight^-0.25 * exp(0.0693 * (temp-20)) * as.numeric(Inputs[11, col])    #Mesozoo Max Respiration
h[6] = 0.02 * as.numeric(Inputs[13, col])                                 #Diatom Min Excretion
h[7] = -0.55 * as.numeric(Inputs[13, col])                                #Diatom Max Excretion
h[8] = 0.02 * as.numeric(Inputs[12, col])                                 #Cyano Min Excretion
h[9] = -0.55 * as.numeric(Inputs[12, col])                                #Cyano Max Excretion
h[16] = -3.6 * Weight^-0.25 * exp(0.0693 * (temp-20)) * as.numeric(Inputs[11, col])  #Mesozoo Max Ingestion

## Setup fractionation parameters
RN2 = 0.0036765
Eps_TL = 3.5
R_TL = Eps_TL / 1000 * RN2 + RN2

## Determine flows for N15 calculations

ToNO3 = 1
FromNH4 = c(4, 11)
ToCya = c(9, 10, 11)
ToDTM = c(2, 3, 4)
ToHNF = c(12, 28)
ToMIC = c(5, 13, 16, 29)
ToMES = c(6, 17, 22, 30, 35)
FromDet = c(28:31)
FromDON = 32

#### Setup sdba vector ####

sdba = ba / 10  ## Default to 10% uncertainty
sdba[1:9] = (R_TL - RN2) / 10  ## 10% uncertainty on N15 equations
sdba[10] = as.numeric(Inputs[2,col]) * (R_TL - RN2) / 10  ## d15NExportNPBalance


## Load N15 flow data from previous model to establish approximate equations
if (file.exists(paste0('N15InverseModel.', model, '.ROutputs.rdata'))) {
    
    load(paste0('N15InverseModel.', model, '.ROutputs.rdata')) ## MCMCmatplain is loaded
    
    tmp = apply(MCMCmatplain, 2, median)
    
    sdba[1] = sum(tmp[ToNO3]) * (R_TL - RN2) / 10   # d15NNO3
    sdba[2] = sum(tmp[FromNH4]) * (R_TL - RN2) / 10   # d15NNH4
    sdba[3] = sum(tmp[ToCya]) * (R_TL - RN2) / 10   # d15NCya
    sdba[4] = sum(tmp[ToDTM]) * (R_TL - RN2) / 10   # d15NDtm
    sdba[5] = sum(tmp[ToHNF]) * (R_TL - RN2) / 10   # d15NHNF
    sdba[6] = sum(tmp[ToMIC]) * (R_TL - RN2) / 10   # d15NMic
    sdba[7] = sum(tmp[ToMES]) * (R_TL - RN2) / 10   # d15NMes
    sdba[8] = sum(tmp[FromDet]) * (R_TL - RN2) / 10   # d15NDet
    sdba[9] = sum(tmp[FromDON]) * (R_TL - RN2) / 10   # d15NDON
    sdba[10] = sum(tmp[ToNO3]) * (R_TL - RN2) / 10  # d15NExportNPBalance
    sdbainputs = 1
}




#### Save output

if (length(ba) > 0) {
    output = list(A=A, Ae=Ae, Aa=Aa, G=G, b=b, be=be, ba=ba, h=h, Inputs=Inputs, sdba=sdba, InputCol=col)
    
} else {
    output = list(A=A, G=G, b=b, h=h)
}
save(file = paste0(model.file, '.', model, '.RW.rdata'), output)
