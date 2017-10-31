%The code communicates with an open excel spreadsheet containing the A,b,G,
%and h matrices for Ax=b & Gx>h.  The new version also allows the option of 
%creating an additional set of matrices Cx~d.  (In this case, the matrices
%A, b, C, d will be named Ae, be, Aa, ba, respectively).  Note that this must 
%be an excel 1997-2003 spreadsheet (not a .xlsx).

clear all
close all

Model='DIAZO.Mesohaline'

% first, get the name of spreadsheet with data
sheetp='N15InverseModelRW';
sheet=[sheetp,'.xls']
workbook='N15'

% now start reading in the data, starting with sizes
datsize = xlsread(sheet,'A2:D2');
neeq = datsize(1);   %Number of exact equalities
naeq = datsize(2);  %Number of approximate equalities
ngt0 = datsize(3);  %Number of inequalities
nvar = datsize(4);  %Number of variables (flows)

ngt = ngt0 + nvar;		% total no. ineqs, including >0
neq = neeq + naeq;      % total number of equalities

rskip = 5;
cskip = 4;
ccl = cskip + 1;	% upper left number to be read
crl = rskip + 1;
crr = rskip + nvar + 1;
ccr = cskip + neq + ngt0 + 2;
crr2 = rskip + nvar + 2;
cread = [char(ExcelCol(ccl)), num2str(crl), ':', char(ExcelCol(ccr)), num2str(crr)]
M0 = xlsread(sheet, workbook, cread);

% % now, sort out
% % do some minimal prcocessing
weight = M0(1:nvar, 1);
M = [M0(:,2:(neq+ngt0+1)),eye(nvar+1 nvar)];

% %setup the weighting
% Awt=diag(1. ./ weight);
% Araw=M(1:nvar,1:neq)';
% A=Araw;    %These steps may have to be changed if we actually use weights
A = M(1:nvar, 1:neq)';
Ae = A(1:neeq, :);
Aa = A(neeq+1:neeq+naeq, :);

% Graw=M(1:nvar,neq+1:neq+ngt)';
% G=Graw;    %These steps may have to be changed if we actually use weights
G = M(1:nvar, neq+1:neq+ngt)';

b = M(nvar+1,1:neq)';

Inputs = xlsread('../../Inputs2.xls','Sheet1','J4:N26');
if strcmp(Model,'DIAZO.Coastal')
    col = 1;
elseif strcmp(Model,'DIAZO.Mesohaline')
    col = 2;
elseif strcmp(Model,'NEMURO.Coastal')
    col = 3;
elseif strcmp(Model,'NEMURO.Offshore')
    col = 4;
elseif strcmp(Model,'Base')
    col = 5;
end
      
b(end) = Inputs(6,col);    %SinkingFlux
b(end-1) = Inputs(2,col);  %New Production (Nitrate uptake - Nitrification)
b(end-2) = Inputs(5,col);    %Mesozoo Grazing
b(end-3) = Inputs(16,col) + Inputs(17,col); %NPP

be = b(1:neeq);
ba = b(neeq+1:neeq+naeq);
h = M(nvar+1, neq+1:neq+ngt)';
Temp = Inputs(14, col);
Weight = 7.5;
h(3) = -1.7 * Weight.^-0.25 * exp(0.0693 * (Temp-20)) * Inputs(10,col);   %Protistan Max Respiration
Weight = 3800000;
h(5) = -14 * Weight.^-0.25 * exp(0.0693 * (Temp-20)) * Inputs(11,col);    %Mesozoo Max Respiration
h(6) = 0.02 * Inputs(13, col);                                 %Diatom Min Excretion
h(7) = -0.55 * Inputs(13, col);                                 %Diatom Max Excretion
h(8) = 0.02 * Inputs(12, col);                                 %Cyano Min Excretion
h(9) = -0.55 * Inputs(12, col);                                %Cyano Max Excretion
h(16) = -3.6 * Weight.^-0.25 * exp(0.0693 * (Temp-20)) * Inputs(11,col);  %Mesozoo Max Ingestion

RN2 = 0.0036765;
Eps_TL = 3.5;
R_TL = Eps_TL/1000 * RN2 + RN2;

ToNO3 = 1;
FromNH4 = [4, 11];
ToCya = [9, 10, 11];
ToDTM = [2, 3, 4];
ToHNF = [12, 28];
ToMIC = [5, 13, 16, 29];
ToMES = [6, 17, 22, 30, 35];
FromDet = [28:31];
FromDON = 32;

sdba = ba/10;
sdba(1:9) = (R_TL-RN2)/10;
sdba(10) = Inputs(2,col) * (R_TL-RN2)/10;  %d15NExportNPBalance
if exist(['N15InverseModel.', Model, '.ROutputs.mat'])>0
    load(['N15InverseModel.', Model, '.ROutputs.mat'], 'MCMCmatplain')
    tmp = median(MCMCmatplain);
    sdba(1) = sum(tmp(ToNO3)) * (R_TL-RN2) / 10;   %d15NNO3
    sdba(2) = sum(tmp(FromNH4)) * (R_TL-RN2) / 10;   %d15NNH4
    sdba(3) = sum(tmp(ToCya)) * (R_TL-RN2) / 10;   %d15NCya
    sdba(4) = sum(tmp(ToDTM)) * (R_TL-RN2) / 10;   %d15NDtm
    sdba(5) = sum(tmp(ToHNF)) * (R_TL-RN2) / 10;   %d15NHNF
    sdba(6) = sum(tmp(ToMIC)) * (R_TL-RN2) / 10;   %d15NMic
    sdba(7) = sum(tmp(ToMES)) * (R_TL-RN2) / 10;   %d15NMes
    sdba(8) = sum(tmp(FromDet)) * (R_TL-RN2) / 10;   %d15NDet
    sdba(9) = sum(tmp(FromDON)) * (R_TL-RN2) / 10;   %d15NDON
    sdba(10) = sum(tmp(ToNO3)) * (R_TL-RN2) / 10;  %d15NExportNPBalance
    sdbainputs = 1
end

% %Just clearing variables that we don't need to see
% clear Araw
% clear Awt
% clear Graw
% clear M
% clear M0
% clear captions
% clear cc2
% clear cc3
% clear ccl
% clear ccr
% clear chsheet
% clear comment
% clear constraintout
% clear constraints
% clear cread
% clear cread2
% clear cread3
% clear cread4
% clear crl
% clear crr
% clear crr2
% clear crr3
% clear crr4
% clear cskip
% clear datsize
% clear equalities
% clear equalitiesout
% clear linsum
% clear neq
% clear naeq
% clear neeq
% clear ngt
% clear ngt0
% clear norm2
% clear npar
% clear nvar
% clear outgate
% clear rc1
% clear cnorm
% clear rr
% clear rrout
% clear rrtrans
% clear rskip
% clear sheet
% clear sols
% clear sqrterror
% clear temp
% clear time1
% clear weight
% clear SD


InputCol = col;

if length(ba)>0
    save([sheetp,'.',Model,'.RW.mat'],'A','Ae','Aa','G','b','be','ba','h','Inputs','sdba','InputCol')
else
    save([sheetp,'.',Model,'.RW.mat'],'A','G','b','h')
end