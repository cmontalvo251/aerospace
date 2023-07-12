function [BfieldNav,pqrNav,ptpNav] = Navigation(BfieldMeasured,pqrMeasured,ptpMeasured)
global BfieldNavPrev pqrNavPrev ptpNavPrev

s = 0.3;

if BfieldNavPrev(1,1) == -99
   BfieldNav = BfieldMeasured;
   pqrNav = pqrMeasured;
   ptpNav = ptpMeasured;
else
    BiasEstimate = [0;0;0];
    BfieldNav = BfieldNavPrev*(1-s) + s*(BfieldMeasured-BiasEstimate);
    pqrBiasEstimate = [0;0;0];
    pqrNav = pqrNavPrev*(1-s) + s*(pqrMeasured-pqrBiasEstimate);
    ptpBiasEstimate = [0;0;0];
    ptpNav = ptpNavPrev*(1-s) + s*(ptpMeasured-ptpBiasEstimate);
end

BfieldNavPrev = BfieldNav;
pqrNavPrev = pqrNav;
ptpNavPrev = ptpNav;

