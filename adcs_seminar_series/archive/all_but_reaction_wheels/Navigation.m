function [BfieldNav,pqrNav] = Navigation(BfieldMeasured,pqrMeasured)
global BfieldNavPrev pqrNavPrev

s = 0.3;

if sum(BfieldNavPrev) + sum(pqrNavPrev) == 0
   BfieldNav = BfieldMeasured;
   pqrNav = pqrMeasured;
else
    BiasEstimate = [0;0;0];
    BfieldNav = BfieldNavPrev*(1-s) + s*(BfieldMeasured-BiasEstimate);
    pqrBiasEstimate = [0;0;0];
    pqrNav = pqrNavPrev*(1-s) + s*(pqrMeasured-pqrBiasEstimate);
end

BfieldNavPrev = BfieldNav;
pqrNavPrev = pqrNav;

