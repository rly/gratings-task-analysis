function [z,se] = fisherTransformRToZ(r,n)
z = atanh(r); % (1/2)*(ln(1+r) - ln(1-r));
se = 1/sqrt(n-3);
