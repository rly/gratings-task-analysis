function x = cleanBySDsFromMean(x, nSDs)
% replaces values that are more than nSDs away from the mean of x with NaN
% if x is not 1-D, collapse x
xall = x(:);
m = mean(xall);
sd = std(xall);
ub = m + nSDs * sd;
lb = m - nSDs * sd;

x(x > ub | x < lb) = NaN;