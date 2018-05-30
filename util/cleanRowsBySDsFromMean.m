function [x,keep] = cleanRowsBySDsFromMean(x, nSDs)
% removes rows with values more than nSDs away from the mean of x (2D)

m = mean(x, 1);
sd = std(x, 0, 1);
ub = m + nSDs * sd;
lb = m - nSDs * sd;

keep = true(size(x, 1), 1);
for i = 1:numel(m)
    keep = keep & x(:,i) >= lb(i) & x(:,i) <= ub(i);
end

x(~keep,:) = [];