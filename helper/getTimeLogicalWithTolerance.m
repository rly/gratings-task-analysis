function isInBounds = getTimeLogicalWithTolerance(t, lb, ub)
% returns a logical vector for which values in t are between lb and ub with
% some tolerance. this ensures 0.025 >= 0.025 even if they are off by 
% miniscule negligible amounts. lb and ub are inclusive
% 
% inputs:
% - t: a vector of time points
% - lb: lower bound, inclusive
% - ub: upper bound, inclusive
%
% OR if second element is a vector of two elements, interpret it as [lb ub]
%
% output:
% - isInBounds: logical (boolean) vector of which elements in t are between
%               lb and ub, inclusive, with some tolerance

if nargin == 2 && numel(lb) == 2
    ub = lb(2);
    lb = lb(1);
end

timeTol = 1e-10;

isInBounds = t >= lb - timeTol & t <= ub + timeTol;