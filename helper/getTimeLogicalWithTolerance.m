function timeLog = getTimeLogicalWithTolerance(timeVec, lb, ub)
% buffer the lb and ub of a range with timeTol so that 0.025 >= 0.025 even
% if they are off by miniscule negligible amounts
% lb and ub are inclusive

if nargin == 2 && numel(lb) == 2
    ub = lb(2);
    lb = lb(1);
end

timeTol = 1e-10;

timeLog = timeVec >= lb - timeTol & timeVec <= ub + timeTol;