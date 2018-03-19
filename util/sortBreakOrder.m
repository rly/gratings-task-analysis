function [sx,sortInd] = sortBreakOrder(x)
% inefficient sorting algorithm that randomly shuffles the order of equal
% elements (usually order is preserved)

assert(isrow(x) | iscolumn(x));

sortInd = nan(size(x));

ux = unique(sort(x));
count = 1;
for i = 1:numel(ux)
    matches = find(x == ux(i));
    perm = randperm(numel(matches));
    newInd = count:count + numel(matches) - 1;
    sortInd(newInd) = matches(perm);
    count = count + numel(matches);
end

sx = x(sortInd);
