function statsByLoc = getNullDistPValByLoc(actualByLoc, nullDist, statsByLoc)
% do two-sided test of actual vs null distribution, split by location
if nargin < 3
    statsByLoc = struct();
end

numRandomizations = numel(nullDist);
nLoc = numel(actualByLoc);
medianNullDist = median(nullDist);

for i = 1:nLoc
    if actualByLoc(i) > medianNullDist
        statsByLoc(i).bootstrap.p = sum(actualByLoc(i) < nullDist) / numRandomizations * 2;
    else
        statsByLoc(i).bootstrap.p = sum(actualByLoc(i) > nullDist) / numRandomizations * 2;
    end
    statsByLoc(i).bootstrap.numRandomizations = numRandomizations;
end