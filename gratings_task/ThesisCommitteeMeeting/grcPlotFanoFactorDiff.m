function [f1,ax1] = grcPlotFanoFactorDiff(firingRates, inRFLocs, exRFLocs, isDPul, isVPul, isSigUnit)

nUnits = numel(firingRates);
fanoFactorInRF = nan(nUnits, 1);
fanoFactorExRF = nan(nUnits, 1);
for i = 1:nUnits
    fanoFactorInRF(i) = computeFanoFactor(firingRates(i), inRFLocs(i));
    fanoFactorExRF(i) = computeFanoFactor(firingRates(i), exRFLocs(i));
end

fprintf('\tMean Fano factor InRF: %0.2f, ExRF %0.2f\n', mean(fanoFactorInRF), mean(fanoFactorExRF));

maxFanoFactor = 5;
toRemove = fanoFactorInRF > maxFanoFactor | fanoFactorExRF > maxFanoFactor;
fanoFactorInRF(toRemove) = [];
fanoFactorExRF(toRemove) = [];
isDPul(toRemove) = [];
isVPul(toRemove) = [];
isSigUnit(toRemove) = [];
fprintf('\tRemoved %d units because of outlier Fano factor (> %0.1f)\n', sum(toRemove), maxFanoFactor);

[f1,ax1] = grcPlotMetricDiff(fanoFactorInRF, fanoFactorExRF, isDPul, isVPul, isSigUnit);
