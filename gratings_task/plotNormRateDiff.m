function [f1, ax1, ax2, ax3] = plotNormRateDiff(firingRates, baselineFiringRates, ...
        inRFNormFactor, exRFNormFactor, inRFLocs, exRFLocs, isDPul, isVPul, isSigUnit)

nUnits = numel(firingRates);
firingInRF = nan(nUnits, 1);
firingExRF = nan(nUnits, 1);
for i = 1:nUnits
    firingInRF(i) = (firingRates(i).byLoc(inRFLocs(i)) - baselineFiringRates(i).byLoc(inRFLocs(i))) / inRFNormFactor(i);
    firingExRF(i) = (firingRates(i).byLoc(exRFLocs(i)) - baselineFiringRates(i).byLoc(exRFLocs(i))) / exRFNormFactor(i);
end

fprintf('\tMean normalized firing rate InRF: %0.1f, ExRF: %0.1f\n', mean(firingInRF), mean(firingExRF));

[f1, ax1, ax2, ax3] = plotMetricDiff(firingInRF, firingExRF, isDPul, isVPul, isSigUnit, 0.1);
