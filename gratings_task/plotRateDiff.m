function [f1, ax1, ax2, ax3] = plotRateDiff(firingRates, inRFLocs, exRFLocs, isDPul, isVPul, isSigUnit)

nUnits = numel(firingRates);
firingInRF = nan(nUnits, 1);
firingExRF = nan(nUnits, 1);
for i = 1:nUnits
    firingInRF(i) = firingRates(i).byLoc(inRFLocs(i));
    firingExRF(i) = firingRates(i).byLoc(exRFLocs(i));
end

fprintf('\tMean firing rate InRF: %0.1f Hz, ExRF: %0.1f Hz\n', mean(firingInRF), mean(firingExRF));

[f1, ax1, ax2, ax3] = plotMetricDiff(firingInRF, firingExRF, isDPul, isVPul, isSigUnit, 1);
