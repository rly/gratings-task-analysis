function [fanoFactorAI,toRemove] = computeFanoFactorAIDiss(firingRates, inRFLocs, exRFLocs)

nUnits = numel(firingRates);
fanoFactorInRF = nan(nUnits, 1);
fanoFactorExRF = nan(nUnits, 1);
for i = 1:nUnits
    fanoFactorInRF(i) = computeFanoFactor(firingRates(i), inRFLocs(i));
    fanoFactorExRF(i) = computeFanoFactor(firingRates(i), exRFLocs(i));
end

maxFanoFactor = 5;
toRemove = fanoFactorInRF > maxFanoFactor | fanoFactorExRF > maxFanoFactor;
fanoFactorInRF(toRemove) = NaN;
fanoFactorExRF(toRemove) = NaN;

fanoFactorAI = (fanoFactorInRF - fanoFactorExRF) ./ (fanoFactorInRF + fanoFactorExRF);
