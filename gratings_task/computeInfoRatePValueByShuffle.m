function infoRateStruct = computeInfoRatePValueByShuffle(...
        eventStruct, averageFiringRatesBySpdfStruct, analysisWindowOffset, numRandomizations)
% skaggs et al 1993 - spatial information rate (bits/sec)
% divide by mean firing rate to get in units of bits/spike
% info rate = sum_over_loc(p(loc) * R(loc) * log2(R(loc) / Rmean))

maxPValueAdaptive = 0.1;
analysisWindowIndices = getTimeLogicalWithTolerance(eventStruct.t, eventStruct.window(1) + analysisWindowOffset);

infoRate = 0;
meanRateAll = averageFiringRatesBySpdfStruct.all;
nLoc = numel(eventStruct.spikeTimesByLoc);
for i = 1:nLoc
    propTrialsLoc = numel(eventStruct.spikeTimesByLoc{i}) / numel(eventStruct.spikeTimes);
    if propTrialsLoc > 0
        infoRate = infoRate + ...
                propTrialsLoc * averageFiringRatesBySpdfStruct.byLoc(i) * ...
                log2(averageFiringRatesBySpdfStruct.byLoc(i) / meanRateAll);
    end
end

shuffledInfoRates = zeros(numRandomizations, 1);
for m = 1:numRandomizations
    % randomly reassign cue location condition to different trials, maintaining
    % original number of trials per condition
    shuffleIndices = randperm(numel(eventStruct.spikeTimes));
    
    startTrialCounter = 1;
    endTrialCounter = numel(eventStruct.spikeTimesByLoc{1});
    for i = 1:nLoc
        propTrialsLoc = numel(eventStruct.spikeTimesByLoc{i}) / numel(eventStruct.spikeTimes);
        if propTrialsLoc > 0
            shuffledSpikeTimes = eventStruct.spikeTimes(shuffleIndices(startTrialCounter:endTrialCounter));
            % if psthT is too large relative to analysisWindowIndices, this
            % could be wasteful
            shuffleSpdfByLoc = fixedPsth(shuffledSpikeTimes, eventStruct.kernelSigma, 0, eventStruct.t); 
            meanFiringRateByLoc = mean(shuffleSpdfByLoc(analysisWindowIndices));

            shuffledInfoRates(m) = shuffledInfoRates(m) + ...
                    propTrialsLoc * meanFiringRateByLoc * ...
                    log2(meanFiringRateByLoc / meanRateAll);
        end
        if i < nLoc
            startTrialCounter = endTrialCounter + 1;
            endTrialCounter = endTrialCounter + numel(eventStruct.spikeTimesByLoc{i+1});
        end
    end
    
    pValue = sum(infoRate < shuffledInfoRates) / numRandomizations;
    % don't care about actual p-value if > max. short circuit loop
    if pValue > maxPValueAdaptive 
        pValue = Inf;
        break;
    end
end
infoRateStruct.infoRate = infoRate;
infoRateStruct.permutation.shuffledInfoRates = shuffledInfoRates;
infoRateStruct.permutation.numRandomizations = numRandomizations;
infoRateStruct.permutation.p = pValue;
