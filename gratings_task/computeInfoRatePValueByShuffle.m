function infoRateStruct = computeInfoRatePValueByShuffle(...
        timeLockedSpikesStruct, averageFiringRatesBySpdfStruct, analysisWindowOffset)
% skaggs et al 1993 - spatial information rate (bits/sec)
% divide by mean firing rate to get in units of bits/spike
% info rate = sum_over_loc(p(loc) * R(loc) * log2(R(loc) / Rmean))

infoRateStruct.numRandomizations = 2;
maxPValueAdaptive = 0.1;
infoRateStruct.shuffledInfoRates = zeros(infoRateStruct.numRandomizations, 1);
analysisWindowIndices = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, timeLockedSpikesStruct.window(1) + analysisWindowOffset);

infoRateStruct.infoRate = 0;
meanRateAll = averageFiringRatesBySpdfStruct.all;
nLoc = numel(timeLockedSpikesStruct.spikeTimesByLoc);
for i = 1:nLoc
    propTrialsLoc = numel(timeLockedSpikesStruct.spikeTimesByLoc{i}) / numel(timeLockedSpikesStruct.spikeTimes);
    if propTrialsLoc > 0
        infoRateStruct.infoRate = infoRateStruct.infoRate + ...
                propTrialsLoc * averageFiringRatesBySpdfStruct.byLoc(i) * ...
                log2(averageFiringRatesBySpdfStruct.byLoc(i) / meanRateAll);
    end
end

for m = 1:infoRateStruct.numRandomizations
    shuffleIndices = randperm(numel(timeLockedSpikesStruct.spikeTimes));
    
    startTrialCounter = 1;
    endTrialCounter = numel(timeLockedSpikesStruct.spikeTimesByLoc{1});
    for i = 1:nLoc
        propTrialsLoc = numel(timeLockedSpikesStruct.spikeTimesByLoc{i}) / numel(timeLockedSpikesStruct.spikeTimes);
        if propTrialsLoc > 0
            shuffledSpikeTimes = timeLockedSpikesStruct.spikeTimes(shuffleIndices(startTrialCounter:endTrialCounter));
            % if psthT is too large relative to analysisWindowIndices, this
            % could be wasteful
            shuffleSpdfByLoc = fixedPsth(shuffledSpikeTimes, timeLockedSpikesStruct.kernelSigma, 0, timeLockedSpikesStruct.t); 
            meanFiringRateByLoc = mean(shuffleSpdfByLoc(analysisWindowIndices));

            infoRateStruct.shuffledInfoRates(m) = infoRateStruct.shuffledInfoRates(m) + ...
                    propTrialsLoc * meanFiringRateByLoc * ...
                    log2(meanFiringRateByLoc / meanRateAll);
        end
        if i < nLoc
            startTrialCounter = endTrialCounter + 1;
            endTrialCounter = endTrialCounter + numel(timeLockedSpikesStruct.spikeTimesByLoc{i+1});
        end
    end
    
    infoRateStruct.infoRatePValueByShuffleSpdf = sum(infoRateStruct.infoRate < infoRateStruct.shuffledInfoRates) / infoRateStruct.numRandomizations;
    % don't care about actual p-value if > max
    if infoRateStruct.infoRatePValueByShuffleSpdf > maxPValueAdaptive 
        infoRateStruct.infoRatePValueByShuffleSpdf = Inf;
        break;
    end
end
