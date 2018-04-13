function attnStats = permutationTestAttn(eventStruct, eventWindowOffset, inRFLoc, exRFLoc, actualRatesByLoc, kernelSigma, numRandomizations)

%% permutation test on cue-target delay period
% randomly reassign InRF / ExRF condition to different trials, maintaining
% original number of trials per condition
shuffleDiff = zeros(numRandomizations, 1);
shuffleAI = zeros(numRandomizations, 1);
windowLogical = getTimeLogicalWithTolerance(eventStruct.t, eventStruct.window(1) + eventWindowOffset);
nTrialsInRF = numel(eventStruct.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(eventStruct.spikeTimesByLoc{exRFLoc});
nTrialsAll = nTrialsInRF + nTrialsExRF;
spikeTimesAll = [eventStruct.spikeTimesByLoc{inRFLoc} eventStruct.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations
    shuffleIndices = randperm(nTrialsAll);
    randSpikeTimesInRF = spikeTimesAll(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = spikeTimesAll(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, eventStruct.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, eventStruct.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(windowLogical));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(windowLogical));
    shuffleDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end

actualDiff = actualRatesByLoc(inRFLoc) - actualRatesByLoc(exRFLoc);
actualSum = actualRatesByLoc(inRFLoc) + actualRatesByLoc(exRFLoc);
actualAI = actualDiff / actualSum;

if actualDiff > median(shuffleDiff)
    attnStats.diff.permutation.p = sum(actualDiff < shuffleDiff) / numRandomizations * 2;
else
    attnStats.diff.permutation.p = sum(actualDiff > shuffleDiff) / numRandomizations * 2;
end
attnStats.diff.actual = actualDiff;
attnStats.diff.permutation.numRandomizations = numRandomizations;

if actualAI > median(shuffleAI)
    attnStats.ai.permutation.p = sum(actualAI < shuffleAI) / numRandomizations * 2;
else
    attnStats.ai.permutation.p = sum(actualAI > shuffleAI) / numRandomizations * 2;
end
attnStats.ai.actual = actualAI;
attnStats.ai.permutation.numRandomizations = numRandomizations;
