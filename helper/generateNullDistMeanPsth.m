function nullDist = generateNullDistMeanPsth(eventStruct, eventWindowOffset, kernelSigma, numRandomizations)
% generate a null distribution of mean psth values in the given
% eventWindowOffset around the event described by eventStruct

nTrials = numel(eventStruct.spikeTimes);
windowLogical = getTimeLogicalWithTolerance(eventStruct.t, eventStruct.window(1) + eventWindowOffset);

nullDist = zeros(numRandomizations, 1);
for m = 1:numRandomizations
    bootRandInd = randi(nTrials, 1, nTrials); % sample with replacement
    bootstrapPsth = fixedPsth(eventStruct.spikeTimes(bootRandInd), kernelSigma, 0, eventStruct.t);
    if ~isempty(bootstrapPsth)
        nullDist(m) = mean(bootstrapPsth(windowLogical));
    end % else, there were no spikes and leave nullDist as zeros
end

