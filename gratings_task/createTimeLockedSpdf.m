function timeLockStruct = createTimeLockedSpdf(spikeTs, eventTimes, eventTimesByLoc, timeLockStruct, kernelSigma)
% timeLockStruct.window and timeLockStruct.spdfWindowOffset must exist
nLoc = numel(eventTimesByLoc);

timeLockStruct.kernelSigma = kernelSigma;
timeLockStruct.t = computeTForSpdf(timeLockStruct.window(1), timeLockStruct.spdfWindowOffset, kernelSigma);

timeLockStruct.spikeTimes = createnonemptydatamatpt(spikeTs, eventTimes, timeLockStruct.window);

[timeLockStruct.spdf,~,timeLockStruct.spdfErr] = fixedPsth(timeLockStruct.spikeTimes, kernelSigma, 2, timeLockStruct.t);

timeLockStruct.spikeTimesByLoc = cell(nLoc, 1);
timeLockStruct.spdfByLoc = nan(nLoc, numel(timeLockStruct.t));
timeLockStruct.spdfErrByLoc = nan(nLoc, numel(timeLockStruct.t));
for i = 1:nLoc
    timeLockStruct.spikeTimesByLoc{i} = createnonemptydatamatpt(spikeTs, eventTimesByLoc{i}, timeLockStruct.window);
    [timeLockStruct.spdfByLoc(i,:),~,timeLockStruct.spdfErrByLoc(i,:)] = fixedPsth(timeLockStruct.spikeTimesByLoc{i}, kernelSigma, 2, timeLockStruct.t);
end

end