function timeLockStruct = createTimeLockedSpdf(spikeTs, eventTimes, eventTimesByLoc, timeLockStruct, kernelSigma, startTime, endTime)
% timeLockStruct.window and timeLockStruct.spdfWindowOffset must exist
nLoc = numel(eventTimesByLoc);

spikeTs = spikeTs(spikeTs >= startTime & spikeTs <= endTime);

timeLockStruct.kernelSigma = kernelSigma;
timeLockStruct.t = computeTForSpdf(timeLockStruct.window(1), timeLockStruct.spdfWindowOffset, kernelSigma);

if ~isempty(eventTimes)
    eventTimesForBoolean = eventTimes;
    eventTimes = eventTimes(((eventTimes > startTime) & (eventTimes <= endTime)));
    for i = 1:nLoc
        timeLockStruct.trialIndicesByLoc{i} = ((eventTimesByLoc{i} > startTime) & (eventTimesByLoc{i} <= endTime));
        eventTimesByLoc{i} = eventTimesByLoc{i}(timeLockStruct.trialIndicesByLoc{i});
    end
    timeLockStruct.trialsBetweenStartEndTimes = ((eventTimesForBoolean > startTime) & (eventTimesForBoolean <= endTime));
end
    
[timeLockStruct.spikeTimes,timeLockStruct.spikeIndices] = createnonemptydatamatpt(spikeTs, eventTimes, timeLockStruct.window);

spikeCount = zeros(size(timeLockStruct.spikeTimes,2),1);
for i = 1:size(timeLockStruct.spikeTimes,2)
    spikeCount(i) = length(timeLockStruct.spikeTimes(i).times);
end
timeLockStruct.spikeCount = spikeCount;

meanSpikeCount = mean(spikeCount);
sdSpikeCount = std(spikeCount);
indices = 1:10:length(spikeCount);
for i = 1:length(indices)
    if i == length(indices)
        meanSpikeCountPerTen(i,:) = mean(spikeCount(indices(i):indices(end)));
        sdSpikeCountPerTen(i,:) = std(spikeCount(indices(i):indices(end)));
    else
        meanSpikeCountPerTen(i,:) = mean(spikeCount(indices(i):indices(i+1)));
        sdSpikeCountPerTen(i,:) = std(spikeCount(indices(i):indices(i+1)));
    end
end

stdFactor = 2;
% spikeCountInd = boolean(spikeCount > mean(spikeCount(spikeCount > 0)) - stdFactor * std(spikeCount(spikeCount > 0)) / sqrt(length(spikeCount)) );
spikeCountInd = boolean(spikeCount > mean(spikeCount) - stdFactor * std(spikeCount) );
timeLockStruct.rateCorrInd = spikeCountInd;

[timeLockStruct.spdf,~,timeLockStruct.spdfErr] = fixedPsth(timeLockStruct.spikeTimes, kernelSigma, 2, timeLockStruct.t);
timeLockStruct.spdfRateCorr = fixedPsth(timeLockStruct.spikeTimes(spikeCountInd), kernelSigma, 2, timeLockStruct.t);

timeLockStruct.spikeTimesByLoc = cell(nLoc, 1);
timeLockStruct.spdfByLoc = nan(nLoc, numel(timeLockStruct.t));
timeLockStruct.spdfErrByLoc = nan(nLoc, numel(timeLockStruct.t));
for i = 1:nLoc
    timeLockStruct.spikeTimesByLoc{i} = createnonemptydatamatpt(spikeTs, eventTimesByLoc{i}, timeLockStruct.window);
    spikeCountByLoc = zeros(size(timeLockStruct.spikeTimesByLoc{i},2),1);
    for ii = 1:size(timeLockStruct.spikeTimesByLoc{i},2)
        spikeCountByLoc(ii) = length(timeLockStruct.spikeTimesByLoc{i}(ii).times);
    end
%     spikeCountIndByLoc = boolean(spikeCountByLoc > mean(spikeCount(spikeCount > 0)) - stdFactor * std(spikeCount(spikeCount > 0)) / sqrt(length(spikeCount)) );
    spikeCountIndByLoc = boolean(spikeCountByLoc > mean(spikeCount(spikeCount > 0)) - stdFactor * std(spikeCount(spikeCount > 0)) );
    clear spikeCountByLoc
    timeLockStruct.rateCorrIndByLoc{i} = spikeCountIndByLoc;

    [timeLockStruct.spdfByLoc(i,:),~,timeLockStruct.spdfErrByLoc(i,:)] = fixedPsth(timeLockStruct.spikeTimesByLoc{i}, kernelSigma, 2, timeLockStruct.t);
    timeLockStruct.spdfRateCorrByLoc(i,:) = fixedPsth(timeLockStruct.spikeTimesByLoc{i}(spikeCountIndByLoc), kernelSigma, 2, timeLockStruct.t);
end

end