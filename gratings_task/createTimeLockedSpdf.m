function timeLockStruct = createTimeLockedSpdf(spikeTs, eventTimes, eventTimesByLoc, timeLockStruct, kernelSigma, ...
        startTime, endTime)
% Create a struct of spike density functions aligned to event times
% Parameters:
% - spikeTs - spike times, in seconds from recording file (0 = start of
% recording)
% - eventTimes - times of the events of interest, in seconds from recording
% file
% - eventTimesByLoc - same as above, but split by cue location
% - timeLockStruct - struct containing:
% --- window and
% --- spdfWindowOffset which specify the window for cutting the data 
% (seconds before, seconds after) and the offset window for computing the 
% spdf
% - kernelSigma - the sigma for the Gaussian kernel for convolving around
% spikes
% - startTime - lower bound for time points to use
% - endTime - upper bound for time points to use
%
% Returns:
%   timeLockStruct with many useful parameters

% TODO was this necessary?
% spikeTs = spikeTs(spikeTs >= startTime & spikeTs <= endTime);

nLoc = numel(eventTimesByLoc);
timeLockStruct.kernelSigma = kernelSigma;

% compute time vector at which spdf values will be computed
timeLockStruct.t = computeTForSpdf(timeLockStruct.window(1), timeLockStruct.spdfWindowOffset, kernelSigma);

% filter event times based on startTime and endTime == "valid"
if ~isempty(eventTimes)
    timeLockStruct.validEvents = (eventTimes > startTime) & (eventTimes <= endTime);
    timeLockStruct.validEventTimes = eventTimes(timeLockStruct.validEvents);
    for i = 1:nLoc
        timeLockStruct.validEventsByLoc{i} = (eventTimesByLoc{i} > startTime) & (eventTimesByLoc{i} <= endTime);
        timeLockStruct.validEventTimesByLoc{i} = eventTimesByLoc{i}(timeLockStruct.validEventsByLoc{i});
    end
end

% align spike times to event (0 = start of window)
[timeLockStruct.spikeTimes,timeLockStruct.spikeIndices] = createnonemptydatamatpt(spikeTs, timeLockStruct.validEventTimes, timeLockStruct.window);

spikeCount = zeros(size(timeLockStruct.spikeTimes,2),1);
for i = 1:size(timeLockStruct.spikeTimes,2)
    spikeCount(i) = length(timeLockStruct.spikeTimes(i).times);
end
timeLockStruct.spikeCount = spikeCount;

% meanSpikeCount = mean(spikeCount);
% sdSpikeCount = std(spikeCount);
% indices = 1:10:length(spikeCount);
% for i = 1:length(indices)
%     if i == length(indices)
%         meanSpikeCountPerTen(i,:) = mean(spikeCount(indices(i):indices(end)));
%         sdSpikeCountPerTen(i,:) = std(spikeCount(indices(i):indices(end)));
%     else
%         meanSpikeCountPerTen(i,:) = mean(spikeCount(indices(i):indices(i+1)));
%         sdSpikeCountPerTen(i,:) = std(spikeCount(indices(i):indices(i+1)));
%     end
% end

stdFactor = 2;
% spikeCountInd = boolean(spikeCount > mean(spikeCount(spikeCount > 0)) - stdFactor * std(spikeCount(spikeCount > 0)) / sqrt(length(spikeCount)) );
spikeCountInd = boolean(spikeCount > mean(spikeCount) - stdFactor * std(spikeCount) );
timeLockStruct.rateCorrInd = spikeCountInd;

% compute the spike density function and bootstrapped error measure
[timeLockStruct.spdf,~,timeLockStruct.spdfErr] = fixedPsth(timeLockStruct.spikeTimes, kernelSigma, 2, timeLockStruct.t);
timeLockStruct.spdfRateCorr = fixedPsth(timeLockStruct.spikeTimes(spikeCountInd), kernelSigma, 2, timeLockStruct.t);

% compute spdf for single trials
timeLockStruct.spdfByEvent = nan(numel(timeLockStruct.validEventTimes), numel(timeLockStruct.t));
for i = 1:numel(timeLockStruct.validEventTimes)
    alignedSpikeTimesThisEvent = createnonemptydatamatpt(spikeTs, timeLockStruct.validEventTimes(i), timeLockStruct.window);
    timeLockStruct.spdfByEvent(i,:) = fixedPsth(alignedSpikeTimesThisEvent, kernelSigma, 0, timeLockStruct.t); % no error
end

timeLockStruct.spikeTimesByLoc = cell(nLoc, 1);
timeLockStruct.spdfByLoc = nan(nLoc, numel(timeLockStruct.t));
timeLockStruct.spdfErrByLoc = nan(nLoc, numel(timeLockStruct.t));
for i = 1:nLoc
    timeLockStruct.spikeTimesByLoc{i} = createnonemptydatamatpt(spikeTs, timeLockStruct.validEventTimesByLoc{i}, timeLockStruct.window);
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
