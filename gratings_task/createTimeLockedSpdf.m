function timeLockStruct = createTimeLockedSpdf(spikeTs, eventTimes, eventTimesByLoc, timeLockStruct, kernelSigma, startTime, endTime)
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
timeLockStruct.spikeTimes = createnonemptydatamatpt(spikeTs, timeLockStruct.validEventTimes, timeLockStruct.window);

% compute the spike density function and bootstrapped error measure
[timeLockStruct.spdf,~,timeLockStruct.spdfErr] = fixedPsth(timeLockStruct.spikeTimes, kernelSigma, 2, timeLockStruct.t);

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
    [timeLockStruct.spdfByLoc(i,:),~,timeLockStruct.spdfErrByLoc(i,:)] = fixedPsth(timeLockStruct.spikeTimesByLoc{i}, kernelSigma, 2, timeLockStruct.t);
end
