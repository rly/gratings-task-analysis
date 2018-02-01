function [diffSpikeTimesInWindow, silentPeriodDurationsPrecedingSpikeInWindow, numSpikesTotalInWindow, spikesInBurst, averageISIWithinBurst] = ...
        burstAnalysis2(spikeTimes, startWindowTimes, endWindowTimes, minSilentPeriod, minISIWithinBurst)


assert(all(size(startWindowTimes) == size(endWindowTimes)));
assert(all(startWindowTimes < endWindowTimes));

numWindows = numel(startWindowTimes);
diffSpikeTimesInWindow = cell(numWindows, 1);
indicesWithPrecedingSilentPeriodInWindow = cell(numWindows, 1);
silentPeriodDurationsPrecedingSpikeInWindow = [];
numSpikesTotalInWindow = 0;
for i = 1:numWindows
    spikeTimesInTrialWindow = spikeTimes(spikeTimes >= startWindowTimes(i) & spikeTimes <= endWindowTimes(i));
    diffSpikeTimesInWindow{i} = diff(spikeTimesInTrialWindow);
    numSpikesTotalInWindow = numSpikesTotalInWindow + numel(spikeTimesInTrialWindow);
    indicesWithPrecedingSilentPeriodInWindow{i} = find(diffSpikeTimesInWindow{i} > minSilentPeriod);
    silentPeriodDurationsPrecedingSpikeInWindow = [silentPeriodDurationsPrecedingSpikeInWindow; ...
            diffSpikeTimesInWindow{i}(indicesWithPrecedingSilentPeriodInWindow{i})];
end

countValidISIs = nan(numWindows, 1); % pre-allocate
averageISIWithinBurst = nan(numWindows, 1); % pre-allocate
countValidISIsIndex = 1;
for i = 1:numWindows
    for j = 1:numel(indicesWithPrecedingSilentPeriodInWindow{i})
        countValidISIs(countValidISIsIndex) = 0;
        startIndex = indicesWithPrecedingSilentPeriodInWindow{i}(j);
        % count the number of diffSpikeTimes <= minISIWithinBurst starting from
        % the one after the spike after a silent period
        while startIndex + countValidISIs(countValidISIsIndex) + 1 < numel(diffSpikeTimesInWindow{i}) && ...
                diffSpikeTimesInWindow{i}(startIndex + countValidISIs(countValidISIsIndex) + 1) <= minISIWithinBurst
            countValidISIs(countValidISIsIndex) = countValidISIs(countValidISIsIndex) + 1;
        end
        averageISIWithinBurst(countValidISIsIndex) = ...
                mean(diffSpikeTimesInWindow{i}(startIndex:startIndex + countValidISIs(countValidISIsIndex)));
        countValidISIsIndex = countValidISIsIndex + 1;
    end
end

countValidISIs(isnan(countValidISIs)) = [];
averageISIWithinBurst(isnan(averageISIWithinBurst)) = [];

% convert ISI count to spike count and remove bursts of size 1
spikesInBurst = countValidISIs + 1;
spikesInBurst(spikesInBurst == 1) = [];

