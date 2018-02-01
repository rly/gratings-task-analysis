function [t,averageResponses] = computeResponsesInWindow(channelData, eventTimes, ...
        periEventWindowOffset, Fs)

nChannels = size(channelData, 1);

% convert flash time to lfp variable index
% slightly more accurate than createdatamatc b/c of rounding after offset
startIndices = round((eventTimes + periEventWindowOffset(1)) * Fs); 
endIndices = startIndices + round(diff(periEventWindowOffset) * Fs) - 1;
t = periEventWindowOffset(1):1/Fs:periEventWindowOffset(2)-1/Fs;
nTime = numel(t);
assert(all(nTime == (endIndices - startIndices + 1)));

nEvents = numel(eventTimes);
averageResponses = nan(nChannels, nTime);

for j = 1:nChannels
    responses = nan(nTime, nEvents);
    for i = 1:nEvents
        responses(:,i) = channelData(j,startIndices(i):endIndices(i));
    end;
    averageResponses(j,:) = mean(responses, 2);
end