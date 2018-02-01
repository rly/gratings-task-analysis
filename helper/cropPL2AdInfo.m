function newAdInfo = cropPL2AdInfo(adInfo, startTime, endTime)
% startTime, endTime in seconds

assert(startTime <= endTime);

% do not run this on data after NaN padding
assert(sum(adInfo.FragCounts) == numel(adInfo.Values));

% indexes into adInfo.FragTs and adInfo.FragCounts
blockStartInd = find(adInfo.FragTs <= startTime, 1, 'last');
blockEndInd = find(adInfo.FragTs <= endTime, 1, 'last');

if isempty(blockStartInd)
    blockStartInd = 1;
    startTime = adInfo.FragTs(1);
end

% indexes into adInfo.Values
endBlockValuesIndex = cumsum(adInfo.FragCounts);
startBlockValuesIndex = [1; endBlockValuesIndex(1:end-1)+1];

newAdInfo.FragTs = adInfo.FragTs(blockStartInd:blockEndInd);
newAdInfo.FragCounts = adInfo.FragCounts(blockStartInd:blockEndInd);
newAdInfo.ADFreq = adInfo.ADFreq;

% find how many samples after the block start is startTime
deltaIndicesStartTimeFromBlockStart = ceil((startTime - adInfo.FragTs(blockStartInd)) * adInfo.ADFreq);

newAdInfo.FragTs(1) = newAdInfo.FragTs(1) + deltaIndicesStartTimeFromBlockStart / adInfo.ADFreq;
newAdInfo.FragCounts(1) = newAdInfo.FragCounts(1) - deltaIndicesStartTimeFromBlockStart;
newStartValuesIndex = startBlockValuesIndex(blockStartInd) + deltaIndicesStartTimeFromBlockStart;

% find how many samples before the block end is endTime
lastTimestampInValues = adInfo.FragTs(end) + (adInfo.FragCounts(end) - 1) / adInfo.ADFreq;
if endTime > lastTimestampInValues
    newEndValuesIndex = numel(adInfo.Values);
else
    deltaIndicesEndTimeFromBlockEnd = floor((endTime - adInfo.FragTs(blockEndInd)) * adInfo.ADFreq);
    newAdInfo.FragCounts(end) = deltaIndicesEndTimeFromBlockEnd + 1;
    newEndValuesIndex = startBlockValuesIndex(blockEndInd) + deltaIndicesEndTimeFromBlockEnd;
end
newAdInfo.Values = adInfo.Values(newStartValuesIndex:newEndValuesIndex);
