function isExtreme = findDataCrossingThreshold(...
        data, numSDsThresh, preExtremeSamples, postExtremeSamples, isUseMAD)
data = makeRowVector(data);

%% compute noise level
% remove waveforms 10k samples pre (250 ms @ 40 kHz) and 30k samples post
% (750 ms @ 40 kHz) threshold crossings

meanData = nanmean(data);
if isUseMAD
    sdData = mad(data, 1);
else
    sdData = nanstd(data);
end

lowerThresh = meanData - numSDsThresh * sdData;
upperThresh = meanData + numSDsThresh * sdData;

extremeInds = find(data < lowerThresh | data > upperThresh);
isExtreme = false(size(data));
for i = 1:numel(extremeInds)
    lb = max(1, extremeInds(i) - preExtremeSamples); % keep inds within range
    ub = min(numel(data), extremeInds(i) + postExtremeSamples);
    isExtreme(lb:ub) = true;
end

