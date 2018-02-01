function [extractedWaveformsLowerThresh,startWaveformLowerThreshInds,extractedWaveformsUpperThresh,startWaveformUpperThreshInds] = findSpikeWaveformsUsingThreshold(...
        data, numSDsThresh, preExtremeSamples, postExtremeSamples, deadSamplesAfterThreshCrossing, isUseMAD, channelID)
% align to trough for lower, peak for upper

data = makeRowVector(data);

%% compute noise level
% remove waveforms 10k samples pre (250 ms @ 40 kHz) and 30k samples post
% (750 ms @ 40 kHz) threshold crossings
initNumSDsThresh = 3; % both + and -
initPreExtremeSamples = 10;
initPostExtremeSamples = 30;

isExtreme = findDataCrossingThreshold(...
        data, initNumSDsThresh, initPreExtremeSamples, initPostExtremeSamples, isUseMAD);

fprintf('Removed %d/%d (%d%%) data points in noise calculation.\n', sum(isExtreme), numel(data), round(sum(isExtreme)/numel(data)*100));

%% find threshold crossings using above noise level
meanAdjData = nanmean(data(~isExtreme));
if isUseMAD
    sdAdjData = mad(data(~isExtreme), 1);
else
    sdAdjData = nanstd(data(~isExtreme));
end

fprintf('Adjusted mean: %0.3f, Adjusted SD: %0.3f\n', meanAdjData, sdAdjData);

adjLowerThresh = meanAdjData - numSDsThresh * sdAdjData;
adjUpperThresh = meanAdjData + numSDsThresh * sdAdjData;

isLowerExtreme = data < adjLowerThresh;
isUpperExtreme = data > adjUpperThresh;

%% extract waveforms
% lower and upper waveforms are computed independently
% dead time is not considered across them
% the same waveform might cross both the upper and lower threshold
fprintf('Lower threshold: %0.3f\n', adjLowerThresh);
[extractedWaveformsLowerThresh,startWaveformLowerThreshInds] = extractWaveforms(data, isLowerExtreme, preExtremeSamples, postExtremeSamples, deadSamplesAfterThreshCrossing, 1);

fprintf('Upper threshold: %0.3f\n', adjUpperThresh);
[extractedWaveformsUpperThresh,startWaveformUpperThreshInds] = extractWaveforms(data, isUpperExtreme, preExtremeSamples, postExtremeSamples, deadSamplesAfterThreshCrossing, 0);

%% plot
figure_tr_inch(6, 10);
subaxis(2, 1, 1);
hold on;
plot(extractedWaveformsUpperThresh');
title(sprintf('Channel %d, Upper Threshold (N=%d)', channelID, size(extractedWaveformsUpperThresh, 1)));

subaxis(2, 1, 2);
hold on;
plot(extractedWaveformsLowerThresh');
title(sprintf('Channel %d, Lower Threshold: (N=%d)', channelID, size(extractedWaveformsLowerThresh, 1)));
