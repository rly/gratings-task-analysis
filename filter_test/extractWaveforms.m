function [extractedWaveforms,startWaveformInds] = extractWaveforms(data, isExtreme, preExtremeSamples, ...
        postExtremeSamples, deadSamplesAfterThreshCrossing, isAlignTrough)
% each row is a waveform, columns are samples
assert(all(size(data) == size(isExtreme)));

alignToExtremeMaxShiftPreThreshSamples = 0;
alignToExtremeMaxShiftPostThreshSamples = 10;

% threshold crossing is when data goes from non-extreme to extreme (0 -> 1)
diffIsExtreme = diff(isExtreme);
threshCrossing = find(diffIsExtreme == 1) + 1; % +1 to capture the extreme not the pre-extreme
fprintf('\tFound %d threshold crossings, ', numel(threshCrossing));

isGoodThreshCrossing = true(size(isExtreme));
extractedWaveforms = nan(numel(isExtreme), preExtremeSamples + postExtremeSamples + 1);
startWaveformInds = nan(numel(isExtreme), 1);
for i = 1:numel(threshCrossing)
    if isGoodThreshCrossing(i)
        % disable any threshold crossings within deadSamplesAfterThreshCrossing after this one
        diffThreshCrossingWithThis = threshCrossing - threshCrossing(i);
        isGoodThreshCrossing(diffThreshCrossingWithThis > 0 & diffThreshCrossingWithThis <= deadSamplesAfterThreshCrossing) = false;
        
        % find trough within sample, up to alignToExtremeMaxShift away
        lb = threshCrossing(i) - alignToExtremeMaxShiftPreThreshSamples;
        ub = threshCrossing(i) + alignToExtremeMaxShiftPostThreshSamples;
        if lb < 1 || ub > numel(data)
            continue; % skip if the waveform extends past data boundaries
        end
        extractedSampleWaveform = data(lb:ub);
        
        if isAlignTrough
            [~,extremeInd] = min(extractedSampleWaveform);
        else
            [~,extremeInd] = max(extractedSampleWaveform);
        end
        extremeGlobalInd = threshCrossing(i) - alignToExtremeMaxShiftPreThreshSamples + extremeInd - 1;
        
        lb = extremeGlobalInd - preExtremeSamples;
        ub = extremeGlobalInd + postExtremeSamples;
        if lb < 1 || ub > numel(data)
            continue; % skip if the waveform extends past data boundaries
        end
        extractedWaveforms(i,:) = data(lb:ub);
        startWaveformInds(i) = lb;
    end
end

extractedWaveforms = trimNanRows(extractedWaveforms);
startWaveformInds = trimNanRows(startWaveformInds);
fprintf('extracted %d waveforms.\n', size(extractedWaveforms, 1));
