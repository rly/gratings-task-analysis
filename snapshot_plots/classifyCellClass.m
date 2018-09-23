function spikeStruct = classifyCellClass(spikeStruct)

nsClassThresh = 0.326/1000; % <= this is Narrow-Spiking
uncClassThresh = 0.326/1000; % > above and <= this is Unclassified

% conditions:
% there must be a peak above 0 V after smoothing
if ~isempty(spikeStruct.peakSmoothedAmps) && ...
        all(spikeStruct.peakSmoothedAmps > 0) && ...
        (~isnan(spikeStruct.fullWidthHalfMax) && spikeStruct.fullWidthHalfMax * 1000 > 0.25) && ... 
        (strcmp(spikeStruct.inflectionPattern, 'tp') || ...
            (strcmp(spikeStruct.inflectionPattern, 'ptp') && ...
            spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
            spikeStruct.peakSmoothedAmps(1) < 1/2 * max(spikeStruct.peakSmoothedAmps)))
    if spikeStruct.troughToPeakTimeFine <= nsClassThresh
        spikeStruct.physClass = 'Narrow-Spiking';
    elseif spikeStruct.troughToPeakTimeFine <= uncClassThresh
        spikeStruct.physClass = 'Unclassified Trough-Peak';
    else
        spikeStruct.physClass = 'Broad-Spiking';
    end
elseif regexp(spikeStruct.inflectionPattern, '^p+t$')
    spikeStruct.physClass = 'Peak-Trough Axonal';
elseif regexp(spikeStruct.inflectionPattern, 'pt+p')
    spikeStruct.physClass = 'Multiphasic Axonal';
else
    spikeStruct.physClass = 'Unknown';
end