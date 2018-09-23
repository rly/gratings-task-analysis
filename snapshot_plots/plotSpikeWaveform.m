function plotSpikeWaveform(allUnitStructs, unitInd)

%% parameters
yBounds = [-0.07 0.07];
unitLineWidth = 4;
otherUnitLineWidth = 1;
unitSDShadingOpacity = 0.4;
otherUnitSDShadingOpacity = 0.1;
numSDsShading = 1;
thresholdCol = [0.5 0.2 0.5]; % purple
nsClassThresh = 0.326; % ms
nsClassThreshCol = [0.5 0.2 0.2]; % dark red
xTick = [-0.4 0 0.4 0.8];

unitStruct = allUnitStructs{unitInd};
nWfTime = numel(unitStruct.meanWf);
spikeFs = unitStruct.Fs;
waveformT = (0:nWfTime-1)/(spikeFs/1000) - unitStruct.thresholdTime * 1000; % start at 0

%% SUA
if ~unitStruct.isMUA
    %% plot threshold and axes
    hold on;
    plot([waveformT(1) waveformT(end)], [0 0], 'Color', 0.5*ones(3, 1));
    plot([waveformT(1) waveformT(end)], [unitStruct.threshold unitStruct.threshold], '--', ...
            'Color', thresholdCol);
    plot([0 0], yBounds, '--', ...
            'Color', thresholdCol);
    plot([nsClassThresh nsClassThresh], yBounds, ':', ...
            'Color', nsClassThreshCol);

    %% plot other waveforms
    allUnitsThisChannel = findAllUnitsSameCh(allUnitStructs, unitInd);
    allCols = lines(numel(allUnitsThisChannel));
    for k = 1:numel(allUnitsThisChannel)
        otherUnitInd = allUnitsThisChannel(k);
        if otherUnitInd == unitInd % the current unit
            currentUnitCol = allCols(k,:);
        else
            otherUnitSpikeStruct = allUnitStructs{otherUnitInd};
            otherUnitSd = std(otherUnitSpikeStruct.wf);
            otherUnitShadingUB = otherUnitSpikeStruct.meanWf + numSDsShading * otherUnitSd;
            otherUnitShadingLB = otherUnitSpikeStruct.meanWf - numSDsShading * otherUnitSd;
            jbfill(waveformT, otherUnitShadingUB, otherUnitShadingLB, allCols(k,:), ...
                    ones(3, 1), otherUnitSDShadingOpacity);
            hold on;

            plot(waveformT, otherUnitSpikeStruct.meanWf, ...
                    'LineWidth', otherUnitLineWidth, 'Color', allCols(k,:));
        end
    end

    %% shading of current waveform
    sd = std(unitStruct.wf);
    shadingUB = unitStruct.meanWf + numSDsShading * sd;
    shadingLB = unitStruct.meanWf - numSDsShading * sd;
    jbfill(waveformT, shadingUB, shadingLB, currentUnitCol, ones(3, 1), unitSDShadingOpacity);
    hold on;

    %% plot current waveform on top of everything else
    plot(waveformT, unitStruct.meanWf, 'LineWidth', unitLineWidth, 'Color', currentUnitCol);
    
else % MUA
    threshold = nanmean(unitStruct.thresholdParams.thresholds);

    %% plot threshold and axes
    hold on;
    plot([waveformT(1) waveformT(end)], [0 0], 'Color', 0.5*ones(3, 1));
    plot([waveformT(1) waveformT(end)], [threshold threshold], '--', ...
            'Color', thresholdCol);
    plot([0 0], yBounds, '--', ...
            'Color', thresholdCol); % TODO this is actually the alignment time
    
    %% plot this waveform on top of everything else
    plot(waveformT, unitStruct.meanWf, 'LineWidth', unitLineWidth, 'Color', lines(1));
end

%% formatting and labels
xlabel('Time (ms)');
xlim([waveformT(1) waveformT(end)]);
set(gca, 'XTick', xTick);
ylim(yBounds);
title('Waveform');
