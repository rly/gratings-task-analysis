function plotSpikeWaveform(D, unitInd, isMUA)

%% parameters
yBounds = [-0.07 0.07];
unitLineWidth = 4;
otherUnitLineWidth = 1;
unitSDShadingOpacity = 0.4;
otherUnitSDShadingOpacity = 0.1;
numSDsShading = 1;

%% setup for single unit
if ~isMUA
    spikeStruct = D.allSpikeStructs{unitInd};
    nWfTime = numel(spikeStruct.meanWf);
    spikeFs = D.timestampFrequency;
    waveformT = (1:nWfTime)/(spikeFs/1000);

    %% plot threshold and axes
    hold on;
    thresholdCol = [0.5 0.2 0.5]; % purple
    plot([waveformT(1) waveformT(end)], [0 0], 'Color', 0.5*ones(3, 1));
    plot([waveformT(1) waveformT(end)], [spikeStruct.threshold spikeStruct.threshold], '--', ...
            'Color', thresholdCol);
    plot([spikeStruct.thresholdTime spikeStruct.thresholdTime]*1000, yBounds, '--', ...
            'Color', thresholdCol);

    %% plot other waveforms
    allUnitsThisChannel = findAllUnitsSameCh(D.allSpikeStructs, unitInd);
    allCols = lines(numel(allUnitsThisChannel));
    for k = 1:numel(allUnitsThisChannel)
        otherUnitInd = allUnitsThisChannel(k);
        if otherUnitInd == unitInd % the current unit
            currentUnitCol = allCols(k,:);
        else
            otherUnitSpikeStruct = D.allSpikeStructs{otherUnitInd};
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

    %% shading
    sd = std(spikeStruct.wf);
    shadingUB = spikeStruct.meanWf + numSDsShading * sd;
    shadingLB = spikeStruct.meanWf - numSDsShading * sd;
    jbfill(waveformT, shadingUB, shadingLB, currentUnitCol, ones(3, 1), unitSDShadingOpacity);
    hold on;

    %% plot this waveform on top of everything else
    plot(waveformT, spikeStruct.meanWf, 'LineWidth', unitLineWidth, 'Color', currentUnitCol);
    
else % MUA
    muaStruct = D.allMUAStructs{unitInd};
    nWfTime = numel(muaStruct.meanWf);
    spikeFs = D.timestampFrequency;
    waveformT = (1:nWfTime)/(spikeFs/1000);

    %% plot threshold and axes
    hold on;
    plot([waveformT(1) waveformT(end)], [0 0], 'Color', 0.5*ones(3, 1));
    
    %% plot this waveform on top of everything else
    plot(waveformT, muaStruct.meanWf, 'LineWidth', unitLineWidth);
end

%% formatting and labels
xlabel('Time (ms)');
xlim([0 waveformT(end)]);
set(gca, 'XTick', [0 waveformT(end)/2 waveformT(end)]);
ylim(yBounds);
title('Waveform');
