function quickSpdfAllEvents(processedDataDir, blockName, spikeStruct, nLoc, ...
        arrayOnsetRel, arrayOnsetHold, ...
        arrayOnsetRelByLoc, arrayOnsetHoldByLoc, targetDim, ...
        targetDimByLoc, ...
        targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
        cueOnset, cueOnsetByLoc)

spikeVar = spikeStruct.ts;
spikeVarName = spikeStruct.name;
filePrefix = sprintf('%s-%s', spikeVarName, blockName);

%% set ylim
origYLim = [0 40]; 
kernelSigma = 0.01;

%% align spikes to CUE
periCueOnsetWindow = [1 1];
cueOnsetSpikeTimes = createnonemptydatamatpt(spikeVar, cueOnset, periCueOnsetWindow);

cueOnsetSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    cueOnsetSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{i}, periCueOnsetWindow);
end

%% spdf
cueOnsetSpdfWindowOffset = [-0.7 0.6];
cueOnsetT = computeTForSpdf(periCueOnsetWindow(1), cueOnsetSpdfWindowOffset, kernelSigma);
cueOnsetSpdf = edpsth_notranspose(cueOnsetSpikeTimes, kernelSigma, 'n', [], 0, cueOnsetT);

figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
hold on;
plot(cueOnsetT - periCueOnsetWindow(1), cueOnsetSpdf, 'LineWidth', 2);
plot([0 0], [0 100], '-', 'Color', 0.3 * ones(3, 1));
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title({spikeVarName, sprintf('Response to Cue Onset (N=%d)', ...
        numel(cueOnset))}, 'Interpreter', 'none');
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 16);
% don't save

%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(cueOnsetSpikeTimesByLoc, nLoc, cueOnsetT, periCueOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
title({spikeVarName, sprintf('Response to Cue Onset (N=%d)', ...
        numel(cueOnset))}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-cueByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop'); 

%% align spikes to ARRAY
periArrayOnsetWindow = [1 1]; % seconds before, after

arrayOnsetRelSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetRel, periArrayOnsetWindow);
arrayOnsetHoldSpikeTimes = createnonemptydatamatpt(spikeVar, arrayOnsetHold, periArrayOnsetWindow);

arrayOnsetRelSpikeTimesByLoc = cell(nLoc,1);
arrayOnsetHoldSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    arrayOnsetRelSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetRelByLoc{i}, periArrayOnsetWindow);
    arrayOnsetHoldSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{i}, periArrayOnsetWindow);
end

%% spdf
arrayOnsetSpdfWindowOffset = [-0.6 0.7];
arrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), arrayOnsetSpdfWindowOffset, kernelSigma);
arrayOnsetRelSpdf = edpsth_notranspose(arrayOnsetRelSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);
arrayOnsetHoldSpdf = edpsth_notranspose(arrayOnsetHoldSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);

figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
hold on;
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetHoldSpdf, 'LineWidth', 2);
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetRelSpdf, 'LineWidth', 2);
plot([0 0], [0 100], '-', 'Color', 0.3 * ones(3, 1));
ylim(origYLim);
xlim([-0.3 0.4]);
xlabel('Time from Array Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title({spikeVarName, sprintf('Response to Array Onset (N=%d)', ...
        numel(arrayOnsetHold) + numel(arrayOnsetRel))}, 'Interpreter', 'none');
legend({sprintf('Target: Hold (N=%d)', numel(arrayOnsetHold)), ...
        sprintf('Target: Release (N=%d)', numel(arrayOnsetRel))}, ...
        'Location', 'NorthWest');
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 16);
% don't save

%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(arrayOnsetHoldSpikeTimesByLoc, nLoc, arrayOnsetT, periArrayOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
title({spikeVarName, sprintf('Response to Array Onset - Hold Target (N=%d)', ...
        numel(arrayOnsetHold))}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-arrayOnsetHoldByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop');
    
%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(arrayOnsetRelSpikeTimesByLoc, nLoc, arrayOnsetT, periArrayOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
title({spikeVarName, sprintf('Response to Array Onset - Release Target (N=%d)', ...
        numel(arrayOnsetRel))}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-arrayOnsetRelByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop'); 

%% align spikes to TARGET DIM
periTargetDimWindow = [1 1];
targetDimSpikeTimes = createnonemptydatamatpt(spikeVar, targetDim, periTargetDimWindow);

targetDimSpikeTimesByLoc = cell(nLoc,1);
targetDimShortHoldDurSpikeTimesByLoc = cell(nLoc,1);
targetDimLongHoldDurSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    targetDimSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimByLoc{i}, periTargetDimWindow);
    targetDimShortHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimShortHoldByLoc{i}, periTargetDimWindow);
    targetDimLongHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimLongHoldByLoc{i}, periTargetDimWindow);
end

%% spdf
targetDimSpdfWindowOffset = [-0.7 0.7];
targetDimT = computeTForSpdf(periTargetDimWindow(1), targetDimSpdfWindowOffset, kernelSigma);
targetDimSpdf = edpsth_notranspose(targetDimSpikeTimes, kernelSigma, 'n', [], 0, targetDimT);

figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
hold on;
plot(targetDimT - periTargetDimWindow(1), targetDimSpdf, 'LineWidth', 2);
plot([0 0], [0 100], '-', 'Color', 0.3 * ones(3, 1));
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title({spikeVarName, sprintf('Response to Target Dimming (N=%d)', ...
        numel(targetDim))}, 'Interpreter', 'none');
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 16);
% don't save

%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(targetDimSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
title({spikeVarName, sprintf('Response to Target Dimming (N=%d)', ...
        numel(targetDim))}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-targetDimByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop'); 
    
%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(targetDimShortHoldDurSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
numTrialsShortHoldAll = sum(cellfun(@numel, targetDimShortHoldDurSpikeTimesByLoc));
title({spikeVarName, sprintf('Response to Target Dimming - SHORT Hold (N=%d)', ...
        numTrialsShortHoldAll)}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-targetDimShortHoldByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop'); 

%%
figure_tr_inch(8, 6);
subaxis(1, 1, 1, 'MB', 0.13, 'MT', 0.12);
quickSpdf(targetDimLongHoldDurSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
numTrialsLongHoldAll = sum(cellfun(@numel, targetDimLongHoldDurSpikeTimesByLoc));
title({spikeVarName, sprintf('Response to Target Dimming - LONG Hold (N=%d)', ...
        numTrialsLongHoldAll)}, 'Interpreter', 'none');

plotFileName = sprintf('%s/%s-targetDimLongHoldByLoc.png', processedDataDir, filePrefix);
export_fig(plotFileName, '-nocrop'); 

%%
close all;
saveFileName = sprintf('%s/%s-processedSpikes.mat', processedDataDir, filePrefix);
save(saveFileName);