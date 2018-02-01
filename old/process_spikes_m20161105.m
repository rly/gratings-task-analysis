function process_spikes_m20161105(sessionName, spikeVar, spikeVarName, nLoc, ...
        arrayOnsetRel, arrayOnsetHold, ...
        arrayOnsetRelByLoc, arrayOnsetHoldByLoc, targetDim, ...
        targetDimByLoc, ...
        targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc, ...
        cueOnset, cueOnsetByLoc)

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
cueOnsetSpdfWindowOffset = [-0.4 0.5];
cueOnsetT = computeTForSpdf(periCueOnsetWindow(1), cueOnsetSpdfWindowOffset, kernelSigma);
cueOnsetSpdf = edpsth_notranspose(cueOnsetSpikeTimes, kernelSigma, 'n', [], 0, cueOnsetT);

figure;
hold on;
plot(cueOnsetT - periCueOnsetWindow(1), cueOnsetSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('%s - %s - Response to Cue Onset (N=%d)', ...
        sessionName, spikeVarName, numel(cueOnset)));

%%
figure;
quickSpdf(cueOnsetSpikeTimesByLoc, nLoc, cueOnsetT, periCueOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(cueOnsetSpdfWindowOffset);
xlabel('Time from Cue Onset (ms)');
title(sprintf('%s - %s - Response to Cue Onset (N=%d)', ...
        sessionName, spikeVarName, numel(cueOnset)));

plotFileName = sprintf('%s-%s-cueByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);

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
arrayOnsetSpdfWindowOffset = [-0.3 1];
arrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), arrayOnsetSpdfWindowOffset, kernelSigma);
arrayOnsetRelSpdf = edpsth_notranspose(arrayOnsetRelSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);
arrayOnsetHoldSpdf = edpsth_notranspose(arrayOnsetHoldSpikeTimes, kernelSigma, 'n', [], 0, arrayOnsetT);

figure;
hold on;
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetHoldSpdf);
plot(arrayOnsetT - periArrayOnsetWindow(1), arrayOnsetRelSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim([-0.3 0.4]);
xlabel('Time from Array Onset (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('%s - Response to Array Onset (N=%d)', sessionName, ...
        numel(arrayOnsetHold) + numel(arrayOnsetRel)));
legend({sprintf('Target: Hold (N=%d)', numel(arrayOnsetHold)), ...
        sprintf('Target: Release (N=%d)', numel(arrayOnsetRel))}, ...
        'Location', 'NorthWest');

%%
figure;
quickSpdf(arrayOnsetHoldSpikeTimesByLoc, nLoc, arrayOnsetT, periArrayOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
title(sprintf('%s - %s - Response to Array Onset - Hold Target (N=%d)', ...
        sessionName, spikeVarName, numel(arrayOnsetHold)));

plotFileName = sprintf('%s-%s-arrayOnsetHoldByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);
    
%%
figure;
quickSpdf(arrayOnsetRelSpikeTimesByLoc, nLoc, arrayOnsetT, periArrayOnsetWindow, kernelSigma);
ylim(origYLim);
xlim(arrayOnsetSpdfWindowOffset);
xlabel('Time from Array Onset (ms)');
title(sprintf('%s - %s - Response to Array Onset - Release Target (N=%d)', ...
        sessionName, spikeVarName, numel(arrayOnsetRel)));

plotFileName = sprintf('%s-%s-arrayOnsetRelByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);

%% align spikes to TARGET DIM
periTargetDimWindow = [1 1];
targetDimSpikeTimes = createnonemptydatamatpt(spikeVar, targetDim, periTargetDimWindow);

targetDimSpikeTimesByLoc = cell(nLoc,1);
targetDimShortHoldDurSpikeTimesByLoc = cell(nLoc,1);
targetDimLongHoldDurSpikeTimesByLoc = cell(nLoc,1);
for i = 1:nLoc
    targetDimSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimByLoc{i}, periTargetDimWindow);
    targetDimShortHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimShortHoldDurByLoc{i}, periTargetDimWindow);
    targetDimLongHoldDurSpikeTimesByLoc{i} = createnonemptydatamatpt(spikeVar, targetDimLongHoldDurByLoc{i}, periTargetDimWindow);
end

%% spdf
targetDimSpdfWindowOffset = [-1 0.7];
targetDimT = computeTForSpdf(periTargetDimWindow(1), targetDimSpdfWindowOffset, kernelSigma);
targetDimSpdf = edpsth_notranspose(targetDimSpikeTimes, kernelSigma, 'n', [], 0, targetDimT);

figure;
hold on;
plot(targetDimT - periTargetDimWindow(1), targetDimSpdf);
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
ylabel('Estimated Firing Rate (Hz)');
title(sprintf('%s - %s - Response to Target Dimming (N=%d)', ...
        sessionName, spikeVarName, numel(targetDim)));

%%
figure;
quickSpdf(targetDimSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
title(sprintf('%s - %s - Response to Target Dimming (N=%d)', ...
        sessionName, spikeVarName, numel(targetDim)));

plotFileName = sprintf('%s-%s-targetDimByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);
    
%%
figure;
quickSpdf(targetDimShortHoldDurSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
numTrialsShortHoldAll = sum(cellfun(@numel, targetDimShortHoldDurSpikeTimesByLoc));
title(sprintf('%s - %s - Response to Target Dimming - SHORT Hold (N=%d)', ...
        sessionName, spikeVarName, numTrialsShortHoldAll));

plotFileName = sprintf('%s-%s-targetDimShortHoldByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);


%%
figure;
quickSpdf(targetDimLongHoldDurSpikeTimesByLoc, nLoc, targetDimT, periTargetDimWindow, kernelSigma);
ylim(origYLim);
xlim(targetDimSpdfWindowOffset);
xlabel('Time from Target Dimming (ms)');
numTrialsLongHoldAll = sum(cellfun(@numel, targetDimLongHoldDurSpikeTimesByLoc));
title(sprintf('%s - %s - Response to Target Dimming - LONG Hold (N=%d)', ...
        sessionName, spikeVarName, numTrialsLongHoldAll));

plotFileName = sprintf('%s-%s-targetDimLongHoldByLoc.png', sessionName, spikeVarName);
export_fig(plotFileName);

