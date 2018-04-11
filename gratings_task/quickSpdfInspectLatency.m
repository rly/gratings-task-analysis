function quickSpdfInspectLatency(ES, blockName, ...
        D, spikeStructInd, spikeStruct, nLoc, isZeroDistractors, plotFileName)

unitName = spikeStruct.name;
nTrials = numel(ES.UE.cueOnset);
cols = lines(6);

%%
f = figure_tr_inch(13, 7.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('Gratings Attention Task: %s (%d trials)', unitName, nTrials);
if isZeroDistractors
    modTitle = [modTitle ' (0 Distractors)'];
end
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
waveformW = 0.1;
spdfW = 0.354; % == 0.098*3 + 0.03*2
waveformH = 0.15;
spdfH = 0.39;

waveformLeft1 = 0.05;
spdfCueOnsetLeft = waveformLeft1 + waveformW + 0.08;
spdfArrayOnsetRelLeft = spdfCueOnsetLeft + spdfW + 0.045;
spdfArrayOnsetHoldLeft = spdfCueOnsetLeft;
spdfTargetDimLeft = spdfArrayOnsetRelLeft;

btm = 0.07;
infoTextTop = 0.02;
infoText2Top = btm + 0.64;
rasterBtm = btm + spdfH + 0.08;
spdfTopRowBtm = rasterBtm;
waveformBtm = 0.76;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D, spikeStructInd, spikeStruct.isMUA);

%% info
writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop, 1);

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft spdfTopRowBtm spdfW spdfH]); 

t = ES.cueOnset.t - ES.cueOnset.window(1);
xBounds = [-0.25 0.25];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.cueOnset.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.cueOnset.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.cueOnset.latencyByLoc(i))
        plot(ES.cueOnset.latencyByLoc(i), ES.cueOnset.spdfByLoc(i,ES.cueOnset.latencyInfoByLoc{i}.latencyTInd), ...
                '+', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if ~isnan(ES.cueOnset.latencyBootByLoc(i))
        plot(ES.cueOnset.latencyBootByLoc(i), ES.cueOnset.spdfByLoc(i,ES.cueOnset.latencyBootInfoByLoc{i}.latencyTInd), ...
                'o', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.cueOnset.spikeTimesByLoc{i}));
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([-0.35 -0.35], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
% fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
%         [0.7 0.7 1], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
%         [0.3 1 0.3], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Cue Onset (s)');
ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset release
axArrayOnsetRelSpdf = axes('Position', [spdfArrayOnsetRelLeft spdfTopRowBtm spdfW spdfH]); 

t = ES.arrayOnsetRel.t - ES.arrayOnsetRel.window(1);
xBounds = [-0.25 0.25];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetRel.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.arrayOnsetRel.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.arrayOnsetRel.latencyByLoc(i))
        plot(ES.arrayOnsetRel.latencyByLoc(i), ES.arrayOnsetRel.spdfByLoc(i,ES.arrayOnsetRel.latencyInfoByLoc{i}.latencyTInd), ...
                '+', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if ~isnan(ES.arrayOnsetRel.latencyBootByLoc(i))
        plot(ES.arrayOnsetRel.latencyBootByLoc(i), ES.arrayOnsetRel.spdfByLoc(i,ES.arrayOnsetRel.latencyBootInfoByLoc{i}.latencyTInd), ...
                'o', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetRel.spikeTimesByLoc{i}));
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([-0.35 -0.35], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
% fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
%         [0.7 0.7 1], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
%         [0.3 1 0.3], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Array Onset (s)');
% ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'Release Trials Only', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset hold
axArrayOnsetHoldSpdf = axes('Position', [spdfArrayOnsetHoldLeft btm spdfW spdfH]); 

t = ES.arrayOnsetHold.t - ES.arrayOnsetHold.window(1);
xBounds = [-0.25 0.25];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetHold.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.arrayOnsetHold.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.arrayOnsetHold.latencyByLoc(i))
        plot(ES.arrayOnsetHold.latencyByLoc(i), ES.arrayOnsetHold.spdfByLoc(i,ES.arrayOnsetHold.latencyInfoByLoc{i}.latencyTInd), ...
                '+', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if ~isnan(ES.arrayOnsetHold.latencyBootByLoc(i))
        plot(ES.arrayOnsetHold.latencyBootByLoc(i), ES.arrayOnsetHold.spdfByLoc(i,ES.arrayOnsetHold.latencyBootInfoByLoc{i}.latencyTInd), ...
                'o', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetHold.spikeTimesByLoc{i}));
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([-0.35 -0.35], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
% fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
%         [0.7 0.7 1], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
%         [0.3 1 0.3], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Array Onset (s)');
ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'Hold Trials Only', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

t = ES.targetDim.t - ES.targetDim.window(1);
xBounds = [-0.25 0.25];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.targetDim.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.targetDim.spdfByLoc(i,:),  ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.targetDim.latencyByLoc(i))
        plot(ES.targetDim.latencyByLoc(i), ES.targetDim.spdfByLoc(i,ES.targetDim.latencyInfoByLoc{i}.latencyTInd), ...
                '+', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
    if ~isnan(ES.targetDim.latencyBootByLoc(i))
        plot(ES.targetDim.latencyBootByLoc(i), ES.targetDim.spdfByLoc(i,ES.targetDim.latencyBootInfoByLoc{i}.latencyTInd), ...
                'o', 'Color', cols(i,:) * 0.9, 'MarkerSize', 10, 'LineWidth', 2);
    end
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([0.28 0.28], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min rt

% shade in baseline and analysis windows
% fillH = jbfill(ES.targetDimDelayWindowOffset, [0 0], [1000 1000], ...
%         [0.7 0.7 1], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% fillH = jbfill(ES.targetDimResponseWindowOffset, [0 0], [1000 1000], ...
%         [0.3 1 0.3], ones(3, 1), 0.1);
% uistack(fillH, 'bottom');
% hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Target Dim (s)');
% ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% adjust all ylim
allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetRelSpdf) ylim(axArrayOnsetHoldSpdf) ylim(axTargetDimSpdf)];
yBounds = [min(allYBounds) max(1, max(allYBounds))];
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetRelSpdf, yBounds);
ylim(axArrayOnsetHoldSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);

%% add response metrics text

textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
text(axBig, -0.03, infoText2Top, {...
        sprintf('Blocks: %s', blockName), ...
        sprintf('Num Randomizations: %d', ES.numRandomizations), ...
        '', ...
        sprintf('{\\color[rgb]{%f,%f,%f}InRF (mean cue): %d} ({\\color[rgb]{%f,%f,%f}peak: %d}, {\\color[rgb]{%f,%f,%f}extreme: %d})', ...
            cols(ES.inRFLoc,:), ES.inRFLoc, ...
            cols(ES.inRFLocByMax,:), ES.inRFLocByMax, ...
            cols(ES.inRFLocByExtreme,:), ES.inRFLocByExtreme), ...
        sprintf('Baseline (SPDF): %0.1f Hz, Max: %0.1f Hz', ES.averageFiringRatesBySpdf.preCueBaseline.all, ES.maxFiringRateBySpdfInclMotor), ...
        sprintf('Cue {\\color[rgb]{%f,%f,%f}InRF}: %0.1f Hz, {\\color[rgb]{%f,%f,%f}ExRF}: %0.1f Hz', ...
            cols(ES.inRFLoc,:), ES.averageFiringRatesBySpdf.cueResponse.byLoc(ES.inRFLoc), ...
            cols(ES.exRFLoc,:), ES.averageFiringRatesBySpdf.cueResponse.byLoc(ES.exRFLoc)), ...
        '', ...
        sprintf('+ : Time to Half Peak'), ...
        sprintf('o : First Significance Over Baseline'), ...
        '', ...
        sprintf('Cue Latency:'), ...
        sprintf('+ : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.cueOnset.latencyByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.cueOnset.latencyByLoc(ES.exRFLoc) * 1000)), ...
        sprintf('o : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.exRFLoc) * 1000)), ...        '', ...
        '', ...
        sprintf('Array Rel Latency:'), ...
        sprintf('+ : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.arrayOnsetRel.latencyByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.arrayOnsetRel.latencyByLoc(ES.exRFLoc) * 1000)) ... 
        sprintf('o : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.arrayOnsetRel.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.arrayOnsetRel.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
        '', ...
        sprintf('Array Hold Latency:'), ...
        sprintf('+ : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.arrayOnsetHold.latencyByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.arrayOnsetHold.latencyByLoc(ES.exRFLoc) * 1000)) ... 
        sprintf('o : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.arrayOnsetHold.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.arrayOnsetHold.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
        '', ...
        sprintf('Dimming Latency:'), ...
        sprintf('+ : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.targetDim.latencyByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.targetDim.latencyByLoc(ES.exRFLoc) * 1000)) ... 
        sprintf('o : {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
        '', ...
        }, ...
        textParams{:});

%         sprintf('Time to Half Peak Method (+): '), ...
%         sprintf('Cue Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.cueOnset.latencyByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.cueOnset.latencyByLoc(ES.exRFLoc) * 1000)), ...
%         sprintf('Array Rel Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.arrayOnsetRel.latencyByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.arrayOnsetRel.latencyByLoc(ES.exRFLoc) * 1000)) ... 
%         sprintf('Array Hold Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.arrayOnsetHold.latencyByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.arrayOnsetHold.latencyByLoc(ES.exRFLoc) * 1000)) ... 
%         sprintf('Dimming Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.targetDim.latencyByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.targetDim.latencyByLoc(ES.exRFLoc) * 1000)) ... 
%         '', ...
%         sprintf('First Significance Over Baseline Method (o):'), ...
%         sprintf('Cue Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.exRFLoc) * 1000)), ...
%         sprintf('Array Rel Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.arrayOnsetRel.latencyBootByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.arrayOnsetRel.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
%         sprintf('Array Hold Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.arrayOnsetHold.latencyBootByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.arrayOnsetHold.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
%         sprintf('Dimming Latency: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
%             cols(ES.inRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.inRFLoc) * 1000), ...
%             cols(ES.exRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
