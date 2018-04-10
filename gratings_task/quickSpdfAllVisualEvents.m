function quickSpdfAllVisualEvents(ES, blockName, ...
        D, spikeStructInd, spikeStruct, nLoc, isZeroDistractors, plotFileName)

unitName = spikeStruct.name;
nTrials = numel(ES.UE.cueOnset);
cols = lines(4);
latencyCol = [0.6 0.05 0.6];

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
rasterW = 0.098;
cueOnsetSpdfW = 0.354; % == 0.098*3 + 0.03*2
arrayOnsetSpdfW = 0.354;
targetDimSpdfW = 0.354;
waveformH = 0.15;
rasterH = 0.39;
spdfH = 0.39;

waveformLeft1 = 0.05;
rasterCueOnsetLeft = waveformLeft1 + waveformW + 0.08;
rasterArrayOnsetLeft = rasterCueOnsetLeft + rasterW + 0.03;
rasterTargetDimLeft = rasterArrayOnsetLeft + rasterW + 0.03;
spdfCueOnsetLeft = rasterTargetDimLeft + rasterW + 0.045;
spdfArrayOnsetLeft = rasterCueOnsetLeft;
spdfTargetDimLeft = spdfCueOnsetLeft;

btm = 0.07;
infoTextTop = 0.02;
infoText2Top = btm + 0.64;
rasterBtm = btm + spdfH + 0.08;
spdfCueOnsetBtm = rasterBtm;
waveformBtm = 0.76;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D, spikeStructInd, spikeStruct.isMUA);

%% info
writeUnitInfo(spikeStruct, axBig, -0.03, infoTextTop, 1);

%% plot raster aligned to cue
axes('Position', [rasterCueOnsetLeft rasterBtm rasterW rasterH]); 
hold on;

xBounds = [-0.25 0.25];
window = ES.cueOnset.window;
data = ES.cueOnset.spikeTimes;
rasterY = 0;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
lineHeight = 5;
for k = 1:numel(data)
    rasterY = rasterY + 1;
    if ~isempty(data(k).times)
        adjTimes = data(k).times - window(1);
        for l = 1:numel(data(k).times)
            % plot mini lines at each spike time
            plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
        end
    end
end

% plot event line marker
yBounds = [0 rasterY+1];
plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1));

xlim(xBounds);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
set(gca, 'XTickLabel', []);
xlabel('Cue Onset');
ylabel('Trial Number');

%% plot raster aligned to array
axes('Position', [rasterArrayOnsetLeft rasterBtm rasterW rasterH]); 
hold on;

xBounds = [-0.25 0.25];
window = ES.arrayOnsetHold.window;
data = ES.arrayOnsetHold.spikeTimes; % TODO show release trials
rasterY = 0;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
lineHeight = 5;
for k = 1:numel(data)
    rasterY = rasterY + 1;
    if ~isempty(data(k).times)
        adjTimes = data(k).times - window(1);
        for l = 1:numel(data(k).times)
            % plot mini lines at each spike time
            plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
        end
    end
end

% plot event line marker
yBounds = [0 rasterY+1];
plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1));

xlim(xBounds);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
set(gca, 'XTickLabel', []);
xlabel('Array Onset');

%% plot raster aligned to target dim
axes('Position', [rasterTargetDimLeft rasterBtm rasterW rasterH]); 
hold on;

xBounds = [-0.25 0.25];
window = ES.targetDim.window;
data = ES.targetDim.spikeTimes; % TODO split by hold/release
rasterY = 0;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
lineHeight = 5;
for k = 1:numel(data)
    rasterY = rasterY + 1;
    if ~isempty(data(k).times)
        adjTimes = data(k).times - window(1);
        for l = 1:numel(data(k).times)
            % plot mini lines at each spike time
            plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
        end
    end
end

% plot event line marker
yBounds = [0 rasterY+1];
plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1));

xlim(xBounds);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
set(gca, 'XTickLabel', []);
xlabel('Target Dimming');

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft spdfCueOnsetBtm cueOnsetSpdfW spdfH]); 

t = ES.cueOnset.t - ES.cueOnset.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.cueOnset.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.cueOnset.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.cueOnset.latencyBootByLoc(i))
        plot(ES.cueOnset.latencyBootByLoc(i), ES.cueOnset.spdfByLoc(i,ES.cueOnset.latencyBootInfoByLoc{i}.latencyTInd), ...
                '.', 'Color', cols(i,:) * 0.9, 'MarkerSize', 18, 'LineWidth', 2);
    end
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.cueOnset.spikeTimesByLoc{i}));
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([-0.35 -0.35], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Cue Onset (s)');
ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm arrayOnsetSpdfW spdfH]); 

t = ES.arrayOnsetHold.t - ES.arrayOnsetHold.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetHold.spdfByLoc(i,:))) % TODO show release trials
        continue;
    end
    plot(t, ES.arrayOnsetHold.spdfByLoc(i,:),  ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.arrayOnsetHold.latencyBootByLoc(i))
        plot(ES.arrayOnsetHold.latencyBootByLoc(i), ES.arrayOnsetHold.spdfByLoc(i,ES.arrayOnsetHold.latencyBootInfoByLoc{i}.latencyTInd), ...
                '.', 'Color', cols(i,:) * 0.9, 'MarkerSize', 18, 'LineWidth', 2);
    end
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetHold.spikeTimesByLoc{i}));
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));

% shade in baseline and analysis windows
fillH = jbfill(ES.cueTargetDelayWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.arrayResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Array Onset (s)');
ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'Hold Trials Only', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm targetDimSpdfW spdfH]); 

t = ES.targetDim.t - ES.targetDim.window(1);
xBounds = [-0.6 0.6];
hold on;
for i = 1:nLoc
    if any(isnan(ES.targetDim.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.targetDim.spdfByLoc(i,:),  ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~isnan(ES.targetDim.latencyBootByLoc(i))
        plot(ES.targetDim.latencyBootByLoc(i), ES.targetDim.spdfByLoc(i,ES.targetDim.latencyBootInfoByLoc{i}.latencyTInd), ...
                '.', 'Color', cols(i,:) * 0.9, 'MarkerSize', 18, 'LineWidth', 2);
    end
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([0.28 0.28], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min rt

% shade in baseline and analysis windows
fillH = jbfill(ES.targetDimDelayWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.targetDimResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Target Dimming (s)');
% ylabel('Estimated Spike Rate (Hz)');

%% adjust all ylim
allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf) ylim(axTargetDimSpdf)];
yBounds = [min(allYBounds) max(1, max(allYBounds))];
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);

%% add response metrics text

redCol = [0.9 0.2 0.05];
greenCol = [0.05 0.5 0.05];

% % just do the most extreme comparison and correct for multiple comparisons
% if min(ES.cueResponseVsBaselineTTestStatsByLoc) < 0.05 / sum(~isnan(ES.cueResponsePValueByBootstrapBaselineSpdfByLoc))
%     cueResponsePValueByBootstrapBaselineSpdfCol = greenCol;
% else
%     cueResponsePValueByBootstrapBaselineSpdfCol = redCol;
% end
% 
% % just do the most extreme comparison and correct for multiple comparisons
% if min(ES.cueTargetDelayPValueByBootstrapBaselineSpdfByLoc) < 0.05 / sum(~isnan(ES.cueTargetDelayPValueByBootstrapBaselineSpdfByLoc))
%     cueTargetDelayPValueByBootstrapBaselineSpdfCol = greenCol;
% else
%     cueTargetDelayPValueByBootstrapBaselineSpdfCol = redCol;
% end
% 
% % just do the most extreme comparison and correct for multiple comparisons
% if min(ES.arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc) < 0.05 / sum(~isnan(ES.arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc))
%     arrayHoldResponsePValueByBootstrapBaselineSpdfCol = greenCol;
% else
%     arrayHoldResponsePValueByBootstrapBaselineSpdfCol = redCol;
% end
% 
% % just do the most extreme comparison and correct for multiple comparisons
% if min(ES.targetDimDelayPValueByBootstrapBaselineSpdfByLoc) < 0.05 / sum(~isnan(ES.targetDimDelayPValueByBootstrapBaselineSpdfByLoc))
%     targetDimDelayPValueByBootstrapBaselineSpdfCol = greenCol;
% else
%     targetDimDelayPValueByBootstrapBaselineSpdfCol = redCol;
% end
% 
% % just do the most extreme comparison and correct for multiple comparisons
% if min(ES.targetDimResponsePValueByBootstrapBaselineSpdfByLoc) < 0.05 / sum(~isnan(ES.targetDimResponsePValueByBootstrapBaselineSpdfByLoc))
%     targetDimResponsePValueByBootstrapBaselineSpdfCol = greenCol;
% else
%     targetDimResponsePValueByBootstrapBaselineSpdfCol = redCol;
% end

if ES.cueResponseVsBaselineSignRankTestStatsByLoc(ES.inRFLoc).p < 0.05 / numel(ES.unusedLocs)
    cueResponseVsBaselinePValueInRFCol = greenCol;
else
    cueResponseVsBaselinePValueInRFCol = redCol;
end

if ES.cueResponseVsBaselineSignRankTestStatsByLoc(ES.exRFLoc).p < 0.05 / numel(ES.unusedLocs)
    cueResponseVsBaselinePValueExRFCol = greenCol;
else
    cueResponseVsBaselinePValueExRFCol = redCol;
end

if ES.arrayHoldResponseVsCueTargetDelaySignRankTestStatsByLoc(ES.inRFLoc).p < 0.05 / numel(ES.unusedLocs)
    arrayHoldResponseVsBaselinePValueInRFCol = greenCol;
else
    arrayHoldResponseVsBaselinePValueInRFCol = redCol;
end

if ES.arrayHoldResponseVsCueTargetDelaySignRankTestStatsByLoc(ES.exRFLoc).p < 0.05 / numel(ES.unusedLocs)
    arrayHoldResponseVsBaselinePValueExRFCol = greenCol;
else
    arrayHoldResponseVsBaselinePValueExRFCol = redCol;
end

if ES.targetDimResponseVsTargetDimDelaySignRankTestStatsByLoc(ES.inRFLoc).p < 0.05 / numel(ES.unusedLocs)
    targetDimResponseVsBaselinePValueInRFCol = greenCol;
else
    targetDimResponseVsBaselinePValueInRFCol = redCol;
end

if ES.targetDimResponseVsTargetDimDelaySignRankTestStatsByLoc(ES.exRFLoc).p < 0.05 / numel(ES.unusedLocs)
    targetDimResponseVsBaselinePValueExRFCol = greenCol;
else
    targetDimResponseVsBaselinePValueExRFCol = redCol;
end

if ES.cueTargetDelayDiffPValueByShuffleSpdf < 0.05
    cueTargetDelayDiffPValueByShuffleSpdfCol = greenCol;
else
    cueTargetDelayDiffPValueByShuffleSpdfCol = redCol;
end

if ES.targetDimDelayDiffPValueByShuffleSpdf < 0.05
    targetDimDelayDiffPValueByShuffleSpdfCol = greenCol;
else
    targetDimDelayDiffPValueByShuffleSpdfCol = redCol;
end

if ES.cueTargetDelayAIPValueByShuffleSpdf < 0.05
    cueTargetDelayAICol = greenCol;
else
    cueTargetDelayAICol = redCol;
end

if ES.targetDimDelayAIPValueByShuffleSpdf < 0.05
    targetDimDelayAICol = greenCol;
else
    targetDimDelayAICol = redCol;
end

if ES.cueResponseInfoRateStruct.infoRatePValueByShuffleSpdf < 0.05
    cueResponseInfoRateCol = greenCol;
else
    cueResponseInfoRateCol = redCol;
end

if ES.cueTargetDelayInfoRateStruct.infoRatePValueByShuffleSpdf < 0.05
    cueTargetDelayInfoRateCol = greenCol;
else
    cueTargetDelayInfoRateCol = redCol;
end

if ES.arrayHoldResponseInfoRateStruct.infoRatePValueByShuffleSpdf < 0.05
    arrayHoldResponseInfoRateCol = greenCol;
else
    arrayHoldResponseInfoRateCol = redCol;
end

if ES.targetDimDelayInfoRateStruct.infoRatePValueByShuffleSpdf < 0.05
    targetDimDelayInfoRateCol = greenCol;
else
    targetDimDelayInfoRateCol = redCol;
end

if ES.targetDimResponseInfoRateStruct.infoRatePValueByShuffleSpdf < 0.05
    targetDimResponseInfoRateCol = greenCol;
else
    targetDimResponseInfoRateCol = redCol;
end

textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
text(axBig, -0.03, infoText2Top, {...
        sprintf('Blocks: %s', blockName), ...
        sprintf('Num Randomizations: %d', ES.numRandomizations), ...
        '', ...
        sprintf('{\\color[rgb]{%f,%f,%f}InRF (mean cue): %d} (peak: %d, extreme: %d)', ...
            cols(ES.inRFLoc,:), ES.inRFLoc, ES.inRFLocByMax, ES.inRFLocByMaxSD), ...
        sprintf('Baseline (SPDF): %0.1f Hz, Max: %0.1f Hz', ES.averageFiringRatesBySpdf.preCueBaseline.all, ES.maxFiringRateBySpdfInclMotor), ...
        sprintf('Cue {\\color[rgb]{%f,%f,%f}InRF}: %0.1f Hz, {\\color[rgb]{%f,%f,%f}ExRF}: %0.1f Hz', ...
            cols(ES.inRFLoc,:), ES.averageFiringRatesBySpdf.cueResponse.byLoc(ES.inRFLoc), ...
            cols(ES.exRFLoc,:), ES.averageFiringRatesBySpdf.cueResponse.byLoc(ES.exRFLoc)), ...
        '', ...
        sprintf('{\\color[rgb]{%f,%f,%f}InRF}: Cue Resp {\\neq} Baseline {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.inRFLoc,:), cueResponseVsBaselinePValueInRFCol, ES.cueResponseVsBaselineSignRankTestStatsByLoc(ES.inRFLoc).p), ...
        sprintf('{\\color[rgb]{%f,%f,%f}ExRF}: Cue Resp {\\neq} Baseline {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.exRFLoc,:), cueResponseVsBaselinePValueExRFCol, ES.cueResponseVsBaselineSignRankTestStatsByLoc(ES.exRFLoc).p), ...
        sprintf('{\\color[rgb]{%f,%f,%f}InRF}: Array Hold Resp {\\neq} CT Delay {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.inRFLoc,:), arrayHoldResponseVsBaselinePValueInRFCol, ES.arrayHoldResponseVsCueTargetDelaySignRankTestStatsByLoc(ES.inRFLoc).p), ...
        sprintf('{\\color[rgb]{%f,%f,%f}ExRF}: Array Hold Resp {\\neq} CT Delay {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.exRFLoc,:), arrayHoldResponseVsBaselinePValueExRFCol, ES.arrayHoldResponseVsCueTargetDelaySignRankTestStatsByLoc(ES.exRFLoc).p), ...
        sprintf('{\\color[rgb]{%f,%f,%f}InRF}: Dimming Resp {\\neq} TD Delay {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.inRFLoc,:), targetDimResponseVsBaselinePValueInRFCol, ES.targetDimResponseVsTargetDimDelaySignRankTestStatsByLoc(ES.inRFLoc).p), ...
        sprintf('{\\color[rgb]{%f,%f,%f}ExRF}: Dimming Resp {\\neq} TD Delay {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cols(ES.exRFLoc,:), targetDimResponseVsBaselinePValueExRFCol, ES.targetDimResponseVsTargetDimDelaySignRankTestStatsByLoc(ES.exRFLoc).p), ...
        '', ...
        sprintf('Cue Resp Lat: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.cueOnset.latencyBootByLoc(ES.exRFLoc) * 1000)), ...
        sprintf('Array Resp Lat: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.arrayOnset.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.arrayOnset.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
        sprintf('Dimming Resp Lat: {\\color[rgb]{%f,%f,%f}InRF}: %d ms, {\\color[rgb]{%f,%f,%f}ExRF}: %d ms', ...
            cols(ES.inRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.inRFLoc) * 1000), ...
            cols(ES.exRFLoc,:), round(ES.targetDim.latencyBootByLoc(ES.exRFLoc) * 1000)) ... 
        '', ...
        sprintf('Cue Resp Info: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cueResponseInfoRateCol, ES.cueResponseInfoRateStruct.infoRate, ...
            cueResponseInfoRateCol, ES.cueResponseInfoRateStruct.infoRatePValueByShuffleSpdf), ...
        sprintf('Cue-Target Delay Info: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cueTargetDelayInfoRateCol, ES.cueTargetDelayInfoRateStruct.infoRate, ...
            cueTargetDelayInfoRateCol, ES.cueTargetDelayInfoRateStruct.infoRatePValueByShuffleSpdf), ...
        sprintf('Array Hold Resp Info: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            arrayHoldResponseInfoRateCol, ES.arrayHoldResponseInfoRateStruct.infoRate, ...
            arrayHoldResponseInfoRateCol, ES.arrayHoldResponseInfoRateStruct.infoRatePValueByShuffleSpdf), ...
        sprintf('Target-Dim Delay Info: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            targetDimDelayInfoRateCol, ES.targetDimDelayInfoRateStruct.infoRate, ...
            targetDimDelayInfoRateCol, ES.targetDimDelayInfoRateStruct.infoRatePValueByShuffleSpdf), ...
        sprintf('Dimming Resp Info: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            targetDimResponseInfoRateCol, ES.targetDimResponseInfoRateStruct.infoRate, ...
            targetDimResponseInfoRateCol, ES.targetDimResponseInfoRateStruct.infoRatePValueByShuffleSpdf), ...
        '', ...
        sprintf('CT Delay Diff: {\\color[rgb]{%f,%f,%f}%0.1f Hz}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cueTargetDelayDiffPValueByShuffleSpdfCol, ES.cueTargetDelayDiff, ...
            cueTargetDelayDiffPValueByShuffleSpdfCol, ES.cueTargetDelayDiffPValueByShuffleSpdf), ...
        sprintf('CT Delay AI: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            cueTargetDelayAICol, ES.cueTargetDelayAI, cueTargetDelayAICol, ES.cueTargetDelayAIPValueByShuffleSpdf), ...
        sprintf('TD Delay Diff: {\\color[rgb]{%f,%f,%f}%0.1f Hz}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            targetDimDelayDiffPValueByShuffleSpdfCol, ES.targetDimDelayDiff, ...
            targetDimDelayDiffPValueByShuffleSpdfCol, ES.targetDimDelayDiffPValueByShuffleSpdf), ...
        sprintf('TD Delay AI: {\\color[rgb]{%f,%f,%f}%0.2f}, Shuf {\\color[rgb]{%f,%f,%f}P = %0.3f}', ...
            targetDimDelayAICol, ES.targetDimDelayAI, targetDimDelayAICol, ES.targetDimDelayAIPValueByShuffleSpdf), ...
        }, ...
        textParams{:});

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
