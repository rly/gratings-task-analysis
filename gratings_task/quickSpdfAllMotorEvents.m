function quickSpdfAllMotorEvents(ES, blockName, ...
        D, unitStructInd, unitStruct, nLoc, isZeroDistractors, plotFileName)

unitName = unitStruct.name;
nTrials = numel(ES.UE.cueOnset);
nValidTrials = numel(ES.cueOnset.validEventTimes);
cols = lines(6);

%%
f = figure_tr_inch(13, 7.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('Gratings Attention Task: %s (%d/%d trials)', unitName, nValidTrials, nTrials);
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
histH = 0.04;

waveformLeft1 = 0.05;
spdfEnterFixationLeft = waveformLeft1 + waveformW + 0.08;
spdfCueOnsetLeft = spdfEnterFixationLeft + spdfW + 0.045;
spdfTargetDimLeft = spdfEnterFixationLeft;
spdfExitFixationLeft = spdfCueOnsetLeft;

btm = 0.07;
infoTextTop = 0.02;
infoText2Top = btm + 0.64;
rasterBtm = btm + spdfH + 0.08;
spdfTopRowBtm = rasterBtm;
waveformBtm = 0.76;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D.allUnitStructs, unitStructInd);

%% info
writeUnitInfo(unitStruct, axBig, -0.03, infoTextTop, 1);

%% spdf for enter fixation
axEnterFixationSpdf = axes('Position', [spdfEnterFixationLeft spdfTopRowBtm spdfW spdfH]); 

t = ES.enterFixation.t - ES.enterFixation.window(1);
xBounds = [-0.6 0.6];
hold on;
for i = 1:nLoc
    if any(isnan(ES.enterFixation.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.enterFixation.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
end
origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([0.35 0.35], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
fillH = jbfill(ES.preEnterFixationWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.postEnterFixationLateWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Enter Fixation (s)');
ylabel('Estimated Spike Rate (Hz)');

axEnterFixationHist = axes('Position', [spdfEnterFixationLeft spdfTopRowBtm spdfW histH]); 
hold on;
histogram(ES.enterFixationToCueOnsetTime, xBounds(1):0.025:xBounds(2), 'FaceColor', [0.25 0.75 0]);
histogram(ES.enterFixationToLeverPressTime, xBounds(1):0.025:xBounds(2), 'FaceColor', [0.75 0.25 0.75]);

xlim(xBounds);
set(gca, 'Color', 'none');
set(gca, 'box', 'off');
set(gca, 'XTick', []);
set(gca, 'YTick', []);

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft spdfTopRowBtm spdfW spdfH]); 

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
% ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

t = ES.targetDimBal.t - ES.targetDimBal.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(nLoc, 1);
for i = 1:nLoc
    if any(isnan(ES.targetDimBal.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.targetDimBal.spdfByLoc(i,:),  ...
            'Color', cols(i,:), 'LineWidth', 2);
    legendEntry{i} = sprintf('{\\color[rgb]{%f,%f,%f}N = %d}', ...
            cols(i,:), numel(ES.targetDimBal.spikeTimesByLoc{i}));
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
xlabel('Time from Target Dim (s)');
ylabel('Estimated Spike Rate (Hz)');

text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

t = ES.exitFixation.t - ES.exitFixation.window(1);
xBounds = [-0.6 0.6];
hold on;
for i = 1:nLoc
    if any(isnan(ES.exitFixation.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.exitFixation.spdfByLoc(i,:),  ...
            'Color', cols(i,:), 'LineWidth', 2);
end
legendEntry = cell(2, 1);
plot(t, ES.exitFixationLeft.spdf,  ...
            'Color', cols(5,:), 'LineWidth', 2);
legendEntry{1} = sprintf('{\\color[rgb]{%f,%f,%f}Saccade Left N = %d}', ...
            cols(5,:), numel(ES.exitFixationLeft.spikeTimes));
plot(t, ES.exitFixationRight.spdf,  ...
            'Color', cols(6,:), 'LineWidth', 2);
legendEntry{2} = sprintf('{\\color[rgb]{%f,%f,%f}Saccade Right N = %d}', ...
            cols(6,:), numel(ES.exitFixationRight.spikeTimes));

origYLim = ylim();
plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));
plot([-0.28 -0.28], [0 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min rt

% shade in baseline and analysis windows
fillH = jbfill(ES.preExitFixationWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Exit Fixation (s)');
% ylabel('Estimated Spike Rate (Hz)');

text(0.72, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'All Trials', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

axExitFixationHist = axes('Position', [spdfExitFixationLeft btm spdfW histH]); 
hold on;
histogram(ES.exitFixationToJuiceEventTime, xBounds(1):0.025:xBounds(2), 'FaceColor', [1 0 0]);
histogram(ES.exitFixationToLeverReleaseTime, xBounds(1):0.025:xBounds(2), 'FaceColor', [0.75 0.25 0.75]);
histogram(-[ES.arrayOnsetRelToExitFixationTime; ES.targetDimToExitFixationTime], xBounds(1):0.025:xBounds(2), 'FaceColor', [1 0.5 0]);

xlim(xBounds);
set(gca, 'Color', 'none');
set(gca, 'box', 'off');
set(gca, 'XTick', []);
set(gca, 'YTick', []);


%% adjust all ylim
allYBounds = [ylim(axEnterFixationSpdf) ylim(axCueOnsetSpdf) ylim(axTargetDimSpdf) ylim(axExitFixationSpdf)];
yBounds = [min(allYBounds) max(1, max(allYBounds))];
ylim(axEnterFixationSpdf, yBounds);
ylim(axCueOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);
ylim(axExitFixationSpdf, yBounds);

%% add response metrics text

textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
text(axBig, -0.03, infoText2Top, {...
        sprintf('Blocks: %s', blockName) ...
        }, ...
        textParams{:});

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
