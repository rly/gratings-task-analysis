function quickSpdfAllVisualEvents(ES, blockName, ...
        D, spikeStructInd, spikeStruct, nLoc, isZeroDistractors, plotFileName)

unitName = spikeStruct.name;
nTrials = numel(ES.UE.cueOnset);
cols = lines(4);

%%
figure_tr_inch(13, 7.5); clf;
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
arrayOnsetRelSpdfW = 0.354;
arrayOnsetHoldSpdfW = 0.354;
waveformH = 0.15;
rasterH = 0.39;
spdfH = 0.39;

waveformLeft1 = 0.05;
rasterCueOnsetLeft = waveformLeft1 + waveformW + 0.08;
rasterArrayOnsetRelLeft = rasterCueOnsetLeft + rasterW + 0.03;
rasterArrayOnsetHoldLeft = rasterArrayOnsetRelLeft + rasterW + 0.03;
spdfCueOnsetLeft = rasterArrayOnsetHoldLeft + rasterW + 0.045;
spdfArrayOnsetRelLeft = rasterCueOnsetLeft;
spdfArrayOnsetHoldLeft = spdfCueOnsetLeft;

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

data = ES.cueOnsetError.spikeTimes;
lineParams = {'Color', [1 0 0], 'LineWidth', 1};
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

%% plot raster aligned to array - release balanced
axes('Position', [rasterArrayOnsetRelLeft rasterBtm rasterW rasterH]); 
hold on;

xBounds = [-0.25 0.25];
window = ES.arrayOnsetRelBal.window;
data = ES.arrayOnsetRelBal.spikeTimes; % TODO show release trials
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

data = ES.arrayOnsetRelBalError.spikeTimes; % TODO show release trials
lineParams = {'Color', [1 0 0], 'LineWidth', 1};
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

%% plot raster aligned to array - hold balanced
axes('Position', [rasterArrayOnsetHoldLeft rasterBtm rasterW rasterH]); 
hold on;

xBounds = [-0.25 0.25];
window = ES.arrayOnsetHoldBal.window;
data = ES.arrayOnsetHoldBal.spikeTimes; % TODO show release trials
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

data = ES.arrayOnsetHoldBalError.spikeTimes; % TODO show release trials
lineParams = {'Color', [1 0 0], 'LineWidth', 1};
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

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft spdfCueOnsetBtm cueOnsetSpdfW spdfH]); 

t = ES.cueOnset.t - ES.cueOnset.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(1, 1); % don't preallocate
legendCount = 0;
for i = 1:nLoc
    if any(isnan(ES.cueOnset.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.cueOnset.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~any(isnan(ES.cueOnset.spdfByLoc(i,:)))
        plot(t, ES.cueOnsetError.spdfByLoc(i,:), '--', ...
                'Color', cols(i,:), 'LineWidth', 2);
    end
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Correct N = %d}', ...
            cols(i,:), numel(ES.cueOnset.spikeTimesByLoc{i}));
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Error N = %d}', ...
            cols(i,:), numel(ES.cueOnsetError.spikeTimesByLoc{i}));
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

text(0.75, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset - release balanced
axArrayOnsetRelSpdf = axes('Position', [spdfArrayOnsetRelLeft btm arrayOnsetRelSpdfW spdfH]); 

t = ES.arrayOnsetRelBal.t - ES.arrayOnsetRelBal.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(1, 1); % don't preallocate
legendCount = 0;
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetRelBal.spdfByLoc(i,:))) % TODO show release trials
        continue;
    end
    plot(t, ES.arrayOnsetRelBal.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~any(isnan(ES.arrayOnsetRelBalError.spdfByLoc(i,:)))
        plot(t, ES.arrayOnsetRelBalError.spdfByLoc(i,:), '--', ...
                'Color', cols(i,:), 'LineWidth', 2);
    end
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Correct N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetRelBal.spikeTimesByLoc{i}));
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Error N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetRelBalError.spikeTimesByLoc{i}));
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

text(0.75, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'Release Trials Balanced Only', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset - hold balanced
axArrayOnsetHoldSpdf = axes('Position', [spdfArrayOnsetHoldLeft btm arrayOnsetHoldSpdfW spdfH]); 

t = ES.arrayOnsetHoldBal.t - ES.arrayOnsetHoldBal.window(1);
xBounds = [-0.6 0.6];
hold on;
legendEntry = cell(1, 1); % don't preallocate
legendCount = 0;
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetHoldBal.spdfByLoc(i,:)))
        continue;
    end
    plot(t, ES.arrayOnsetHoldBal.spdfByLoc(i,:), ...
            'Color', cols(i,:), 'LineWidth', 2);
    if ~any(isnan(ES.arrayOnsetHoldBal.spdfByLoc(i,:)))
        plot(t, ES.arrayOnsetHoldBalError.spdfByLoc(i,:), '--', ...
                'Color', cols(i,:), 'LineWidth', 2);
    end
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Correct N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetHoldBal.spikeTimesByLoc{i}));
    legendCount = legendCount + 1;
    legendEntry{legendCount} = sprintf('{\\color[rgb]{%f,%f,%f}Error N = %d}', ...
            cols(i,:), numel(ES.arrayOnsetHoldBalError.spikeTimesByLoc{i}));
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

text(0.75, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);
text(0.02, 1, 'Hold Trials Balanced Only', 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% adjust all ylim
allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetRelSpdf) ylim(axArrayOnsetHoldSpdf)];
yBounds = [min(allYBounds) max(1, max(allYBounds))];
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetRelSpdf, yBounds);
ylim(axArrayOnsetHoldSpdf, yBounds);

%% add response metrics text

textParams = {'Units', 'normalized', 'FontSize', 8, 'VerticalAlignment', 'top'};
text(axBig, -0.03, infoText2Top, {}, ...
        textParams{:});

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
