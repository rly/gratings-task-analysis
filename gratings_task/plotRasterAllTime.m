function plotRasterAllTime(ES, D, unitStructInd, unitStruct, isZeroDistractors, plotFileName)

unitName = unitStruct.name;
nTrials = numel(ES.UE.cueOnset);
nValidTrials = numel(ES.cueOnset.validEventTimes);

%%
figure_tr_inch(13, 7.5); clf;
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
rasterW = 0.77;
waveformH = 0.15;
rasterH = 0.86;

waveformLeft1 = 0.05;
rasterLeft = waveformLeft1 + waveformW + 0.05;

btm = 0.07;
infoTextTop = 0.02;
infoText2Top = btm + 0.64;
rasterBtm = btm;
waveformBtm = 0.76;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
plotSpikeWaveform(D.allUnitStructs, unitStructInd);

%% info
writeUnitInfo(unitStruct, axBig, -0.03, infoTextTop, 1);

%% plot raster aligned to cue
axes('Position', [rasterLeft rasterBtm rasterW rasterH]); 
hold on;

%% plot 20 seconds per line
window = [0 20];
rowStartTimes = ES.spikeTs(1):window(2):ES.spikeTs(end);
% plot out of bounds data
spikeTsOutBounds = ES.spikeTs < ES.startTime | ES.spikeTs >= ES.endTime;
data = createnonemptydatamatpt(ES.spikeTs(spikeTsOutBounds), rowStartTimes, window);
rasterY = 0;
lineParams = {'Color', [1 0 0], 'LineWidth', 1};
lineHeight = 2;
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

xlim(window);
ylim(yBounds);
set(gca, 'YDir', 'reverse');
set(gca, 'XTickLabel', window(1):1:window(2));
xlabel('Time (s)');
ylabel('');

% plot in bounds data
spikeTsInBounds = ES.spikeTs >= ES.startTime & ES.spikeTs < ES.endTime & ~ES.isSpikeTsInTrial;
data = createnonemptydatamatpt(ES.spikeTs(spikeTsInBounds), rowStartTimes, window);
rasterY = 0;
lineParams = {'Color', [0 0 0], 'LineWidth', 1};
lineHeight = 2;
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

% plot trial spikes
spikeTsInBounds = ES.spikeTs >= ES.startTime & ES.spikeTs < ES.endTime & ES.isSpikeTsInTrial;
data = createnonemptydatamatpt(ES.spikeTs(spikeTsInBounds), rowStartTimes, window);
rasterY = 0;
lineParams = {'Color', [0 1 0.5], 'LineWidth', 1};
lineHeight = 2;
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

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
