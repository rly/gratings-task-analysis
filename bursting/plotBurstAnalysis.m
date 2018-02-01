function plotBurstAnalysis(silentPeriodDurationsPrecedingSpike, numSpikeTimes, ...
        spikesInBurst, diffSpikeTimes, Fs, minSilentPeriod, minISIWithinBurst, ...
        spikeStructInd, allSpikeStructs, figTitle, saveFile)

figure_tr_inch(14, 8); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

minPlotPowTen = 0;
maxPlotPowTen = 4; % 10000

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.90], 'Visible', 'off');
set(get(axBig,'Title'), 'Visible', 'on')

titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(figTitle, 'FontSize', 18, titleParams{:});

%% location params
waveformW = 0.1;
rasterByTimeW = 0.1;
isiHistW = 0.2;
isiIsiW = 0.6;
% miniHeatW = 0.18;
waveformH = 0.15;
rasterByTimeH = 0.33;
isiHistH = 0.83;
isiIsiH = 0.83;
% miniHeatH = 0.385;

waveformLeft1 = 0.05;
rasterByTimeLeft1 = waveformLeft1;
isiHistLeft1 = rasterByTimeLeft1 + rasterByTimeW + 0.06;
isiIsiLeft = isiHistLeft1 + isiHistW + 0.02;
% miniHeatLeft1 = heatLeft + heatW + 0.005;
% miniHeatLeft2 = miniHeatLeft1 + miniHeatW + 0.01;
btm = 0.08;
rasterByTimeBtm = btm + 0.25;
waveformBtm = rasterByTimeBtm + rasterByTimeH + 0.1;
% miniHeatBtm2 = btm + miniHeatH + 0.08;

%% plot spike waveform
axes('Position', [waveformLeft1 waveformBtm waveformW waveformH]); 
hold on;

spikeStruct = allSpikeStructs{spikeStructInd};

nWfTime = numel(spikeStruct.meanWf);
spikeFs = 40000;
waveformT = (1:nWfTime)/(spikeFs/1000);
yBounds = [-0.08 0.08];
plot([waveformT(1) waveformT(end)], [0 0], 'Color', 0.5*ones(3, 1));
plot([waveformT(1) waveformT(end)], [spikeStruct.threshold spikeStruct.threshold], '--', 'Color', [0.5 0.2 0.5]);
plot([spikeStruct.thresholdTime spikeStruct.thresholdTime]*1000, yBounds, '--', 'Color', [0.5 0.2 0.5]);
allUnitsThisChannel = findAllUnitsSameCh(allSpikeStructs, spikeStructInd);
allCols = lines(6);
for k = 1:numel(allUnitsThisChannel)
    unitInd = allUnitsThisChannel(k);
    if unitInd == spikeStructInd % the current unit
        currentUnitCol = allCols(k,:);
    end
    plot(waveformT, allSpikeStructs{unitInd}.meanWf, 'LineWidth', 2, 'Color', allCols(k,:));
end
plot(waveformT, allSpikeStructs{spikeStructInd}.meanWf, 'LineWidth', 4, 'Color', currentUnitCol);
sd = std(spikeStruct.wf);
jbfill(waveformT, spikeStruct.meanWf + sd, spikeStruct.meanWf - sd, currentUnitCol, ones(3, 1), 1, 0.5);
xlabel('Time (ms)');
xlim([0 waveformT(end)]);
set(gca, 'XTick', [0 waveformT(end)/2 waveformT(end)]);
ylim(yBounds);
title('Waveform');

%% plot raster
axes('Position', [rasterByTimeLeft1 rasterByTimeBtm rasterByTimeW rasterByTimeH]); 
hold on;

window = [0.25 0.25]; % secs before, secs after
% allAlignedSpikeTimes = createdatamatpt(spikeTimes.ts, flashOnsets, window);
% 
% data = allAlignedSpikeTimes;
% rasterY = 0;
% % plotParams = {'Color', [0 0 0], 'MarkerSize', 1, 'MarkerFaceColor', [0 0 0]};
% lineParams = {'Color', [0 0 0], 'LineWidth', 1};
% lineHeight = 5;
% % plot square at each spike, one row (y-coord) per trial
% for k = 1:numel(data)
%     rasterY = rasterY + 1;
%     if ~isempty(data(k).times)
%         adjTimes = data(k).times - window(1);
%         for l = 1:numel(data(k).times)
%             plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
%         end
% %         plot(data(k).times - window(1), rasterY*ones(1, length(data(k))),...
% %                 'd', plotParams{:});
%     end
% end


% yBounds = [0 rasterY+1];
% plot([0 0], yBounds, 'k');
% plot([analysisWindowOffset(1) analysisWindowOffset(1)], yBounds, '-', 'Color', [0 0.7 0]);
% plot([analysisWindowOffset(2) analysisWindowOffset(2)], yBounds, '-', 'Color', [0 0.7 0]);
% xlim([-0.1 0.25]);
% ylim(yBounds);
% set(gca, 'YDir', 'reverse');
% xlabel('Time from flash onset (s)'); 
% ylabel('Flash number');
% titleH = title('Raster by time', 'Interpreter', 'None');
% set(gca, 'Layer', 'top');


%% info
fsClassThresh = 0.35/1000;
if startsWith(spikeStruct.inflectionPattern, 'tp')
    if spikeStruct.troughToPeakTime <= fsClassThresh
        fsClass = 'Narrow';
    else
        fsClass = 'Broad';
    end
else
    fsClass = '';
end

textParams = {'Units', 'normalized', 'FontSize', 8, 'Interpreter', 'none'};
text(axBig, -0.03, 0.11, {spikeStruct.name, ...
        sprintf('Session Name: %s', spikeStruct.sessionName), ...
        sprintf('Channel ID: %d', spikeStruct.channelID), ...
        sprintf('Num Spikes: %d', numel(spikeStruct.ts)), ...
        sprintf('Threshold: %0.2f', spikeStruct.threshold), ...
        sprintf('Mean SD: %0.3f', mean(std(spikeStruct.wf))), ...
        sprintf('Trough to Peak Time: %0.2f ms', spikeStruct.troughToPeakTime*1000), ...
        sprintf('Num Inflections: %d', spikeStruct.numInflections), ...
        sprintf('Inflection Pattern: %s', spikeStruct.inflectionPattern), ...
        sprintf('Putative classification: %s', fsClass), ...
        sprintf('SPKC High-pass Filter: 300 Hz')}, ...
        textParams{:});

%% plot isi x isi distribution
axes('Position', [isiHistLeft1 btm isiHistW isiHistH]); 

histogram(diffSpikeTimes * Fs, 10.^(minPlotPowTen:0.1:maxPlotPowTen));
xlim(10.^[minPlotPowTen maxPlotPowTen]);
set(gca, 'XScale', 'log') 
set(gca, 'XTick', 10.^(minPlotPowTen:maxPlotPowTen));
set(gca, 'XTickLabel', 10.^(minPlotPowTen:maxPlotPowTen));
xlabel('ISI (ms)');
ylabel('Frequency');
origYLim = ylim;
hold on;
set(gca, 'FontSize', 12);
% title(sprintf('%d%% (%d/%d) of spikes are \npreceded by a silent period > %0.2fs', ...
%         round(numel(silentPeriodDurationsPrecedingSpike) / numSpikeTimes * 100), ...
%         numel(silentPeriodDurationsPrecedingSpike), numSpikeTimes, minSilentPeriod));

% shade in a green area to represent spikes that start a burst
xFillStartBurst = (minSilentPeriod * Fs):10^4;
jbfill(xFillStartBurst, zeros(size(xFillStartBurst)), ...
        ones(size(xFillStartBurst)) * origYLim(2), ...
        [0 1 0], [0 1 0], 1, 0.2);

    
% shade in an orange area to represent spikes that may be within a burst
xFillWithinBurst = 1:(minISIWithinBurst * Fs);
jbfill(xFillWithinBurst, zeros(size(xFillWithinBurst)), ...
        ones(size(xFillWithinBurst)) * origYLim(2), ...
        [1 0.4 0], [1 0.4 0], 1, 0.3);
    

% subplot(1,4,2);
% if ~isempty(spikesInBurst)
%     histogram(spikesInBurst, 2:6);
%     xlabel('Number of Spikes in a Burst');
%     ylabel('Frequency');
%     set(gca, 'FontSize', 16);
% end

%% plot previous ISI vs following ISI, without first or last spike
isiIsiH = axes('Position', [isiIsiLeft btm isiIsiW isiIsiH]);

loglog(diffSpikeTimes(1:end-1) * Fs, diffSpikeTimes(2:end) * Fs, '.');
axis(isiIsiH, 'square');
xlim(10.^[minPlotPowTen maxPlotPowTen]);
ylim(10.^[minPlotPowTen maxPlotPowTen]);
xlabel('Previous ISI (ms)');
ylabel('Following ISI (ms)');
set(gca, 'XTickLabel', 10.^(minPlotPowTen:maxPlotPowTen));
set(gca, 'YTickLabel', 10.^(minPlotPowTen:maxPlotPowTen));
set(gca, 'FontSize', 12);

% shade in a green area to represent spikes that start a burst
xFillStartBurst = (minSilentPeriod * Fs):10^maxPlotPowTen;
jbfill(xFillStartBurst, ones(size(xFillStartBurst)), ...
        ones(size(xFillStartBurst)) * minISIWithinBurst * Fs, ...
        [0 1 0], [0 1 0], 1, 0.2);

% shade in an orange area to represent spikes that may be within a burst
xFillWithinBurst = 1:(minISIWithinBurst * Fs);
jbfill(xFillWithinBurst, ones(size(xFillWithinBurst)), ...
        ones(size(xFillWithinBurst)) * minISIWithinBurst * Fs, ...
        [1 0.4 0], [1 0.4 0], 1, 0.3);
    
%% title over all subplots
% suptitle(figTitle, 'FontSize', 22);

%% save
if nargin >= 11
	export_fig(saveFile, '-nocrop'); %, '-r300');
end