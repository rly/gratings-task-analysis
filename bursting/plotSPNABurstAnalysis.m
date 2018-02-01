function plotSPNABurstAnalysis(spikeTimesInWindow, spikeFs, lagT, maxLag, ...
        meanAutoCorr, meanCrossCorr, stdCrossCorr, SPNA, BRI, figTitle, saveFile)


figure_tr_inch(18, 9);
set(gcf,'Color','white');
set(gcf,'renderer','painters');

% raster plot by time
subplot(1, 4, 1);
hold on;
rasterY = 0;
plotParams = {'Color', [0 0 0], 'MarkerSize', 1.5, 'MarkerFaceColor', [0 0 0]};

% plot diamond at each spike, one row (y-coord) per trial
for j = 1:numel(spikeTimesInWindow)
    rasterY = rasterY + 1;
    if ~isempty(spikeTimesInWindow{j})
        plot(spikeTimesInWindow{j}, ...
                rasterY*ones(1, length(spikeTimesInWindow{j})),...
                'd', plotParams{:});
    end
end
title('Raster plot by time');
xlabel('Time (ms)');
ylabel('Trial number');
ylim([0 rasterY + 1]);
set(gca, 'FontSize', 16);

subplot(1, 4, 2);
diffSpikeTimesInTrial = cellfun(@diff, spikeTimesInWindow, 'UniformOutput', 0);
diffSpikeTimesInTrialAll = vertcat(diffSpikeTimesInTrial{:});
histogram(diffSpikeTimesInTrialAll, [10.^(0:0.1:1) 10.^(1.1:0.1:3)]);
xlim([1 10^3]);
set(gca, 'XScale', 'log') 
set(gca, 'XTick', 10.^(0:3));
set(gca, 'XTickLabel', 10.^(0:3));
xlabel('ISI (ms)');
ylabel('Frequency');
set(gca, 'FontSize', 16);

subplot(1, 4, 3);
hold on;
plot(lagT, meanAutoCorr(maxLag+2:end), '-o', 'LineWidth', 3);
plot(lagT, meanCrossCorr(maxLag+2:end), '-r', 'LineWidth', 3);
plot(lagT, meanCrossCorr(maxLag+2:end) + stdCrossCorr(maxLag+2:end), '--r', 'LineWidth', 3);
title('Auto-correlation');
xlabel('Lag (ms)');
ylabel('Correlation units');
xlim([1 maxLag]);
legend({'Auto-correlation', 'Mean cross-correlation', 'Mean + SD cross-correlation'});
set(gca, 'FontSize', 16);

subplot(1, 4, 4);
hold on;
plot(lagT, SPNA(maxLag+2:end), '-o', 'LineWidth', 3);
xlabel('Lag (ms)');
ylabel('Shuffle Predictor SD Units');
title('Shuffle Predictor Normalized Autocorrelation');
xlim([1 maxLag]);

% shade and show BRI over first 4 ms
fillT = 1:1000/spikeFs:4;
briThresh = 1;
plot(fillT, ones(size(fillT)) * BRI, '-g', 'LineWidth', 3);
plot(fillT, ones(size(fillT)) * briThresh, '-m', 'LineWidth', 3);

currYLim = ylim();
jbfill(fillT, ones(size(fillT)) * currYLim(1), ...
        ones(size(fillT)) * currYLim(2), ...
        [0 1 0], [0 1 0], 1, 0.2);
ylim(currYLim);

set(gca, 'FontSize', 16);

legend({'SPNA', sprintf('BRI = %0.1f', BRI), sprintf('BRI threshold = %0.1f', briThresh)});

%% title over all subplots
suptitle(figTitle, 'FontSize', 22);

%% save
if nargin >= 11
	export_fig(saveFile, '-nocrop'); %, '-r300');
end