function [f1, ax1, ax2, ax3] = plotMetricDiff(mInRF, mExRF, isDPul, isVPul, histBinStep)
% works for firing rate and fano factor and whatever else

%% stat test
mDiff = mInRF - mExRF;

meanMDiff = mean(mDiff);
medianMDiff = median(mDiff);
p = signrank(mDiff);
fprintf('\tAll: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d\n', ...
        meanMDiff, medianMDiff, p, numel(mDiff));
    
mDiffDPul = mDiff(isDPul);
meanMDiff = mean(mDiffDPul);
medianMDiff = median(mDiffDPul);
p = signrank(mDiffDPul);
fprintf('\tDPul: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d\n', ...
        meanMDiff, medianMDiff, p, numel(mDiffDPul));

mDiffVPul = mDiff(isVPul);
meanMDiff = mean(mDiffVPul);
medianMDiff = median(mDiffVPul);
p = signrank(mDiffVPul);
fprintf('\tVPul: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d\n', ...
        meanMDiff, medianMDiff, p, numel(mDiffVPul));

%% plot parameters
frBounds = [0 max([max(mInRF) max(mExRF)])];
maxAbsDiffFR = max(abs(mDiff));

histXBounds = [-ceil(maxAbsDiffFR / histBinStep) ceil(maxAbsDiffFR / histBinStep)] * histBinStep;
histBinEdges = histXBounds(1):histBinStep:histXBounds(2);

cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

%% plot
f1 = figure_tr_inch(13, 4.5);

%% scatter plot
ax1 = subaxis(1, 3, 1);
hold on;
plot(frBounds, frBounds, 'Color', 0.3*ones(3, 1)); 
plotParams = {'MarkerSize', 20};
plot(mInRF, mExRF, '.', plotParams{:}); % plot under as indexing sanity check
h1 = plot(mInRF(isDPul), mExRF(isDPul), '.', 'Color', dPulCol, plotParams{:});
h2 = plot(mInRF(isVPul), mExRF(isVPul), '.', 'Color', vPulCol, plotParams{:});
axis equal;
xlim(frBounds);
ylim(frBounds);
box off;
title('All Units');
legend([h1 h2], {'dPul', 'vPul'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

%% histogram of differences dorsal pulvinar
ax2 = subaxis(1, 3, 2); 
hold on;
histH = histogram(mDiffDPul, histBinEdges);
histH.FaceColor = dPulCol;
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
ylabel('Number of Units');
title(sprintf('Dorsal Pulvinar (N=%d)', sum(isDPul)));
box off;
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

%% histogram of differences ventral pulvinar
ax3 = subaxis(1, 3, 3); 
hold on;
histH = histogram(mDiffVPul, histBinEdges);
histH.FaceColor = vPulCol;
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
ylabel('Number of Units');
title(sprintf('Ventral Pulvinar (N=%d)', sum(isVPul)));
box off;
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

