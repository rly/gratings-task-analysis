function ax = plotMetricDiffHist(mInRF, mExRF, col, isSigUnit, histBinStep)
% works for firing rate and fano factor and whatever else

%% stat test
mDiff = mInRF - mExRF;

meanMDiff = mean(mDiff);
medianMDiff = median(mDiff);
p = signrank(mDiff);
fprintf('\tAll: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
        meanMDiff, medianMDiff, p, numel(mDiff), round(sum(isSigUnit) / numel(isSigUnit) * 100));

%% plot parameters
maxAbsDiffFR = max(abs(mDiff));

histXBounds = [-ceil(maxAbsDiffFR / histBinStep) ceil(maxAbsDiffFR / histBinStep)] * histBinStep;
histBinEdges = histXBounds(1):histBinStep:histXBounds(2);

%% plot
f1 = figure_tr_inch(5, 5);

%% histogram of differences ventral pulvinar
ax = subaxis(1, 1, 1, 'MB', 0.18, 'MT', 0.05); 
hold on;
histH = histogram(mDiff, histBinEdges);
histH.FaceColor = col;
histH = histogram(mDiff(isSigUnit), histBinEdges);
histH.FaceColor = zeros(3, 1);

origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
ylabel('Number of Units');
box off;
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

