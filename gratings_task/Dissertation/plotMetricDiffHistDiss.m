function ax = plotMetricDiffHistDiss(mInRF, mExRF, col, isSigUnit, histBinStep)
% works for firing rate and fano factor and whatever else

%% stat test
mDiff = mInRF - mExRF;

meanMDiff = mean(mDiff);
medianMDiff = median(mDiff);
p = signrank(mDiff);
fprintf('\tAll: Mean diff = %0.2f, median diff = %0.2f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
        meanMDiff, medianMDiff, p, numel(mDiff), round(sum(isSigUnit) / numel(isSigUnit) * 100));

%% plot parameters
maxAbsDiffFR = max(abs(mDiff));

histXBounds = [-1 1];%[-ceil(maxAbsDiffFR / histBinStep) ceil(maxAbsDiffFR / histBinStep)] * histBinStep;
histBinEdges = -1:0.1:1;%histXBounds(1):histBinStep:histXBounds(2);

%% plot
f1 = figure_tr_inch(5, 4);

%% histogram of differences ventral pulvinar
ax = subaxis(1, 1, 1, 'MB', 0.2, 'MT', 0.05, 'ML', 0.15); 
hold on;
histH = histogram(mDiff, histBinEdges);
histH.FaceColor = col;
histH = histogram(mDiff(isSigUnit), histBinEdges);
histH.FaceColor = zeros(3, 1); % overlay (def alpha is 0.6)

origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
ylabel('Number of Units');
box off;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
text(gca, 1, 1, sprintf('N=%d', numel(mDiff)), 'Units', 'normalized', 'FontSize', 14, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right'); 

if origYLim(2) == 2
    set(gca, 'YTick', 0:2);
end
