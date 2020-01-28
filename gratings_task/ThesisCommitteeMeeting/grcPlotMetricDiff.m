function [f1,ax1] = grcPlotMetricDiff(mInRF, mExRF, isDPul, isVPul, isSigUnit)
% works for firing rate and fano factor and whatever else

%% stat test
mDiff = mInRF - mExRF;

meanMDiff = mean(mDiff);
medianMDiff = median(mDiff);
p = signrank(mDiff);
fprintf('\tAll: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
        meanMDiff, medianMDiff, p, numel(mDiff), round(sum(isSigUnit) / numel(isSigUnit) * 100));

if any(isDPul)
    mDiffDPul = mDiff(isDPul);
    meanMDiff = mean(mDiffDPul);
    medianMDiff = median(mDiffDPul);
    p = signrank(mDiffDPul);
    fprintf('\tDPul: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
            meanMDiff, medianMDiff, p, numel(mDiffDPul), round(sum(isSigUnit(isDPul)) / numel(isDPul) * 100));
end

if any(isVPul)
    mDiffVPul = mDiff(isVPul);
    meanMDiff = mean(mDiffVPul);
    medianMDiff = median(mDiffVPul);
    p = signrank(mDiffVPul);
    fprintf('\tVPul: Mean diff = %0.1f, median diff = %0.1f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
            meanMDiff, medianMDiff, p, numel(mDiffVPul), round(sum(isSigUnit(isVPul)) / numel(isVPul) * 100));
end

%% plot parameters
mBounds = [min([min(mInRF) min(mExRF)]) max([max(mInRF) max(mExRF)])];
if mBounds(1) > 0
    mBounds(1) = 0;
end

cols = lines(6);
dPulCol = cols(4,:);
vPulCol = cols(5,:);

%% plot
f1 = figure_tr_inch(5, 5);

%% scatter plot
ax1 = subaxis(1, 1, 1, 'MB', 0.15);
hold on;
plot(mBounds, mBounds, 'Color', 0.3*ones(3, 1)); 
plotParams = {'MarkerSize', 6, 'LineWidth', 1};
plot(mInRF, mExRF, '.', plotParams{:}); % plot under as indexing sanity check
h1 = plot(mInRF(isDPul), mExRF(isDPul), 'o', 'MarkerFaceColor', 0.8*ones(3, 1), 'MarkerEdgeColor', dPulCol, plotParams{:});
h2 = plot(mInRF(isVPul), mExRF(isVPul), 'o', 'MarkerFaceColor', 0.8*ones(3, 1), 'MarkerEdgeColor', vPulCol, plotParams{:});
plot(mInRF(isDPul & isSigUnit), mExRF(isDPul & isSigUnit), 'o', 'MarkerFaceColor', zeros(3, 1), 'MarkerEdgeColor', dPulCol, plotParams{:});
plot(mInRF(isVPul & isSigUnit), mExRF(isVPul & isSigUnit), 'o', 'MarkerFaceColor', zeros(3, 1), 'MarkerEdgeColor', vPulCol, plotParams{:});
axis equal;
xlim(mBounds);
ylim(mBounds);
box off;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
