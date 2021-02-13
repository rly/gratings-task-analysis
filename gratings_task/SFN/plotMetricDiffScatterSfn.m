function ax = plotMetricDiffScatterSfn(mInRF, mExRF, col, isSigUnit)

%% plot
f1 = figure_tr_inch(5, 4.5);

%% scatter plot of InRF fano factor vs ExRF
ax = subaxis(1, 1, 1, 'MB', 0.17, 'MT', 0.05, 'ML', 0.15); 
hold on;
plot(mBounds, mBounds, 'Color', 'k', 'LineWidth', 2); 
plotParams = {'MarkerSize', 6, 'LineWidth', 1.5};
plot(mInRF, mExRF, 'o', 'MarkerFaceColor', 0.8*ones(3, 1), 'MarkerEdgeColor', col, plotParams{:});
plot(mInRF(isSigUnit), mExRF(isSigUnit), 'o', 'MarkerFaceColor', zeros(3, 1), 'MarkerEdgeColor', col, plotParams{:});
axis equal;
xlim(mBounds);
ylim(mBounds);
box off;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');
