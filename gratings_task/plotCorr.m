function plotCorr(xInRF, xExRF, yInRF, yExRF, isDPul, isVPul)

xDiff = xInRF - xExRF;
yDiff = yInRF - yExRF;

[r, p] = corr(xDiff, yDiff, 'type', 'Spearman');
fprintf('\tAll: Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, numel(xDiff));

[r, p] = corr(xDiff(isDPul), yDiff(isDPul), 'type', 'Spearman');
fprintf('\tDPul: Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, numel(xDiff(isDPul)));
    
[r, p] = corr(xDiff(isVPul), yDiff(isVPul), 'type', 'Spearman');
fprintf('\tVPul: Spearman rho = %0.2f, p = %0.5f, N = %d\n', ...
        r, p, numel(xDiff(isVPul)));
%% plot parameters
maxAbsDiffX = max(abs(xDiff));
xDiffBounds = maxAbsDiffX * [-1 1];

maxAbsDiffY = max(abs(yDiff));
yDiffBounds = maxAbsDiffY * [-1 1];

maxDiffBounds = max([maxAbsDiffX maxAbsDiffY]) * [-1 1];

cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

%% plot
f1 = figure_tr_inch(5, 5);
hold on;
plot(maxDiffBounds, maxDiffBounds, 'Color', 0.3*ones(3, 1)); 
plot(maxDiffBounds, [0 0], 'Color', zeros(3, 1)); 
plot([0 0], maxDiffBounds, 'Color', zeros(3, 1)); 
plotParams = {'MarkerSize', 20};
plot(xDiff, yDiff, '.', 'MarkerSize', 20, plotParams{:});
h1 = plot(xDiff(isDPul), yDiff(isDPul), '.', 'Color', dPulCol, plotParams{:});
h2 = plot(xDiff(isVPul), yDiff(isVPul), '.', 'Color', vPulCol, plotParams{:});
xlim(xDiffBounds);
ylim(yDiffBounds);
legend([h1 h2], {'dPul', 'vPul'}, 'Location', 'SouthEast');

set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);