function quickSpdf(spikeTimesByLoc, nLoc, t, periEventWindow, kernelSigma)

hold on;

spdfByLoc = cell(nLoc,1);
legendEntry = cell(4,1);
fHandles = nan(nLoc, 1);
for i = 1:nLoc
    spdfByLoc{i} = edpsth_notranspose(spikeTimesByLoc{i}, kernelSigma, 'n', [], 0, t);
    if isempty(spdfByLoc{i}) || any(isnan(spdfByLoc{i}))
        continue;
    end
    legendEntry{i} = sprintf('P%d (N=%d)', i, numel(spikeTimesByLoc{i}));
    fHandles(i) = plot(t - periEventWindow(1), spdfByLoc{i}, 'LineWidth', 2);
end
legendEntry(isnan(fHandles)) = [];
fHandles(isnan(fHandles)) = [];

plot([0 0], [0 100], '-', 'Color', 0.3 * ones(3, 1));
ylabel('Estimated Firing Rate (Hz)');
legend(fHandles, legendEntry, 'Location', 'NorthWest');
set(gcf, 'Color', 'white');
set(gca, 'FontSize', 16);