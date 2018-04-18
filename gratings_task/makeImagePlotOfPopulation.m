function makeImagePlotOfPopulation(R, t, groupName, eventName, isDiff, plotFileName)


nUnits = size(R, 1);

figure_tr_inch(8, 8); clf;
subaxis(1, 1, 1, 'MB', 0.1, 'ML', 0.1, 'MR', 0.07);
hold on;
imagesc(t, 1:nUnits, R);
plot([0 0], [0 nUnits + 1], '-', 'Color', 0.3 * ones(3, 1));
hcb = colorbar;

xlim(t([1 end]));
ylim([0 nUnits + 1]);

xlabel(sprintf('Time from %s (s)', eventName));
ylabel('Unit Number');

titleText = sprintf('%s: %s (%d Units)', groupName, eventName, nUnits);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(titleText, 'FontSize', 14, titleParams{:});

set(gca, 'box', 'off');
set(gca, 'FontSize', 14);

if isDiff
    clabel = 'InRF - ExRF Norm. Firing Rate';
else
    clabel = 'InRF Norm. Firing Rate';
end
ylabel(hcb, clabel, 'Rotation', 270, 'FontSize', 14, ...
        'Units', 'Normalized', 'Position', [3.4 0.5 0]); % colorbar label

%% save
if ~isempty(plotFileName)
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
