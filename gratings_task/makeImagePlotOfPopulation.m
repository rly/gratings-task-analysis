function [hImage, hColorBar] = makeImagePlotOfPopulation(R, t, eventName, isDiff)

nUnits = size(R, 1);

hold on;
hImage = imagesc(t, 1:nUnits, R);
plot([0 0], [0 nUnits + 1], '-', 'Color', 0.3 * ones(3, 1));
hColorBar = colorbar;

xlim(t([1 end]));
ylim([0.5 nUnits + 0.5]);
caxis([-1 1]);

xlabel(sprintf('Time from %s (s)', eventName));

titleText = eventName;
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(titleText, 'FontSize', 14, titleParams{:});

set(gca, 'box', 'off');
set(gca, 'FontSize', 16);
set(gca, 'FontName', 'Calibri');
set(gca, 'LineWidth', 2);

if isDiff
    clabel = 'InRF - ExRF Norm. Firing Rate';
else
    clabel = 'InRF Norm. Firing Rate';
end
ylabel(hColorBar, clabel, 'Rotation', 270, 'FontSize', 14, ...
        'Units', 'Normalized', 'Position', [3.4 0.5 0]); % colorbar label
    
set(hColorBar, 'Visible', 'off');
