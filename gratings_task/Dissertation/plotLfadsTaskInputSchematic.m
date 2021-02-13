
x1 = [zeros(100, 1); ones(100, 1); zeros(550, 1); zeros(750, 1); zeros(300, 1); zeros(200, 1)];
x2 = [zeros(100, 1); zeros(100, 1); zeros(550, 1); zeros(750, 1); zeros(300, 1); zeros(200, 1)];
x3 = [zeros(100, 1); zeros(100, 1); zeros(550, 1); ones(750, 1); ones(300, 1); zeros(200, 1)];
x4 = [zeros(100, 1); zeros(100, 1); zeros(550, 1); -1*ones(750, 1); -1*ones(300, 1); zeros(200, 1)];
x5 = [zeros(100, 1); zeros(100, 1); zeros(550, 1); ones(750, 1); 1/2*ones(300, 1); zeros(200, 1)];
x6 = [zeros(100, 1); zeros(100, 1); zeros(550, 1); -1*ones(750, 1); -1*ones(300, 1); zeros(200, 1)];

xTicks = [100 200 750 1500 1800];
xBounds = [0 numel(x1)];
yBounds = [-1 1];

figure_tr_inch(7, 5);

subaxis(6, 1, 1, 'MT', 0.05, 'ML', 0.30, 'MB', 0.13);
plot(x1, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Cue at Loc. 1:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);

subaxis(6, 1, 2);
plot(x2, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Cue at Loc. 3:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);

subaxis(6, 1, 3);
plot(x3, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Shape at Loc. 1:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);


subaxis(6, 1, 4);
plot(x4, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Shape at Loc. 2:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);

subaxis(6, 1, 5);
plot(x5, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Shape at Loc. 3:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);

subaxis(6, 1, 6);
plot(x6, 'LineWidth', 2);
box off;
xlim(xBounds);
ylim(yBounds);
set(gca, 'LineWidth', 1);
set(gca, 'FontSize', 14);
set(gca, 'XTick', xTicks);
set(gca, 'XTickLabel', '');
set(gca, 'TickDir', 'out');
text(-850, 0, {'Shape at Loc. 4:'}, 'HorizontalAlignment', 'Left', 'FontSize', 16);

text(100+25, -2.5, {'Cue', 'Onset'}, 'HorizontalAlignment', 'Right', 'FontSize', 14);
text(200-25, -2.5, {'Cue', 'Offset'}, 'HorizontalAlignment', 'Left', 'FontSize', 14);
text(750, -2.5, {'Array', 'Onset'}, 'HorizontalAlignment', 'Center', 'FontSize', 14);
text(1500+50, -2.5, {'Target', 'Dim'}, 'HorizontalAlignment', 'Right', 'FontSize', 14);
text(1800, -2.5, {'Response', ''}, 'HorizontalAlignment', 'Center', 'FontSize', 14);

processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
outputDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_LFADS/Dissertation');
saveFileName = sprintf('%s/lfadsTaskInputSchematic.png', outputDir);
export_fig(saveFileName, '-nocrop');