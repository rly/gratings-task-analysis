function plotLfpPower(relPower1, relPower2, fAxis, xBounds, yBounds, col1, col2, label1, label2)

nChannels1 = size(relPower1, 1);
nChannels2 = size(relPower2, 1);

meanPower1 = mean(relPower1);
meanPower2 = mean(relPower2);
sePower1 = std(relPower1) / sqrt(nChannels1);
sePower2 = std(relPower2) / sqrt(nChannels2);

figure_tr_inch(7, 5); 
subaxis(1, 1, 1, 'MB', 0.14, 'ML', 0.15);
hold on;
fillH = jbfill(fAxis, meanPower2 - sePower2, meanPower2 + sePower2, ...
        col2, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanPower1 - sePower1, meanPower1 + sePower1, ...
        col1, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
h1 = plot(fAxis, meanPower1, 'LineWidth', 2, 'Color', col1);
h2 = plot(fAxis, meanPower2, 'LineWidth', 2, 'Color', col2);
legend([h1 h2], {sprintf(' %s (N=%d)', label1, nChannels1), sprintf(' %s (N=%d)', label2, nChannels2)}, ...
        'box', 'off', 'Location', 'SouthEast');
xlim(xBounds);
ylim(yBounds);
xlabel('Frequency (Hz)');
ylabel('Percent Change in Power Rel. to Baseline');
set(gca, 'box', 'off');
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);