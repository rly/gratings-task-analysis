function plotLfpPower4(power1a, power1b, power2a, power2b, fAxis, xBounds, yBounds, col1, col2, ...
        label1a, label1b, label2a, label2b)

nChannels1 = size(power1a, 1);
nChannels2 = size(power2a, 1);

meanPower1a = mean(power1a);
meanPower1b = mean(power1b);
meanPower2a = mean(power2a);
meanPower2b = mean(power2b);
sePower1a = std(power1a) / sqrt(nChannels1);
sePower1b = std(power1b) / sqrt(nChannels1);
sePower2a = std(power2a) / sqrt(nChannels2);
sePower2b = std(power2b) / sqrt(nChannels2);

figure_tr_inch(7, 5); 
subaxis(1, 1, 1, 'MB', 0.14, 'ML', 0.15);
hold on;
fillH = jbfill(fAxis, meanPower2b - sePower2b, meanPower2b + sePower2b, ...
        col2, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanPower2a - sePower2a, meanPower2a + sePower2a, ...
        col2, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanPower1b - sePower1b, meanPower1b + sePower1b, ...
        col1, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
fillH = jbfill(fAxis, meanPower1a - sePower1a, meanPower1a + sePower1a, ...
        col1, ones(3, 1), 0.3);
uistack(fillH, 'bottom');
hold on;
plot([-100 100], [0 0], 'k');
h1a = plot(fAxis, meanPower1a, '--', 'LineWidth', 2, 'Color', col1);
h1b = plot(fAxis, meanPower1b, 'LineWidth', 2, 'Color', col1);
h2a = plot(fAxis, meanPower2a, '--', 'LineWidth', 2, 'Color', col2);
h2b = plot(fAxis, meanPower2b, 'LineWidth', 2, 'Color', col2);
legend([h1a h1b h2a h2b], {sprintf(' %s (N=%d)', label1a, nChannels1), ...
        sprintf(' %s (N=%d)', label1b, nChannels1), ...
        sprintf(' %s (N=%d)', label2a, nChannels2), ...
        sprintf(' %s (N=%d)', label2b, nChannels2)}, ...
        'box', 'off', 'Location', 'SouthEast');
xlim(xBounds);
ylim(yBounds);
xlabel('Frequency (Hz)');
ylabel('Power');
set(gca, 'box', 'off');
set(gca, 'FontSize', 16);
set(gca, 'LineWidth', 2);