function plotLatencyDiff(arrayHoldBalLatencyInRF, arrayHoldBalLatencyExRF, goodUnits, isDPul, isVPul)

maxLatency = 0.125;
minLatency = 0.025;
goodUnits = ~isnan(arrayHoldBalLatencyInRF) & ~isnan(arrayHoldBalLatencyExRF) & ...
        arrayHoldBalLatencyInRF <= maxLatency & arrayHoldBalLatencyExRF <= maxLatency & ...
        (arrayHoldBalLatencyInRF >= minLatency | arrayHoldBalLatencyExRF >= minLatency) & ...
        goodUnits;
arrayOnsetHoldLatencyDiff = arrayHoldBalLatencyInRF(goodUnits) - arrayHoldBalLatencyExRF(goodUnits);
fprintf('Mean array onset latency Attend-RF: %0.3f s (median: %0.3f s)\n', mean(arrayHoldBalLatencyInRF(goodUnits)), median(arrayHoldBalLatencyInRF(goodUnits)));
fprintf('Mean array onset latency Attend-Away: %0.3f s (median: %0.3f s)\n', mean(arrayHoldBalLatencyExRF(goodUnits)), median(arrayHoldBalLatencyExRF(goodUnits)));
fprintf('Number of units with latency reduction >= 10 ms: %d (%d%%)\n', ...
        sum(arrayOnsetHoldLatencyDiff <= -0.01), round(sum(arrayOnsetHoldLatencyDiff <= -0.01)/numel(arrayOnsetHoldLatencyDiff) * 100));
fprintf('Number of units with latency increase >= 10 ms: %d (%d%%)\n', ...
        sum(arrayOnsetHoldLatencyDiff >= 0.01), round(sum(arrayOnsetHoldLatencyDiff >= 0.01)/numel(arrayOnsetHoldLatencyDiff) * 100));

goodUnitsDPul = isDPul & goodUnits;
goodUnitsVPul = isVPul & goodUnits;

%%
latBounds = [0 0.13];
maxAbsDiffFR = max(abs(arrayOnsetHoldLatencyDiff));

binStep = 0.005;
histXBounds = [-ceil(maxAbsDiffFR / binStep) ceil(maxAbsDiffFR / binStep)] * binStep;
histBinEdges = histXBounds(1):binStep:histXBounds(2);

cols = lines(6);
colDPul = cols(3,:);
colVPul = cols(5,:);

%%
figure_tr_inch(16, 5);
subaxis(1, 3, 1);
hold on;
% plot(arrayHoldBalLatencyInRF(goodUnits), arrayHoldBalLatencyExRF(goodUnits), '.', 'MarkerSize', 20);
h1 = plot(arrayHoldBalLatencyInRF(goodUnitsDPul), arrayHoldBalLatencyExRF(goodUnitsDPul), '.', 'MarkerSize', 20, 'Color', colDPul);
h2 = plot(arrayHoldBalLatencyInRF(goodUnitsVPul), arrayHoldBalLatencyExRF(goodUnitsVPul), '.', 'MarkerSize', 20, 'Color', colVPul);
plot([0 1], [0 1], 'Color', 0.3*ones(3, 1)); 
axis equal;
xlim(latBounds); 
ylim(latBounds);
xlabel('Array Onset Latency Attend-RF (s)');
ylabel('Array Onset Latency Attend-Away (s)');
box off;
% legend([h1 h2], {'dPul', 'vPul'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

subaxis(1, 3, 2); 
hold on;
histH = histogram(arrayOnsetHoldLatencyDiff, histBinEdges);
histH.FaceColor = cols(4,:);
origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
xlabel('Array Onset Latency Difference (s)');
ylabel('Number of Units');
box off;
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);

subaxis(1, 3, 3);
hold on;
corrP3DPul = cdfplot(arrayHoldBalLatencyInRF(goodUnits));
corrP1DPul = cdfplot(arrayHoldBalLatencyExRF(goodUnits));
set(corrP3DPul, 'LineWidth', 3);
set(corrP1DPul, 'LineWidth', 3);
xlim(latBounds); 
xlabel('Array Onset Latency (s)');
ylabel('Proportion of Units');
legend({'Attend-RF', 'Attend-Away'}, 'Location', 'SouthEast');
set(gca, 'FontSize', 14);
set(gca, 'LineWidth', 1);


meanArrayOnsetHoldLatencyDiff = mean(arrayOnsetHoldLatencyDiff);
medianArrayOnsetHoldLatencyDiff = median(arrayOnsetHoldLatencyDiff);
p = signrank(arrayOnsetHoldLatencyDiff);
fprintf('mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d\n', ...
        meanArrayOnsetHoldLatencyDiff, medianArrayOnsetHoldLatencyDiff, p, sum(goodUnits));