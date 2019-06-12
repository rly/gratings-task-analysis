function grcPlotLfpPower2(powerLF, powerHF, fAxisLF, fAxisHF, conds, varargin)
% powerLF, powerHF - channel x frequency
% conds - logical, channel x condition
xBounds = [fAxisLF([1 end]); fAxisHF([1 end])];
yBounds = [-1 1];
cols = lines(size(conds, 2));
lineLabels = cell(size(conds, 2), 1);
ylabelText = 'Power (dB/Hz)';
titleText = 'Power';
doDB = 1;
overridedefaults(who, varargin);

%%
nCond = size(conds, 2); % each column is sep condition of row indices into powerLF

for i = 1:nCond
    fprintf('%s: %d channels\n', lineLabels{i}, sum(conds(:,i)));
end

meanPowerLF = nan(nCond, numel(fAxisLF));
sePowerLF = nan(nCond, numel(fAxisLF));
seLBPowerLF = nan(nCond, numel(fAxisLF));
seUBPowerLF = nan(nCond, numel(fAxisLF));
meanPowerHF = nan(nCond, numel(fAxisHF));
sePowerHF = nan(nCond, numel(fAxisHF));
seLBPowerHF = nan(nCond, numel(fAxisHF));
seUBPowerHF = nan(nCond, numel(fAxisHF));
for i = 1:nCond
    powerCond = powerLF(conds(:,i),:);
    meanPowerLF(i,:) = mean(powerCond, 1);
    sePowerLF(i,:) = std(powerCond) / sqrt(sum(conds(:,i)));
    seLBPowerLF(i,:) = meanPowerLF(i,:) - sePowerLF(i,:);
    seUBPowerLF(i,:) = meanPowerLF(i,:) + sePowerLF(i,:);
    
    powerCond = powerHF(conds(:,i),:);
    meanPowerHF(i,:) = mean(powerCond, 1);
    sePowerHF(i,:) = std(powerCond) / sqrt(sum(conds(:,i)));
    seLBPowerHF(i,:) = meanPowerHF(i,:) - sePowerHF(i,:);
    seUBPowerHF(i,:) = meanPowerHF(i,:) + sePowerHF(i,:);
end

if doDB
    meanPowerLF = pow2db(meanPowerLF);
    meanPowerHF = pow2db(meanPowerHF);
    seLBPowerLF = pow2db(seLBPowerLF);
    seUBPowerLF = pow2db(seUBPowerLF);
    seLBPowerHF = pow2db(seLBPowerHF);
    seUBPowerHF = pow2db(seUBPowerHF);
end
clear powerCond;

%% plot positionining params
% these params are optimized for dB units
plotLFW = diff(xBounds(1,:)) / 102;
plotHFW = diff(xBounds(2,:)) / 102;
plotH = 0.8;

plotLFLeft = 0.14;
plotHFLeft = plotLFLeft + plotLFW + 0.12;

btm = 0.16;

figure_tr_inch(6.75, 4.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.14 0.92 0.80], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(titleText, 'FontSize', 18, titleParams{:});

xlabel('Frequency (Hz)', 'Visible', 'on', 'FontSize', 18);

%% LF
axes('Position', [plotLFLeft btm plotLFW plotH]); 
hold on;
for i = nCond:-1:1
    fillH = jbfill(fAxisLF, seLBPowerLF(i,:), seUBPowerLF(i,:), cols(i,:), ones(3, 1), 0.5);
    uistack(fillH, 'bottom');
end
hold on;
plot([-100 100], [0 0], 'k');

plotHs = nan(nCond, 1);
legendLabels = cell(nCond, 1);
for i = 1:nCond
    plotHs(i) = plot(fAxisLF, meanPowerLF(i,:), 'LineWidth', 4, 'Color', cols(i,:));
    legendLabels{i} = sprintf(' %s (N=%d)', lineLabels{i}, sum(conds(:,i)));
end
% legend(plotHs, legendLabels, 'box', 'off', 'Location', 'SouthEast');
xlim(xBounds(1,:));
ylim(yBounds(1,:));
% xlabel('Frequency (Hz)');
ylabel(ylabelText);
if diff(yBounds(1,:)) <= 0.04
    set(gca, 'YTick', 0:0.01:1);
end
set(gca, 'box', 'off');
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
grid on;

%% HF
axes('Position', [plotHFLeft btm plotHFW plotH]); 
hold on;
for i = nCond:-1:1
    fillH = jbfill(fAxisHF, seLBPowerHF(i,:), seUBPowerHF(i,:), cols(i,:), ones(3, 1), 0.5);
    uistack(fillH, 'bottom');
end
hold on;
plot([-100 100], [0 0], 'k');

plotHs = nan(nCond, 1);
legendLabels = cell(nCond, 1);
for i = 1:nCond
    plotHs(i) = plot(fAxisHF, meanPowerHF(i,:), 'LineWidth', 4, 'Color', cols(i,:));
    legendLabels{i} = sprintf(' %s', lineLabels{i});
end
% legend(plotHs, legendLabels, 'box', 'off', 'Location', 'NorthEast');
xlim(xBounds(2,:));
ylim(yBounds(2,:));
% xlabel('Frequency (Hz)');
% ylabel('Power');
% set(gca, 'YTickLabel', []);
set(gca, 'box', 'off');
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
grid on;
