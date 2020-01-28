function quickSpdfAllEvents4InARowPopMean(...
        cueOnsetSpdfInRFNorm, cueOnsetSpdfExRFNorm, ...
        arrayOnsetHoldSpdfInRFNorm, arrayOnsetHoldSpdfExRFNorm, targetDimSpdfInRFNorm, ...
        targetDimSpdfExRFNorm, exitFixationSpdfInRFNorm, exitFixationSpdfExRFNorm, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, yBounds, isShowLabels, plotFileName)

% check size assert
nSessions = size(cueOnsetSpdfInRFNorm, 1);
assert(all([size(cueOnsetSpdfExRFNorm, 1);
        size(arrayOnsetHoldSpdfInRFNorm, 1);
        size(arrayOnsetHoldSpdfExRFNorm, 1);
        size(targetDimSpdfInRFNorm, 1);
        size(targetDimSpdfExRFNorm, 1); 
        size(exitFixationSpdfInRFNorm, 1);
        size(exitFixationSpdfExRFNorm, 1)] == nSessions))

meanCueOnsetSpdfInRFNorm = mean(cueOnsetSpdfInRFNorm);
meanCueOnsetSpdfExRFNorm = mean(cueOnsetSpdfExRFNorm);
meanArrayOnsetHoldSpdfInRFNorm = mean(arrayOnsetHoldSpdfInRFNorm);
meanArrayOnsetHoldSpdfExRFNorm = mean(arrayOnsetHoldSpdfExRFNorm);
meanTargetDimSpdfInRFNorm = mean(targetDimSpdfInRFNorm);
meanTargetDimSpdfExRFNorm = mean(targetDimSpdfExRFNorm);
meanExitFixationSpdfInRFNorm = mean(exitFixationSpdfInRFNorm);
meanExitFixationSpdfExRFNorm = mean(exitFixationSpdfExRFNorm);

seCueOnsetSpdfInRFNorm = std(cueOnsetSpdfInRFNorm)/sqrt(nSessions);
seCueOnsetSpdfExRFNorm = std(cueOnsetSpdfExRFNorm)/sqrt(nSessions);
seArrayOnsetHoldSpdfInRFNorm = std(arrayOnsetHoldSpdfInRFNorm)/sqrt(nSessions);
seArrayOnsetHoldSpdfExRFNorm = std(arrayOnsetHoldSpdfExRFNorm)/sqrt(nSessions);
seTargetDimSpdfInRFNorm = std(targetDimSpdfInRFNorm)/sqrt(nSessions);
seTargetDimSpdfExRFNorm = std(targetDimSpdfExRFNorm)/sqrt(nSessions);
seExitFixationSpdfInRFNorm = std(exitFixationSpdfInRFNorm)/sqrt(nSessions);
seExitFixationSpdfExRFNorm = std(exitFixationSpdfExRFNorm)/sqrt(nSessions);

inRFCol = [0.9 0 0];
exRFCol = [0 0 0.9];

%%
f = figure_tr_inch(16, 3.75); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% location params
spdfW = 0.22;
spdfH = 0.74;

spdfCueOnsetLeft = 0.05;
spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.015;
spdfTargetDimLeft = spdfArrayOnsetLeft + spdfW + 0.015;
spdfExitFixationLeft = spdfTargetDimLeft + spdfW + 0.015;

btm = 0.18;

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
hold on;
plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min enter fixation/lever press
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(cueOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(cueOnsetT, ...
        meanCueOnsetSpdfExRFNorm + seCueOnsetSpdfExRFNorm, ...
        meanCueOnsetSpdfExRFNorm - seCueOnsetSpdfExRFNorm, ...
        exRFCol, exRFCol, 0.5);
plot(cueOnsetT, meanCueOnsetSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);
jbfill(cueOnsetT, ...
        meanCueOnsetSpdfInRFNorm + seCueOnsetSpdfInRFNorm, ...
        meanCueOnsetSpdfInRFNorm - seCueOnsetSpdfInRFNorm, ...
        inRFCol, inRFCol, 0.5);
plot(cueOnsetT, meanCueOnsetSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

xlim(xBounds);
if isShowLabels
	xlabel('Time from Cue Onset (s)');
    ylabel('Normalized Firing Rate');
end
set(gca, 'FontSize', 16);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'FontName', 'Arial');
set(gca, 'LineWidth', 2);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
hold on;

plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(arrayOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(arrayOnsetT, ...
        meanArrayOnsetHoldSpdfExRFNorm + seArrayOnsetHoldSpdfExRFNorm, ...
        meanArrayOnsetHoldSpdfExRFNorm - seArrayOnsetHoldSpdfExRFNorm, ...
        exRFCol, exRFCol, 0.5);
plot(arrayOnsetT, meanArrayOnsetHoldSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);

jbfill(arrayOnsetT, ...
        meanArrayOnsetHoldSpdfInRFNorm + seArrayOnsetHoldSpdfInRFNorm, ...
        meanArrayOnsetHoldSpdfInRFNorm - seArrayOnsetHoldSpdfInRFNorm, ...
        inRFCol, inRFCol, 0.5);
plot(arrayOnsetT, meanArrayOnsetHoldSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

xlim(xBounds);
if isShowLabels
    xlabel('Time from Array Onset (s)');
end
set(gca, 'FontSize', 16);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'FontName', 'Arial');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
hold on;

plot([0.28 0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min rt
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(targetDimT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(targetDimT, ...
        meanTargetDimSpdfExRFNorm + seTargetDimSpdfExRFNorm, ...
        meanTargetDimSpdfExRFNorm - seTargetDimSpdfExRFNorm, ...
        exRFCol, exRFCol, 0.5);
plot(targetDimT, meanTargetDimSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);

jbfill(targetDimT, ...
        meanTargetDimSpdfInRFNorm + seTargetDimSpdfInRFNorm, ...
        meanTargetDimSpdfInRFNorm - seTargetDimSpdfInRFNorm, ...
        inRFCol, inRFCol, 0.5);
plot(targetDimT, meanTargetDimSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

xlim(xBounds);
if isShowLabels
    xlabel('Time from Target Dimming (s)');
end
set(gca, 'FontSize', 16);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'FontName', 'Arial');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);

%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
hold on;

plot([-0.28 -0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min rt
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(exitFixationT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(exitFixationT, ...
        meanExitFixationSpdfExRFNorm + seExitFixationSpdfExRFNorm, ...
        meanExitFixationSpdfExRFNorm - seExitFixationSpdfExRFNorm, ...
        exRFCol, exRFCol, 0.5);
plot(exitFixationT, meanExitFixationSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);

jbfill(exitFixationT, ...
        meanExitFixationSpdfInRFNorm + seExitFixationSpdfInRFNorm, ...
        meanExitFixationSpdfInRFNorm - seExitFixationSpdfInRFNorm, ...
        inRFCol, inRFCol, 0.5);
plot(exitFixationT, meanExitFixationSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

xlim(xBounds);
if isShowLabels
    xlabel('Time from Break Fixation (s)');
end
set(gca, 'FontSize', 16);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'FontName', 'Arial');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);

%% adjust all ylim
allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf) ylim(axTargetDimSpdf) ylim(axExitFixationSpdf)];
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);
ylim(axExitFixationSpdf, yBounds);

%% save
if ~isempty(plotFileName)
    export_fig(plotFileName, '-nocrop');
end