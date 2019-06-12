cueWindowOffset = [-0.1 0.5];
cueWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, cueWindowOffset);
arrayWindowOffset = [-0.1 0.4];
arrayWindowIndices = getTimeLogicalWithTolerance(arrayOnsetT, arrayWindowOffset);
targetDimWindowOffset = [-0.2 0.2];
targetDimWindowIndices = getTimeLogicalWithTolerance(targetDimT, targetDimWindowOffset);
periSaccadeWindowOffset = [-0.2 0.2];
periSaccadeWindowIndices = getTimeLogicalWithTolerance(exitFixationT, periSaccadeWindowOffset);

cueOnsetSpdfDPul = spdfInfo.cueOnsetSpdfInRFNorm(isInDPulvinar & isSignificantCueResponseInc,cueWindowIndices);
meanCueOnsetSpdfDPul = mean(cueOnsetSpdfDPul);
seCueOnsetSpdfDPul = std(cueOnsetSpdfDPul) / sqrt(size(cueOnsetSpdfDPul, 1));
cueOnsetSpdfVPul = spdfInfo.cueOnsetSpdfInRFNorm(isInVPulvinar & isSignificantCueResponseInc,cueWindowIndices);
meanCueOnsetSpdfVPul = mean(cueOnsetSpdfVPul);
seCueOnsetSpdfVPul = std(cueOnsetSpdfVPul) / sqrt(size(cueOnsetSpdfVPul, 1));

arrayOnsetHoldSpdfDPul = spdfInfo.arrayOnsetHoldSpdfInRFNorm(isInDPulvinar & isSignificantCueResponseInc,arrayWindowIndices);
meanArrayOnsetHoldSpdfDPul = mean(arrayOnsetHoldSpdfDPul);
seArrayOnsetHoldSpdfDPul = std(arrayOnsetHoldSpdfDPul) / sqrt(size(arrayOnsetHoldSpdfDPul, 1));
arrayOnsetHoldSpdfVPul = spdfInfo.arrayOnsetHoldSpdfInRFNorm(isInVPulvinar & isSignificantCueResponseInc,arrayWindowIndices);
meanArrayOnsetHoldSpdfVPul = mean(arrayOnsetHoldSpdfVPul);
seArrayOnsetHoldSpdfVPul = std(arrayOnsetHoldSpdfVPul) / sqrt(size(arrayOnsetHoldSpdfVPul, 1));

targetDimSpdfDPul = spdfInfo.targetDimSpdfInRFNorm(isInDPulvinar & isSignificantCueResponseInc,targetDimWindowIndices);
meanTargetDimSpdfDPul = mean(targetDimSpdfDPul);
seTargetDimSpdfDPul = std(targetDimSpdfDPul) / sqrt(size(targetDimSpdfDPul, 1));
targetDimSpdfVPul = spdfInfo.targetDimSpdfInRFNorm(isInVPulvinar & isSignificantCueResponseInc,targetDimWindowIndices);
meanTargetDimSpdfVPul = mean(targetDimSpdfVPul);
seTargetDimSpdfVPul = std(targetDimSpdfVPul) / sqrt(size(targetDimSpdfVPul, 1));

exitFixationSpdfDPul = spdfInfo.exitFixationSpdfInRFNorm(isInDPulvinar & isSignificantCueResponseInc,periSaccadeWindowIndices);
meanExitFixationSpdfDPul = mean(exitFixationSpdfDPul);
seExitFixationSpdfDPul = std(exitFixationSpdfDPul) / sqrt(size(exitFixationSpdfDPul, 1));
exitFixationSpdfVPul = spdfInfo.exitFixationSpdfInRFNorm(isInVPulvinar & isSignificantCueResponseInc,periSaccadeWindowIndices);
meanExitFixationSpdfVPul = mean(exitFixationSpdfVPul);
seExitFixationSpdfVPul = std(exitFixationSpdfVPul) / sqrt(size(exitFixationSpdfVPul, 1));

size(cueOnsetSpdfDPul, 1)
size(cueOnsetSpdfVPul, 1)

isShowLabels = 1;
yBounds = [-0.13 0.31];

outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\SUA_MUA_GRATINGS_SUMMARY\';
plotFileName = sprintf('%s/allSessions-meanSpdfsCueArrayExitFix-dpul-vpul-sfn-v%d.png', outputDir, v);

%%
f = figure_tr_inch(13, 4); 
clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% location params
spdfW = 0.22;
spdfH = 0.776;
spdfCueOnsetW = 0.35;

spdfCueOnsetLeft = 0.08;
% spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.04;
spdfTargetDimLeft = spdfCueOnsetLeft + spdfCueOnsetW + 0.05;
spdfExitFixationLeft = spdfTargetDimLeft + spdfW + 0.05;

btm = 0.2;

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfCueOnsetW spdfH]); 

hold on;
plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min enter fixation/lever press
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(cueOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(cueOnsetT(cueWindowIndices), ...
        meanCueOnsetSpdfVPul + seCueOnsetSpdfVPul, ...
        meanCueOnsetSpdfVPul - seCueOnsetSpdfVPul, ...
        vPulCol, vPulCol, 0.5);
plot(cueOnsetT(cueWindowIndices), meanCueOnsetSpdfVPul, 'Color', vPulCol, 'LineWidth', 4);
jbfill(cueOnsetT(cueWindowIndices), ...
        meanCueOnsetSpdfDPul + seCueOnsetSpdfDPul, ...
        meanCueOnsetSpdfDPul - seCueOnsetSpdfDPul, ...
        dPulCol, dPulCol, 0.5);
plot(cueOnsetT(cueWindowIndices), meanCueOnsetSpdfDPul, 'Color', dPulCol, 'LineWidth', 4);

xlim(cueWindowOffset);
if isShowLabels
    xlabel('Time from Cue Onset (s)');
    ylabel('Normalized Firing Rate');
end
set(gca, 'FontSize', 18);
set(gca, 'XTick', -0.3:0.1:0.5);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');

%% spdf for array onset
% axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 
% 
% hold on;
% 
% plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
% plot(arrayOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
% 
% jbfill(arrayOnsetT(arrayWindowIndices), ...
%         meanArrayOnsetHoldSpdfVPul + seArrayOnsetHoldSpdfVPul, ...
%         meanArrayOnsetHoldSpdfVPul - seArrayOnsetHoldSpdfVPul, ...
%         vPulCol, vPulCol, 0.5);
% plot(arrayOnsetT(arrayWindowIndices), meanArrayOnsetHoldSpdfVPul, 'Color', vPulCol, 'LineWidth', 4);
% 
% jbfill(arrayOnsetT(arrayWindowIndices), ...
%         meanArrayOnsetHoldSpdfDPul + seArrayOnsetHoldSpdfDPul, ...
%         meanArrayOnsetHoldSpdfDPul - seArrayOnsetHoldSpdfDPul, ...
%         dPulCol, dPulCol, 0.5);
% plot(arrayOnsetT(arrayWindowIndices), meanArrayOnsetHoldSpdfDPul, 'Color', dPulCol, 'LineWidth', 4);
% 
% xlim(arrayWindowOffset);
% if isShowLabels
%     xlabel('Time from Array Onset (s)');
% end
% set(gca, 'FontSize', 18);
% set(gca, 'XTick', -0.3:0.1:0.3);
% set(gca, 'YTickLabel', []);
% set(gca, 'LineWidth', 1);
% set(gca, 'TickDir', 'out');

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

hold on;

plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(targetDimT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(targetDimT(targetDimWindowIndices), ...
        meanTargetDimSpdfVPul + seTargetDimSpdfVPul, ...
        meanTargetDimSpdfVPul - seTargetDimSpdfVPul, ...
        vPulCol, vPulCol, 0.5);
plot(targetDimT(targetDimWindowIndices), meanTargetDimSpdfVPul, 'Color', vPulCol, 'LineWidth', 4);

jbfill(targetDimT(targetDimWindowIndices), ...
        meanTargetDimSpdfDPul + seTargetDimSpdfDPul, ...
        meanTargetDimSpdfDPul - seTargetDimSpdfDPul, ...
        dPulCol, dPulCol, 0.5);
plot(targetDimT(targetDimWindowIndices), meanTargetDimSpdfDPul, 'Color', dPulCol, 'LineWidth', 4);

xlim(targetDimWindowOffset);
if isShowLabels
    xlabel('Time from Target Dim (s)');
end
set(gca, 'FontSize', 18);
set(gca, 'XTick', -0.3:0.1:0.3);
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');

%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

hold on;

% plot([-0.28 -0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min rt
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot(exitFixationT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

jbfill(exitFixationT(periSaccadeWindowIndices), ...
        meanExitFixationSpdfVPul + seExitFixationSpdfVPul, ...
        meanExitFixationSpdfVPul - seExitFixationSpdfVPul, ...
        vPulCol, vPulCol, 0.5);
plot(exitFixationT(periSaccadeWindowIndices), meanExitFixationSpdfVPul, 'Color', vPulCol, 'LineWidth', 4);

jbfill(exitFixationT(periSaccadeWindowIndices), ...
        meanExitFixationSpdfDPul + seExitFixationSpdfDPul, ...
        meanExitFixationSpdfDPul - seExitFixationSpdfDPul, ...
        dPulCol, dPulCol, 0.5);
plot(exitFixationT(periSaccadeWindowIndices), meanExitFixationSpdfDPul, 'Color', dPulCol, 'LineWidth', 4);

xlim(periSaccadeWindowOffset);
if isShowLabels
    xlabel('Time from End Fixation (s)');
end
set(gca, 'FontSize', 18);
set(gca, 'XTick', -0.3:0.1:0.3);
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');

%% adjust all ylim
%     allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf) ylim(axExitFixationSpdf)];
ylim(axCueOnsetSpdf, yBounds);
% ylim(axArrayOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);
ylim(axExitFixationSpdf, yBounds);

%% save
if ~isempty(plotFileName)
    export_fig(plotFileName, '-nocrop');
end