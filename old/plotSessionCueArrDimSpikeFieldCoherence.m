function plotSessionCueArrDimSpikeFieldCoherence(sessionName, name, ...
        cueInRFC, cueExRFC, arrInRFC, arrExRFC, dimInRFC, dimExRFC, ...
        tCue, tArr, tDim, f, conditionsInfo, ...
        isZeroDistractors, saveFile, varargin)

isVisible = 0;
posCLim = [0 0.2];
diffCLim = [-0.1 0.1];
overridedefaults(who, varargin);

periCueWindow = conditionsInfo.periCueOnsetWindow;
periArrWindow = conditionsInfo.periArrayOnsetWindow;
periDimWindow = conditionsInfo.periTargetDimWindow;

load('params.mat');

figure_tr_inch(16,8); clf;
set(gcf, 'Color', 'white');
% set(gcf, 'InvertHardCopy','off');
set(gcf, 'Renderer', 'painters');
if ~isVisible
    set(gcf, 'Visible', 'off');
end

%% location params

cohgramW = 0.3; % width of a single coherogram plot
cohgramH = 0.23; % height of a single coherogram plot

diffCohgramBtm = 0.07; % far bottom
exRFCohgramBtm = diffCohgramBtm + cohgramH + 0.07;
inRFCohgramBtm = exRFCohgramBtm + cohgramH + 0.07;

col1Left = 0.05; % far left
col2Left = col1Left + cohgramW + 0.03;
col3Left = col2Left + cohgramW + 0.03;

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig,'Title'), 'Visible', 'on')

% format: Session L110524, LIP Neuron B (140 trials)
zeroDAppend = '';
if isZeroDistractors
    zeroDAppend = ' (0 Distractors)'; 
end

modTitle = sprintf('Session %s%s - %s - Spike Field Coherence Around Cue Onset, Array Onset, and Target Dimming', ...
        sessionName, zeroDAppend, name);
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 15, titleParams{:});

%% inRF
axInRFCohgramBig = axes('Position', [col1Left inRFCohgramBtm cohgramW cohgramH]);
set(get(axInRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tCue - periCueWindow(1), f, cueInRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
ylabel('Frequency (Hz)');

title(sprintf('InRF P%d Cue SFC (%d trials)', conditionsInfo.inRFLoc, conditionsInfo.numTrials.inRF));

%% exRF
axExRFCohgramBig = axes('Position', [col1Left exRFCohgramBtm cohgramW cohgramH]);
set(get(axExRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tCue - periCueWindow(1), f, cueExRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
ylabel('Frequency (Hz)');

title(sprintf('ExRF P%d Cue SFC (%d trials)', conditionsInfo.exRFLoc, conditionsInfo.numTrials.exRF));

%% diff
axDiffCohgramBig = axes('Position', [col1Left diffCohgramBtm cohgramW cohgramH]);
set(get(axDiffCohgramBig, 'Title'), 'Visible', 'on')

cueDiffC = cueInRFC - cueExRFC;

imagesc(tCue - periCueWindow(1), f, cueDiffC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(diffCLim);
colormap(gca, getCoolWarmMap());
colorbar;
xlabel('Time from event onset (s)');
ylabel('Frequency (Hz)');

title('InRF - ExRF Cue SFC');

%% inRF
axInRFCohgramBig = axes('Position', [col2Left inRFCohgramBtm cohgramW cohgramH]);
set(get(axInRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tArr - periArrWindow(1), f, arrInRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
% ylabel('Frequency (Hz)');

title(sprintf('InRF P%d Array SFC (%d trials)', conditionsInfo.inRFLoc, conditionsInfo.numTrials.inRF));

%% exRF
axExRFCohgramBig = axes('Position', [col2Left exRFCohgramBtm cohgramW cohgramH]);
set(get(axExRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tArr - periArrWindow(1), f, arrExRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
% ylabel('Frequency (Hz)');

title(sprintf('ExRF P%d Array SFC (%d trials)', conditionsInfo.exRFLoc, conditionsInfo.numTrials.exRF));

%% diff
axDiffCohgramBig = axes('Position', [col2Left diffCohgramBtm cohgramW cohgramH]);
set(get(axDiffCohgramBig, 'Title'), 'Visible', 'on')

arrDiffC = arrInRFC - arrExRFC;

imagesc(tArr - periArrWindow(1), f, arrDiffC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(diffCLim);
colormap(gca, getCoolWarmMap());
colorbar;
xlabel('Time from event onset (s)');
% ylabel('Frequency (Hz)');

title('InRF - ExRF Array SFC');

%% inRF
axInRFCohgramBig = axes('Position', [col3Left inRFCohgramBtm cohgramW cohgramH]);
set(get(axInRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tDim - periDimWindow(1), f, dimInRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
% ylabel('Frequency (Hz)');

title(sprintf('InRF P%d Target Dim SFC (%d trials)', conditionsInfo.inRFLoc, conditionsInfo.numTrials.inRF));

%% exRF
axExRFCohgramBig = axes('Position', [col3Left exRFCohgramBtm cohgramW cohgramH]);
set(get(axExRFCohgramBig, 'Title'), 'Visible', 'on')

imagesc(tDim - periDimWindow(1), f, dimExRFC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(posCLim);
colorbar;
% ylabel('Frequency (Hz)');

title(sprintf('ExRF P%d Target Dim SFC (%d trials)', conditionsInfo.exRFLoc, conditionsInfo.numTrials.exRF));

%% diff
axDiffCohgramBig = axes('Position', [col3Left diffCohgramBtm cohgramW cohgramH]);
set(get(axDiffCohgramBig, 'Title'), 'Visible', 'on')

dimDiffC = dimInRFC - dimExRFC;

imagesc(tDim - periDimWindow(1), f, dimDiffC')
set(gca, 'YDir', 'normal');
hold on;
plot([0 0], [-10 max(f)+10], 'k-', 'LineWidth', 2); % plot line at window with center at event onset
caxis(diffCLim);
colormap(gca, getCoolWarmMap());
colorbar;
xlabel('Time from event onset (s)');
% ylabel('Frequency (Hz)');

title('InRF - ExRF Target Dim SFC');

%% alternately, just pass inRFC, exRFC, etc.

%% plot histogram of shuffled 

%% save
if ~isVisible
    set(gcf, 'Visible', 'off');
end
if ~isempty(saveFile)
    % note: running export_fig in parfor loop leads to cut off figs
	export_fig(saveFile, '-nocrop', '-r300');
end
if ~isVisible
    close;
end
