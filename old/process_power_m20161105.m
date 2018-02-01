% evt5 = cue onset
% evt1,2,3,4 = cue offset
% evt6 = array onset and also release shape resp win start
% evt7 = target dim
% evt8 = juice

% 1105 g3-g6, ch65-69 above brain
% 1104 g3-g6, ch65-66 above brain
% 1102 g2-g9, ch65-68 above brain
% 1101 g2-g3, ch65-67 above brain
clear;

sessionRange = 8;
numSessions = numel(sessionRange);

% dim 1: session name
% dim 2: log indices
% dim 3: linear array channel offset from dura
sessions = {...
        '20161027', 3:6, 0; % also g6
        '20161028', 4:7, 0;
        '20161031', 1:4, 0;
        '20161101', 2:3, 3;
        '20161102', 2:9, 4;
        '20161104', 3:6, 2;
        '20161105', 3:6, 5;
        '20170127', 4:7, 0;
        };

lfpVarName = 'FP001';

numFreqInS256To511 = 92;
numFreqInS128To255 = 46;
SPreArrayInRFAll = nan(numSessions, numFreqInS256To511);
SPreArrayExRFAll = nan(numSessions, numFreqInS256To511);
SPreDimInRFAll = nan(numSessions, numFreqInS256To511);
SPreDimExRFAll = nan(numSessions, numFreqInS256To511);
SPreDimShortHoldInRFAll = nan(numSessions, numFreqInS256To511);
SPreDimShortHoldExRFAll = nan(numSessions, numFreqInS256To511);
SPreDimLongHoldInRFAll = nan(numSessions, numFreqInS256To511);
SPreDimLongHoldExRFAll = nan(numSessions, numFreqInS256To511);
SPreCueInRFAll = nan(numSessions, numFreqInS256To511);
SPreCueAllTrialsAll = nan(numSessions, numFreqInS256To511);

nChannels = 32;
startChannel = 1;

inRFLoc = 3;
exRFLoc = 1;

nLoc = 4;

outputDir = 'processed_data';

%%
for m = 1:numSessions
    l = sessionRange(m);
    sessionName = sessions{l,1};
    logIndices = sessions{l,2};
    xChannelAlignmentOffset = sessions{l,3};
    fprintf('Processing %s...\n', sessionName);

    dataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-data/%s/', ...
            sessionName);
    dataFile = sprintf('%s/%s-gmerged-fp.mat', dataDir, sessionName);
    load(dataFile);
    load('params.mat')
    
%     channelRange = xChannelAlignmentOffset:
    
%     for m = channelRange

%     deepChannelOrig = 27; % or 4
%     deepChannel = deepChannelOrig + xChannelAlignmentOffset;
%     lfpVarName = sprintf('FP%03d', startChannel + deepChannel - 1); % overwrite

    fprintf('\tUsing LFP var: %s\n', lfpVarName);

%% process events and sort them into different conditions
usefulEvents = getUsefulEvents(dataDir, logIndices, nLoc, ...
        EVT02, EVT04, EVT06, EVT07, EVT14, EVT15, EVT16, EVT08);

% unload these variables to workspace
% (cueOnset, cueOnsetByLoc, ...
%         arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
%         arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%         arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%         targetDim, targetDimByLoc, ...
%         targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
%         nTrialShortHold, nTrialLongHold, rt);
struct2var(usefulEvents);

%% LFPs
lfpTs = eval([lfpVarName '_ts']);
lfpTsStep = eval([lfpVarName '_ts_step']);
lfpInd = eval([lfpVarName '_ind']);
Fs = round(1/lfpTsStep);
adjLfp = padNaNsToAdjustLfpOffset(eval(lfpVarName), lfpTs, lfpInd, Fs);
% there are nan's in the data at certain events!!

%%
cueOnsetLfpWindow = [1 1.5];
cueOnsetLfpT = -cueOnsetLfpWindow(1) * Fs : cueOnsetLfpWindow(2) * Fs - 1;
lfpAroundCueOnset = createdatamatc(adjLfp, cueOnset, Fs, cueOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundCueOnset))));

arrayOnsetLfpWindow = [1.5 1];
arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
arrayOnsetLfpT = -arrayOnsetLfpWindow(1) * Fs : arrayOnsetLfpWindow(2) * Fs - 1;
lfpAroundArrayOnset = createdatamatc(adjLfp, arrayOnset, Fs, arrayOnsetLfpWindow);
assert(~any(any(isnan(lfpAroundArrayOnset))));

targetDimLfpWindow = [1.5 1];
targetDimLfpXLim = [-0.7 0.4] * Fs;
targetDimLfpT = -targetDimLfpWindow(1) * Fs : targetDimLfpWindow(2) * Fs - 1;
lfpAroundTargetDim = createdatamatc(adjLfp, targetDim, Fs, targetDimLfpWindow);
assert(~any(any(isnan(lfpAroundTargetDim))));

lfpAroundCueOnsetByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLoc = cell(nLoc, 1);
lfpAroundArrayOnsetRelByLoc = cell(nLoc, 1);
lfpAroundTargetDimByLoc = cell(nLoc, 1);
lfpAroundTargetDimShortHoldByLoc = cell(nLoc, 1);
lfpAroundTargetDimLongHoldByLoc = cell(nLoc, 1);

for i = 1:nLoc
    lfpAroundCueOnsetByLoc{i} = createdatamatc(adjLfp, cueOnsetByLoc{i}, Fs, cueOnsetLfpWindow);
    lfpAroundArrayOnsetHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundArrayOnsetShortHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetShortHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundArrayOnsetLongHoldByLoc{i} = createdatamatc(adjLfp, arrayOnsetLongHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundArrayOnsetRelByLoc{i} = createdatamatc(adjLfp, arrayOnsetRelByLoc{i}, Fs, arrayOnsetLfpWindow);
    lfpAroundTargetDimByLoc{i} = createdatamatc(adjLfp, targetDimByLoc{i}, Fs, targetDimLfpWindow);
    lfpAroundTargetDimShortHoldByLoc{i} = createdatamatc(adjLfp, targetDimShortHoldByLoc{i}, Fs, targetDimLfpWindow);
    lfpAroundTargetDimLongHoldByLoc{i} = createdatamatc(adjLfp, targetDimLongHoldByLoc{i}, Fs, targetDimLfpWindow);
end

%%
load('params.mat');
params.fpass = [5 50];
params.pad = 2;
params.tapers = [1 1];

%%
cueOnsetLfpSelection = cueOnsetLfpT >= -300 & cueOnsetLfpT <= 0;
[SPreCueInRF,fSpectrumCue,SErr] = mtspectrumc(lfpAroundCueOnsetByLoc{inRFLoc}(cueOnsetLfpSelection,:), params); % unused
% pre-cue response for all trials, all locations, all attend conditions
[SPreCueAllTrials,~,~] = mtspectrumc(lfpAroundCueOnset(cueOnsetLfpSelection,:), params);

% InRF
arrayOnsetLfpSelection = arrayOnsetLfpT >= -400 & arrayOnsetLfpT <= 0;
[SPreArrayInRF,fSpectrum,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{inRFLoc}(arrayOnsetLfpSelection,:), params);

arrayOnsetLfpSelection = arrayOnsetLfpT >= 300 & arrayOnsetLfpT <= 700;
[SPostArrayInRF,~,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{inRFLoc}(arrayOnsetLfpSelection,:), params); % unused

targetDimLfpSelection = targetDimLfpT >= -400 & targetDimLfpT <= 0;
[SPreDimInRF,~,SErr] = mtspectrumc(lfpAroundTargetDimByLoc{inRFLoc}(targetDimLfpSelection,:), params);
[SPreDimShortHoldInRF,~,SErr] = mtspectrumc(lfpAroundTargetDimShortHoldByLoc{inRFLoc}(targetDimLfpSelection,:), params);
[SPreDimLongHoldInRF,~,SErr] = mtspectrumc(lfpAroundTargetDimLongHoldByLoc{inRFLoc}(targetDimLfpSelection,:), params);

% ExRF
arrayOnsetLfpSelection = arrayOnsetLfpT >= -400 & arrayOnsetLfpT <= 0;
[SPreArrayExRF,~,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{exRFLoc}(arrayOnsetLfpSelection,:), params);

arrayOnsetLfpSelection = arrayOnsetLfpT >= 300 & arrayOnsetLfpT <= 700;
[SPostArrayExRF,~,SErr] = mtspectrumc(lfpAroundArrayOnsetHoldByLoc{exRFLoc}(arrayOnsetLfpSelection,:), params); % unused

targetDimLfpSelection = targetDimLfpT >= -400 & targetDimLfpT <= 0;
[SPreDimExRF,~,SErr] = mtspectrumc(lfpAroundTargetDimByLoc{exRFLoc}(targetDimLfpSelection,:), params);
[SPreDimShortHoldExRF,~,SErr] = mtspectrumc(lfpAroundTargetDimShortHoldByLoc{exRFLoc}(targetDimLfpSelection,:), params);
[SPreDimLongHoldExRF,~,SErr] = mtspectrumc(lfpAroundTargetDimLongHoldByLoc{exRFLoc}(targetDimLfpSelection,:), params);

%%
if strcmp(lfpVarName, 'FP098')
    yAxis = [-65 -45];
else
    yAxis = [-35 -10];
end

figure;
hold on;
plot_vector(SPreCueAllTrials, fSpectrumCue, 'l', [], 'k');
plot_vector(SPreArrayInRF, fSpectrum, 'l', [], 'r');
plot_vector(SPreArrayExRF, fSpectrum, 'l', [], 'b');
plot_vector(SPreDimInRF, fSpectrum, 'l', [], 'm');
plot_vector(SPreDimExRF, fSpectrum, 'l', [], 'c');
xlim(params.fpass);
ylim(yAxis);
legend({'Pre-Cue Baseline', ...
        'Pre-Stimulus Delay, Attend In RF', ...
        'Pre-Stimulus Delay, Attend Out RF', ...
        'Post-Stimulus Delay, Attend In RF', ...
        'Post-Stimulus Delay, Attend Out RF', ...
        }, 'Location', 'NorthEast');
plotFileName = sprintf('%s-%s-powerComparison-%0.1f-%0.1fHz.png', ...
        sessionName, lfpVarName, params.fpass);
export_fig(plotFileName);

SPreArrayInRFAll(m,:) = SPreArrayInRF;
SPreArrayExRFAll(m,:) = SPreArrayExRF;
SPreDimInRFAll(m,:) = SPreDimInRF;
SPreDimExRFAll(m,:) = SPreDimExRF;
SPreDimShortHoldInRFAll(m,:) = SPreDimShortHoldInRF;
SPreDimShortHoldExRFAll(m,:) = SPreDimShortHoldExRF;
SPreDimLongHoldInRFAll(m,:) = SPreDimLongHoldInRF;
SPreDimLongHoldExRFAll(m,:) = SPreDimLongHoldExRF;
SPreCueInRFAll(m,:) = SPreCueInRF;
SPreCueAllTrialsAll(m,:) = SPreCueAllTrials;

end % end numSessions

%% summary figs across sessions
figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
h0 = plot(fSpectrumCue, 10*log10(mean(SPreCueAllTrialsAll, 1)), ...
        'Color', [0 0 0], ...
        'LineWidth', 4);
h1 = plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
% plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
h2 = plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
% plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) + std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) - std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
h3 = plot(fSpectrum, 10*log10(mean(SPreDimInRFAll, 1)), ...
        'Color', [1 0 1], ...
        'LineWidth', 4);
h4 = plot(fSpectrum, 10*log10(mean(SPreDimExRFAll, 1)), ...
        'Color', [0 1 1], ...
        'LineWidth', 4);
% h1 = plot(fSpectrum, zeros(size(fSpectrum)), 'Color', [1 1 1]);
% h2 = plot(fSpectrum, zeros(size(fSpectrum)), 'Color', [1 1 1]);
% h3 = plot(fSpectrum, zeros(size(fSpectrum)), 'Color', [1 1 1]);
% h4 = plot(fSpectrum, zeros(size(fSpectrum)), 'Color', [1 1 1]);

xlim(params.fpass);
ylim(yAxis);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend([h0 h1 h2 h3 h4], ...
        {sprintf('\\color[rgb]{%s}Pre-Cue Baseline', num2str(h0.Color)), ...
        sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend In RF', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend Out RF', num2str(h2.Color)), ...
        sprintf('\\color[rgb]{%s}Post-Stimulus, Attend In RF', num2str(h3.Color)), ...
        sprintf('\\color[rgb]{%s}Post-Stimulus, Attend Out RF', num2str(h4.Color))}, ...
        'Position', [0.5 0.78 0.4 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s/meanSessions-%s-powerComparison-%0.1f-%0.1fHz-n%d.png', ...
        outputDir, lfpVarName, params.fpass, numel(sessionRange));
export_fig(plotFileName);

%%
yAxisPercent = [-100 100];

figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
baselinePower = mean(SPreCueAllTrialsAll, 1);
plot(params.fpass, [0 0], 'k-');
h1 = plot(fSpectrum, (mean(SPreArrayInRFAll, 1) - baselinePower)./baselinePower*100, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
% plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
h2 = plot(fSpectrum, (mean(SPreArrayExRFAll, 1) - baselinePower)./baselinePower*100, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
% plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) + std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) - std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
h3 = plot(fSpectrum, (mean(SPreDimInRFAll, 1) - baselinePower)./baselinePower*100, ...
        'Color', [1 0 1], ...
        'LineWidth', 4);
h4 = plot(fSpectrum, (mean(SPreDimExRFAll, 1) - baselinePower)./baselinePower*100, ...
        'Color', [0 1 1], ...
        'LineWidth', 4);
xlim(params.fpass);
ylim(yAxisPercent);
xlabel('Frequency (Hz)');
ylabel('Relative Power Change (%)');
legend([h1 h2 h3 h4], ...
        {sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend In RF', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend Out RF', num2str(h2.Color)), ...
        sprintf('\\color[rgb]{%s}Post-Stimulus, Attend In RF', num2str(h3.Color)), ...
        sprintf('\\color[rgb]{%s}Post-Stimulus, Attend Out RF', num2str(h4.Color))}, ...
        'Position', [0.5 0.78 0.4 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s/meanSessions-%s-powerComparisonPercentChange-%0.1f-%0.1fHz-n%d.png', ...
        outputDir, lfpVarName, params.fpass, numel(sessionRange));
export_fig(plotFileName);

%%
% close all
% save(sprintf('%s/%sworkspace-%s-powerComparison-%0.1f-%0.1fHz-n%d.mat', ...
%         outputDir, datestr(now, 'yyyymmdd'), lfpVarName, params.fpass, numel(sessionRange)));
