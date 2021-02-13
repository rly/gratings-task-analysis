% function lfpVEPMappingSummary(processedDataRootDir, recordingInfoFileName, sessionInds, ref)

clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
ref = 'RAW';
sessionInds = [1:16 19:20 23];

fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Summary - LFP\n');
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Reference: %s\n', ref);
fprintf('Session index: %d\n', sessionInds);
fprintf('------------------------\n');

v = 12;

%%
recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

outputDir = sprintf('%s/%s/', processedDataRootDir, 'LFP_GRATINGS_SUMMARY');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

% TEMP until it is saved
periFlashWindowOffset = [-0.25 0.3]; % seconds around flash
Fs = 1000;

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session VEP Mapping Summary\n');
fprintf('Reference: %s\n', ref);

meanResponseAll = cell(nSessions, 1);

params.tapers = [7 13];
params.pad = 2;
params.Fs = 1000;
params.fpass = [70 200];
params.err = [0 0];
params.trialave = 1;
movingWin = [0.1 0.025];

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    lfpChannelsToLoad = R.lfpChannelsToLoad;
    blockName = strjoin(R.blockNames(R.vepmIndices), '-');
    processedDataDir = sprintf('%s/%s/LFP_VEPM/', processedDataRootDir, sessionName);
    fileNamePrefix = sprintf('%s-ind%d-%s-ch%d-ch%d-%s', sessionName, sessionInd, areaName, lfpChannelsToLoad([1 end]), blockName);
    saveFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDir, fileNamePrefix, ref, v);
    fprintf('Loading file %s...\n', saveFileName);
    RR = load(saveFileName);

    % channel x time
    meanResponseAll{i} = nan(size(RR.responses, 1), 19);
    
    % get mean response across flashes
    for j = 1:numel(lfpChannelsToLoad)
        [S,t,f] = mtspecgramc(squeeze(RR.responses(j,:,:)), movingWin, params);
        meanS = mean(S, 2);
        meanResponseAll{i}(j,:) = meanS;
    end
end

%% common average across all recordings and channels
t = t + periFlashWindowOffset(1);
baselineInd = t < 0;
commonAverageAll = nan(nSessions, numel(t));
for i = 1:nSessions
    assert(size(meanResponseAll{i}, 1) == 32);
    commonAverageAll(i,:) = mean(meanResponseAll{i}, 1);
end

superCommonAverage = mean(commonAverageAll, 1);

% note: subtracting the super common average later can flip the sign for
% some sessions and keep the sign for others... this won't matter for the
% fitting procedure but for visualization

%% brute-force search against template

% sessionInd 9 is very noisy and unusual
% sessionInd 12, 14, 15, 16, 17, 19 are very unusual
% especially 16

% positive shift means test (dim1) is deeper than template (dim2)
shiftBest = zeros(nSessions, nSessions);
shiftBest(2,1) = 15;
shiftBest(2,2) = 0;
shiftBest(2,3) = -9;
shiftBest(2,4) = -9;
shiftBest(2,5) = 0;
shiftBest(2,6) = 3;
shiftBest(2,7) = 12;
shiftBest(2,8) = -3;
shiftBest(2,9) = 12; % shift might be less - look at pre-cue activity
shiftBest(2,10) = 7;
shiftBest(2,11) = 5;
shiftBest(2,12) = 0; % not sure
shiftBest(2,13) = 5;
shiftBest(2,14) = -7; % not sure
shiftBest(2,15) = -2;
shiftBest(2,16) = 0; % don't know
shiftBest(2,17) = 11;
shiftBest(2,18) = 6;
shiftBest(2,19) = 11;

templateSessionIndInd = 2; % index into sessionInds (compare to s)

nTime = numel(t);
nChannels = size(meanResponseAll{templateSessionIndInd}, 1);

%%
sessionIndsIncl = [15 2 5 11 10];
nSessionsIncl = numel(sessionIndsIncl);

figure_tr_inch(12, 7);
plotHs = nan(nSessionsIncl, 1);
maxCAxis = -Inf;

minShift = min(shiftBest(templateSessionIndInd,sessionIndsIncl));
maxShift = max(shiftBest(templateSessionIndInd,sessionIndsIncl));

[~,sessionIndsOrder] = sort(shiftBest(templateSessionIndInd,sessionIndsIncl));
assert(numel(sessionIndsOrder) == nSessionsIncl);

for s = 1:nSessionsIncl
    so = sessionIndsOrder(s);
    sessionInd = sessionIndsIncl(s);
    
    shift = shiftBest(templateSessionIndInd,sessionIndsIncl(so));
    
    % plot shifted channel data side by side
    newY = (1+minShift):(32+maxShift);
    nChannelsNew = numel(newY);
    shiftResponse = nan(nChannelsNew, nTime);
    [~,newYInd] = intersect(newY, (1:nChannels)+shift);
    shiftResponse(newYInd,:) = meanResponseAll{sessionIndsIncl(so)} - mean(meanResponseAll{1}(:,baselineInd), 2) - superCommonAverage;  
    shiftResponse(isnan(shiftResponse)) = -Inf; % make nans Inf for colormap
    plotHs(s) = subaxis(1, nSessionsIncl, s, 'SH', 0.004, 'ML', 0.04, 'MR', 0.1, 'MB', 0.11, 'MT', 0.06);
    hold on;

    imagesc(t, newY, shiftResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], newYInd([1 end]) + [-2.5 -1.5], '-', 'Color', 0.3*ones(3, 1));
    xlim([-0.05 0.25]);
    ylim(newY([1 end]) + [-0.5 0.5]);
    set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
    maxCAxis = max(abs(caxis)); % scale larger b/c max will show as white with nan=white mapping
    caxis([-0.25 0.015] * 1e-3);
%     colormap([1 1 1; colormap('parula')]); % make nans appear as white
%     colormap([1 1 1; getCoolWarmMap()]);
colormap(getCoolWarmMap());
%     colorbar;
    set(gca, 'FontSize', 18);
    title(sprintf('Session %d', sessionIndsIncl(so)));
    
    set(gca, 'XTickLabel', []);
    if s == 5
        ax2 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'FontSize', 18);
        set(ax2, 'XLim', get(plotHs(s), 'XLim'), 'YLim', get(plotHs(s), 'YLim'));
        set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right', 'XAxisLocation', 'bottom', 'TickDir', 'out');
        ax3 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
%         maxCAxis = max(abs(caxis(plotHs(s))))*0.75; %max(abs(caxis))*0.5 % scale larger b/c max will show as white with nan=white mapping
        caxis([-0.25 0.015] * 1e-3);
        cax = colorbar('FontSize', 18);
        set(cax, 'YTick', -1:0.5:1);
        cpos = get(cax, 'Position');
        set(cax, 'Position', [0.91 cpos(2:4)]);
        yl = ylabel(cax, 'LFP Voltage (Baseline Z-Scored, Grand Avg Ref)', 'Rotation', -90, 'FontSize', 18);
        ylpos = get(yl, 'Position');
        set(yl, 'Position', ylpos + [1.5 0 0]);
    end
    if s == 3
        xlabel({'','Time from Flash Onset (s)'});
    end
    if s == 1
        ylabel('(bottom)    --------    Channel Index    --------    (top)');
    end
end

plotFileName = sprintf('%s/SFN/ex-bhg-vepm-alignment-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');