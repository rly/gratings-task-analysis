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
t = periFlashWindowOffset(1):1/Fs:periFlashWindowOffset(2)-1/Fs;

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session VEP Mapping Summary\n');
fprintf('Reference: %s\n', ref);

meanResponseAll = cell(nSessions, 1);

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
    
    % get mean response across flashes
    meanResponseAll{i} = mean(RR.responses, 3);
    % channel x time
    
    % don't do baseline correction or anything -- just the data as is
end

%% common average across all recordings and channels
commonAverageAll = nan(nSessions, numel(t));
for i = 1:nSessions
    assert(size(meanResponseAll{i}, 1) == 32);
    commonAverageAll(i,:) = mean(meanResponseAll{i}, 1);
end

superCommonAverage = mean(commonAverageAll, 1);

figure;
hold on;
plot(t, commonAverageAll);
plot(t, superCommonAverage, 'k', 'LineWidth', 3);
xlim(periFlashWindowOffset);

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
for s = 1:nSessions
    sessionInd = sessionInds(s);
    templateSessionInd = sessionInds(templateSessionIndInd);
    if s == templateSessionIndInd
        continue;
    end

    fprintf('\n');
    fprintf('Session %d vs Template Session %d:\n', sessionInd, templateSessionIndInd);
    
    % subtract superCommonAverage for visualization purposes
    % BUT it might need to be rescaled per session
    x1 = meanResponseAll{templateSessionIndInd};% - superCommonAverage;
    x2 = meanResponseAll{s};% - superCommonAverage;
    nTime = numel(t);
    nChannels = size(x1, 1);
    meanSDTime = (std(x1) + std(x2))/2;
    
    %% plot original
    figure_tr_inch(16, 9);

    origResponses = {x1 x2};
    plotHs = nan(2, 1);
    maxCAxis = -Inf;
    for k = 1:2
        plotHs(k) = subaxis(1, 2, k, 'MB', 0.1);
        hold on;

        imagesc(t, 1:nChannels, origResponses{k});
        set(gca, 'YDir', 'reverse');
        plot([0 0], [0 nChannels + 1], '-', 'Color', 0.3*ones(3, 1));
        xlim(periFlashWindowOffset);
        ylim([0.5 nChannels+0.5]);
        set(gca, 'YTickLabel', []);
        xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
        maxCAxis = max(abs(caxis));
        caxis([-maxCAxis maxCAxis]);
%         colormap(getCoolWarmMap());
        colorbar;
        set(gca, 'FontSize', 16);
        if k == 1
            title(sprintf('Session %d', templateSessionInd));
        elseif k == 2
            title(sprintf('Session %d', sessionInd));
        end
    end

%     for k = 1:2
%         caxis(plotHs(k), [-maxCAxis maxCAxis]);
%     end

    %% find shift that minimizes SSE
    minShift = -nChannels+1;
    maxShift = nChannels-1;
    shiftToTest = minShift:maxShift;
    nShiftToTest = numel(shiftToTest);
    nRandRefDiffsToTest = 1000;
    sigmaRefDiffToTest = meanSDTime; % add noise with variance over time equal to the mean across channels variance TODO check reasonable
    sseBest = Inf * ones(nShiftToTest, 1);
    nMatchesBest = nan(nShiftToTest, 1);
    refDiffBest = nan(nShiftToTest, nTime);
    for j = 1:nShiftToTest
        shift = shiftToTest(j);
        for k = 1:nRandRefDiffsToTest
            if k == 1
                refDiff = zeros(1, nTime); % make sure to test the 0 diff case
            else
                refDiff = sigmaRefDiffToTest .* randn(1, nTime); % 0 
            end
            sse = 0;
            nMatches = 0;
            for i = 1:nChannels
                if i - shift >= 1 && i - shift <= nChannels
                    sse = sse + sum((x1(i,:) - x2(i - shift,:) - refDiff).^2);
                    nMatches = nMatches + 1;
                end
            end
            if sse < sseBest(j)
                sseBest(j) = sse;
                nMatchesBest(j) = nMatches;
                refDiffBest(j,:) = refDiff;
            end
        end
        fprintf('j = %d, shift = %d, nMatches = %d, sse = %0.2f\n', j, shift, nMatchesBest(j), sseBest(j));
    end

    figure_tr_inch(8, 8);
    plot(shiftToTest, sseBest);
    xlabel('shift');
    ylabel('sse');
    title(sprintf('Session %d vs %d', templateSessionInd, sessionInd));

    %% plot shift
    shift = shiftBest(templateSessionIndInd,s);
    
    % plot ref diff for best shift
%     figure;
%     plot(t, refDiffBest(shiftToTest == shift,:));
%     xlim(periFlashWindowOffset);
    
    % plot shifted channel data side by side
    shiftResponse = cell(2, 1);
    newY = 1:nChannels+abs(shift);
    nChannelsNew = numel(newY);
    shiftResponse{1} = nan(nChannelsNew, nTime);
    shiftResponse{2} = nan(nChannelsNew, nTime);
    if shift <= 0 % a is deeper (larger number) than b
        shiftResponse{1}((1:nChannels)+abs(shift),:) = x1;
        shiftResponse{2}(1:nChannels,:) = x2;
    elseif shift > 0 % a is shallower (smaller number) than b
        shiftResponse{1}(1:nChannels,:) = x1;
        shiftResponse{2}((1:nChannels)+abs(shift),:) = x2;
    end
    

    figure_tr_inch(16, 9);
    
    plotHs = nan(2, 1);
    maxCAxis = -Inf;
    for k = 1:2
        plotHs(k) = subaxis(1, 2, k, 'MB', 0.1);
        hold on;

        imagesc(t, newY, shiftResponse{k});
        set(gca, 'YDir', 'reverse');
        plot([0 0], [0 nChannelsNew + 1], '-', 'Color', 0.3*ones(3, 1));
        xlim(periFlashWindowOffset);
        ylim([0.5 nChannelsNew+0.5]);
        set(gca, 'YTickLabel', []);
        xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
        maxCAxis = max(abs(caxis));
        caxis([-maxCAxis maxCAxis]);
%         colormap(getCoolWarmMap());
        colorbar;
        set(gca, 'FontSize', 16);
        if k == 1
            title(sprintf('Session %d', templateSessionInd));
        elseif k == 2
            title(sprintf('Session %d', sessionInd));
        end
    end

%     for k = 1:2
%         caxis(plotHs(k), [-maxCAxis maxCAxis]);
%     end
    
    drawnow;
end

%% compute average
% figure_tr_inch(20, 9);
maxCAxis = -Inf;
minShift = min(shiftBest(templateSessionIndInd,:));
maxShift = max(shiftBest(templateSessionIndInd,:));

[~,sessionIndsOrder] = sort(shiftBest(2,:));%[3 4 14 8 15 
assert(numel(sessionIndsOrder) == nSessions);

sumShiftResponses = zeros(maxShift - minShift + 32, numel(t));
countShiftResponses = zeros(maxShift - minShift + 32, 1);

for s = 1:nSessions
    so = sessionIndsOrder(s);
    sessionInd = sessionInds(s);
    
    shift = shiftBest(templateSessionIndInd,so);
    
    % plot shifted channel data side by side
    newY = (1+minShift):(32+maxShift);
    nChannelsNew = numel(newY);
    shiftResponse = nan(nChannelsNew, nTime);
    [~,newYInd] = intersect(newY, (1:nChannels)+shift);
    
    sumShiftResponses(newYInd,:) = sumShiftResponses(newYInd,:) + meanResponseAll{so} - superCommonAverage;
    countShiftResponses(newYInd) = countShiftResponses(newYInd) + 1;
end

meanShiftResponses = sumShiftResponses ./ countShiftResponses;

figure_tr_inch(4, 9);
subaxis(1, 1, 1, 'ML', 0.1);
hold on;
imagesc(t, newY, meanShiftResponses);
set(gca, 'YDir', 'reverse');
plot([0 0], newY([1 end]) + [-1 1], '-', 'Color', 0.3*ones(3, 1));
xlim(periFlashWindowOffset);
ylim(newY([1 end]) + [-0.5 0.5]);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
maxCAxis = 1/2 * max(abs(caxis));
caxis([-maxCAxis maxCAxis]);
%     colormap(getCoolWarmMap());
    colorbar;
    set(gca, 'FontSize', 16);
title(sprintf('Mean'));

%%
figure_tr_inch(20, 9);
plotHs = nan(nSessions, 1);
maxCAxis = -Inf;
minShift = min(shiftBest(templateSessionIndInd,:));
maxShift = max(shiftBest(templateSessionIndInd,:));

[~,sessionIndsOrder] = sort(shiftBest(templateSessionIndInd,:));%[3 4 14 8 15 
assert(numel(sessionIndsOrder) == nSessions);

for s = 1:nSessions
    so = sessionIndsOrder(s);
    sessionInd = sessionInds(s);
    
    shift = shiftBest(templateSessionIndInd,so);
    
    % plot shifted channel data side by side
    newY = (1+minShift):(32+maxShift);
    nChannelsNew = numel(newY);
    shiftResponse = nan(nChannelsNew, nTime);
    [~,newYInd] = intersect(newY, (1:nChannels)+shift);
    shiftResponse(newYInd,:) = meanResponseAll{so} - superCommonAverage;  
    
    plotHs(s) = subaxis(1, nSessions, s, 'SH', 0.001, 'ML', 0.02, 'MR', 0.02);
    hold on;

    imagesc(t, newY, shiftResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], newY([1 end]) + [-1 1], '-', 'Color', 0.3*ones(3, 1));
    xlim(periFlashWindowOffset);
    ylim(newY([1 end]) + [-0.5 0.5]);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
    maxCAxis = max(abs(caxis));
    caxis([-maxCAxis maxCAxis]);
%     colormap(getCoolWarmMap());
%     colorbar;
%     set(gca, 'FontSize', 16);
    title(sprintf('Sess %d', so));
end

%%
sessionIndsIncl = [15 2 5 11 10];
nSessionsIncl = numel(sessionIndsIncl);

figure_tr_inch(10, 7);
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
    shiftResponse(newYInd,:) = meanResponseAll{sessionIndsIncl(so)} - superCommonAverage;  
    shiftResponse(isnan(shiftResponse)) = -Inf; % make nans Inf for colormap
    plotHs(s) = subaxis(1, nSessionsIncl, s, 'SH', 0.004, 'ML', 0.04, 'MR', 0.02, 'MB', 0.11, 'MT', 0.06);
    hold on;

    imagesc(t, newY, shiftResponse);
    set(gca, 'YDir', 'reverse');
    newYInd([1 end])
    plot([0 0], newYInd([1 end]) + [-2.5 -1.5], '-', 'Color', 0.3*ones(3, 1));
    xlim([-0.05 0.25]);
    ylim(newY([1 end]) + [-0.5 0.5]);
    set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
    maxCAxis = max(abs(caxis))*1.1; % scale larger b/c max will show as white with nan=white mapping
    caxis([-maxCAxis maxCAxis]);
    colormap([1 1 1; colormap('parula')]); % make nans appear as white
%     colormap(getCoolWarmMap());
%     colorbar;
    set(gca, 'FontSize', 16);
    title(sprintf('Session %d', sessionIndsIncl(so)));
    
    set(gca, 'XTickLabel', []);
    if s == 5
        ax2 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'FontSize', 14);
        set(ax2, 'XLim', get(plotHs(s), 'XLim'), 'YLim', get(plotHs(s), 'YLim'));
        set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right', 'XAxisLocation', 'bottom', 'TickDir', 'out');
    end
    if s == 3
        xlabel({'','Time from Flash Onset (s)'});
    end
    if s == 1
        ylabel('(bottom)     ----------     Channel Index     ----------     (top)');
    end
end

plotFileName = sprintf('%s/ex-lfp-vepm-alignment-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');

%% plots in similar areas 
% one probe is 4.65 mm
% session 2, 3 channel 9 has anomalous low amplitude. is this a preprocessing issue or a hardware
% impedance issue??
% LFP data is preprocessed and outlier events removed. no CAR. this is RAW data. 
% session 5 also has a bunch of low amplitude signals

% sessionInd 4, 5, 16 next to 2, 3
% 16 was much more shallow (2.7 or 3.3 mm) than 4, 5
% sessionInd 7, 11, 15, 14, 20, 4, 5, 16 next to 2, 3. should have similar profiles. 
% need a mapping from first wave to second wave of recordings

%%
sessionIndsIncl = [4 5 16 2 3];
nSessionsIncl = numel(sessionIndsIncl);

figure_tr_inch(10, 7);
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
    shiftResponse(newYInd,:) = meanResponseAll{sessionIndsIncl(so)};% - superCommonAverage;  
    
%         [~,strongestNegativeCh] = min(min(shiftResponse, [], 2));
%         strongestNegativeCh = 29
%         meanAroundNegativeCh = mean(shiftResponse((strongestNegativeCh-1):(strongestNegativeCh+1),:));
%         shiftResponse = shiftResponse - meanAroundNegativeCh;
    
    shiftResponse(isnan(shiftResponse)) = -Inf; % make nans Inf for colormap
    
    plotHs(s) = subaxis(1, nSessionsIncl, s, 'SH', 0.004, 'ML', 0.04, 'MR', 0.02, 'MB', 0.11, 'MT', 0.06);
    hold on;

    imagesc(t, newY, shiftResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], newYInd([1 end]) + [-2.5 -1.5], '-', 'Color', 0.3*ones(3, 1));
    xlim([-0.05 0.25]);
    ylim(newY([1 end]) + [-0.5 0.5]);
    set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
    maxCAxis = max(abs(caxis))*1.1; % scale larger b/c max will show as white with nan=white mapping
    caxis([-maxCAxis maxCAxis]);
    colormap([1 1 1; colormap('parula')]); % make nans appear as white
%     colormap(getCoolWarmMap());
%     colorbar;
    set(gca, 'FontSize', 16);
    title(sprintf('Session %d', sessionIndsIncl(so)));
    
    set(gca, 'XTickLabel', []);
    if s == 5
        ax2 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'FontSize', 14);
        set(ax2, 'XLim', get(plotHs(s), 'XLim'), 'YLim', get(plotHs(s), 'YLim'));
        set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right', 'XAxisLocation', 'bottom', 'TickDir', 'out');
    end
    if s == 3
        xlabel({'','Time from Flash Onset (s)'});
    end
    if s == 1
        ylabel('(bottom)     ----------     Channel Index     ----------     (top)');
    end
end

plotFileName = sprintf('%s/ex-lfp-vepm-alignment-test-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');




