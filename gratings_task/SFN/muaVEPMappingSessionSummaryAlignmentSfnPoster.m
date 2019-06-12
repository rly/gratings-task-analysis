clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
sessionInds = [1:16 19:20 23];
v = 11;

fprintf('\n-------------------------------------------------------\n');
fprintf('VEP Mapping Summary - MUA\n');
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Session index: %d\n', sessionInds);
fprintf('------------------------\n');

%%
recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session VEP Mapping Summary\n');

meanResponseAll = cell(nSessions, 1);

%% session loop
allMUAStructs = cell(nSessions, 1);
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    muaChannelsToLoad = R.muaChannelsToLoad;
    blockName = strjoin(R.blockNames(R.vepmIndices), '-');
    processedDataDir = sprintf('%s/%s/MUA_VEPM/', processedDataRootDir, sessionName);
    allMUAStructs{i} = cell(numel(muaChannelsToLoad), 1);
    for j = 1:numel(muaChannelsToLoad)
        c = muaChannelsToLoad(j);
        saveFileName = sprintf('%s/%s_%s_%dM-%s-vepm-v%d.mat', processedDataDir, sessionName, areaName, c, blockName, v);
        
        fprintf('Loading file %s...\n', saveFileName);
        RR = load(saveFileName);
        allMUAStructs{i}{j} = RR.muaStruct;
    end
    spikeStruct = allMUAStructs{i}{1};
    t = spikeStruct.vepmPsthParams.t;
    meanResponseAll{i} = nan(numel(muaChannelsToLoad), numel(spikeStruct.vepmPsthParams.normPsthResponse));
    for j = 1:numel(muaChannelsToLoad)
        meanResponseAll{i}(j,:) = allMUAStructs{i}{j}.vepmPsthParams.normPsthResponse;
    end
end

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
    shiftResponse(newYInd,:) = meanResponseAll{sessionIndsIncl(so)};  
    shiftResponse(isnan(shiftResponse)) = -Inf; % make nans Inf for colormap
    plotHs(s) = subaxis(1, nSessionsIncl, s, 'SH', 0.004, 'ML', 0.04, 'MR', 0.1, 'MB', 0.13, 'MT', 0.06);
    hold on;

    imagesc(t, newY, shiftResponse);
    set(gca, 'YDir', 'reverse');
    plot([0 0], newYInd([1 end]) + [-2.5 -1.5], '-', 'Color', 0.3*ones(3, 1));
    xlim([-0.05 0.25]);
    ylim(newY([1 end]) + [-0.5 0.5]);
    set(gca, 'YTickLabel', []);
%     xlabel('Time from Flash Onset (s)');
%         maxCAxis = max([maxCAxis max(abs(caxis))]);
    maxCAxis = 18; %max(abs(caxis))*0.5 % scale larger b/c max will show as white with nan=white mapping
    caxis([-maxCAxis maxCAxis]);
%     colormap([1 1 1; colormap('parula')]); % make nans appear as white
    colormap([1 1 1; getCoolWarmMap()]);
%     colorbar;
    set(gca, 'FontSize', 18);
    title(sprintf('Session %d', sessionIndsIncl(so)));
    
    set(gca, 'XTickLabel', []);
    if s == 5
        ax2 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'FontSize', 18);
        set(ax2, 'XLim', get(plotHs(s), 'XLim'), 'YLim', get(plotHs(s), 'YLim'));
        set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right', 'XAxisLocation', 'bottom', 'TickDir', 'out');
        ax3 = axes('Position', get(plotHs(s), 'Position'), 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
        maxCAxis = 18; %max(abs(caxis))*0.5 % scale larger b/c max will show as white with nan=white mapping
        caxis([-maxCAxis maxCAxis]);
        cax = colorbar('FontSize', 18);
        cpos = get(cax, 'Position');
        set(cax, 'Position', [0.91 cpos(2:4)]);
        yl = ylabel(cax, 'Firing Rate (Baseline Z-Scored)', 'Rotation', -90, 'FontSize', 18);
        ylpos = get(yl, 'Position');
        set(yl, 'Position', ylpos + [1.7 0 0]);
    end
    if s == 3
        xlabel({'','Time from Flash Onset (s)'});
    end
    if s == 1
        ylabel('(bottom)   ----   Channel Index   ----   (top)');
    end
end

plotFileName = sprintf('%s/SFN/ex-mua-vepm-alignment-v%d.png', processedDataRootDir, v);
export_fig(plotFileName, '-nocrop');
