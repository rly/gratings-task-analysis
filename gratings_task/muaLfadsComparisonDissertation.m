%%
% clear;

% load LFADS output
load('C:\Users\Ryan\Documents\MATLAB\gratings-task-data\run180614_twoLocs_sixDays_170127to170311.mat');
lfads = data;
clear data;

unitInds = [4 7 11 16 27 30];

%%
sessionName = 'M20170127';
sessionInd = 1;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = unitInds;
muaChannelsToLoad = unitInds;
lfpChannelsToLoad = [];
lfpChannels = lfpChannelsToLoad;
numRandomizations = 2;
isLoadSortedSua = 0;
isLoadMua = 1;

v = 12;

outputDir = sprintf('%s/%s/', processedDataRootDir, 'MUA_LFADS/Dissertation');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

lfadsOut = lfads.lfadsOut{sessionInd};
lfadsWhichNeurons = lfads.whichNeurons{sessionInd};
clear lfads;

%%
taskName = 'GRATINGS';
scriptName = 'MUA_GRATINGS';
isZeroDistractors = 0;
[R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, taskName, scriptName, isLoadSortedSua, isLoadMua, 0, 0, 0);
sessionName = R.sessionName;
areaName = R.areaName;


%% plot example multi-units around cue onset
% TODO only do for P3
nUnits = numel(D.allUnitStructs);

binWidth = 0.008;
binBounds = -0.8:binWidth:0.8;
binTInds = binBounds > -0.1 & binBounds <= 0.5;

lfadsBinWidth = 0.004;
lfadsWindowOffset = [-0.1 0.5];
lfadsBinBounds = lfadsWindowOffset(1):lfadsBinWidth:lfadsWindowOffset(2);

loc = 3;

unitNames = cell(nUnits, 1);
saveFileNames = cell(nUnits, 1);
for j = 1:nUnits
    unitStruct = D.allUnitStructs{j};
    unitName = unitStruct.name;
    spikeTimes = unitStruct.ts;
    unitNames{j} = unitName;

    saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            processedDataDir, unitName, blockName, v);
    fprintf('Processing %s...\n', saveFileName);

    ES = load(saveFileName);
    
    nTrial = numel(ES.cueOnset.spikeTimesByLoc{loc});
    spikeMat = zeros(nTrial, numel(binBounds));
    for k = 1:nTrial
        spikeTs = ES.cueOnset.spikeTimesByLoc{3}(k).times;
        for l = 1:numel(spikeTs)
            spikeMat(k,ceil(spikeTs(l)/binWidth)) = spikeMat(k,ceil(spikeTs(l)/binWidth)) + 1;
        end
    end
    
    locTrialInds = find(ES.UE.cueLoc == loc);
    assert(numel(ES.cueOnset.spikeTimes) == numel(lfadsOut));
    assert(nTrial == numel(locTrialInds));
    
    m = find(lfadsWhichNeurons == unitInds(j));
    lfadsRateMat = zeros(nTrial, numel(lfadsBinBounds));
    for k = 1:nTrial
        l = locTrialInds(k);
        lfadsRateMat(k,:) = lfadsOut(l).rates(m, (lfadsOut(l).cueOnset + lfadsWindowOffset(1)*1000/4) : (lfadsOut(l).cueOnset + lfadsWindowOffset(2)*1000/4));
    end
    
    cols = lines(4);
    
    psthEst = mean(spikeMat(:,binTInds)) * 1/binWidth;
    psthRange = max(psthEst) - min(psthEst);
    
    xBounds = binBounds(binTInds);
    xBounds = xBounds([1 end]);
    yBounds = [min(psthEst)-psthRange/10 max(psthEst)+psthRange/10];
    
    figure_tr_inch(4, 4);
    ax1 = subaxis(2, 2, 1, 'SH', 0.08, 'MT', 0.08, 'MR', 0.03, 'ML', 0.23, 'MB', 0.17);
    plot(binBounds(binTInds), mean(spikeMat(:,binTInds)) * 1/binWidth, 'LineWidth', 1, 'Color', cols(1,:));
    xlim(xBounds);
    ylim(yBounds);
    ylabel({'Estimated', 'Firing Rate (Hz)'});
    title('Actual Data');
    set(gca, 'FontSize', 14);
    set(gca, 'XTickLabel', []);
    set(gca, 'TickDir', 'out');
    
    ax2 = subaxis(2, 2, 2);
    hold on;
    plot(binBounds(binTInds), mean(spikeMat(:,binTInds)) * 1/binWidth, 'LineWidth', 1, 'Color', cols(1,:));
    plot(lfadsBinBounds, mean(lfadsRateMat), 'LineWidth', 3, 'Color', cols(2,:));
    xlim(xBounds);
    ylim(yBounds);
    title('LFADS');
    set(gca, 'FontSize', 14);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickDir', 'out');

    subaxis(2, 2, 3);
    imagesc(binBounds(binTInds), 1:nTrial, spikeMat(:,binTInds));
    caxis([0 1]); % values may be >1 but for visualization just show binary
    xlim(xBounds);
    ylabel('Trial Number');
    xlabel('                          Time from Cue Onset (s)');
    set(gca, 'FontSize', 14);
    set(gca, 'YTick', 0:40:nTrial);
    set(gca, 'TickDir', 'out');

    subaxis(2, 2, 4);
    imagesc(lfadsBinBounds, 1:nTrial, lfadsRateMat);
    xlim(xBounds);
%     ylabel('Trial Number');
%     xlabel('Time from Cue Onset (s)');
    set(gca, 'FontSize', 14);
    set(gca, 'YTick', 0:40:nTrial);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickDir', 'out');
    
%     yBounds = [min([ylim(ax1) ylim(ax2)]) max([ylim(ax1) ylim(ax2)])];
%     ylim(ax1, yBounds);
%     ylim(ax2, yBounds);
    
%     suptitle(sprintf('LFADS Estimated Firing for Unit %s', unitNames{j}), 'FontSize', 18);
    
    % save as pdf
    saveFileNames{j} = sprintf('%s/%s-P%d-rawVsLfads-v%d.png', outputDir, unitNames{j}, loc, v);
    export_fig(saveFileNames{j}, '-nocrop');
end

%% collate all pdfs into one
% appendSaveFileName = sprintf('%s/%s-all-P%d-rawVsLfads-v%d.pdf', outputDir, sessionName, loc, v);
% deleteAndAppendPdfs(appendSaveFileName, saveFileNames);

