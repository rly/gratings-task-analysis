clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

sessionInd = 1;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% setup and load data
fprintf('\n-------------------------------------------------------\n');
fprintf('LFP Correlation Structure Analysis\n');
fprintf('Loading %s...\n', pl2FileName);

tic;
isLoadSpikes = 0;
isLoadLfp = 1;
isLoadSpkc = 0;
isLoadDirect = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

%% correlate every channel with every other channel
nChannels = size(D.adjLfps, 1);
corrs = nan(nChannels, nChannels);
nanTimes = any(isnan(D.adjLfps), 1);
adjLfpsNoNan = D.adjLfps(:,~nanTimes);
for j = 1:nChannels - 1
    fprintf('Correlating channel %d with x...\n', j);
    for k = j + 1:nChannels
        corrs(j,k) = corr(adjLfpsNoNan(j,:)', adjLfpsNoNan(k,:)');
    end
end

corrsBothSides = triu(corrs, 1) + triu(corrs, 1)';
corrsBothSides(eye(size(corrsBothSides)) == 1) = NaN;

saveFileName = sprintf('%s/lfpCorrelationAnalysis-ch%d-to-ch%d-all.mat', ...
        processedDataDir, lfpChannelsToLoad([1 end]));
save(saveFileName, 'corrsBothSides');
% removing contamination from spikes does almost nothing at this level

%% plot
figure_tr_inch(12, 6);

subaxis(1, 2, 1);
bar(nanmean(corrsBothSides), 'FaceColor', lines(1));
xlim([0 size(corrsBothSides, 1) + 1]);
ylim([0 1]);
xlabel('Starting Channel Number');
ylabel('Mean Pearson Correlation With Other Channels');
set(gca, 'FontSize', 12);

subaxis(1, 2, 2);
imagesc(corrsBothSides');
xlabel('Channel Number');
ylabel('Channel Number');
set(gca, 'FontSize', 12);

suptitle(sprintf('%s - LFP Correlation Analysis - All Time Points', sessionName), ...
        'FontSize', 16);

plotFileName = sprintf('%s/lfpCorrelationAnalysis-ch%d-to-ch%d-all.png', ...
        processedDataDir, lfpChannelsToLoad([1 end]));
export_fig(plotFileName, '-nocrop');