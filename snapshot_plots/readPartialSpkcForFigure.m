sessionInd = 9;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
tic;
isLoadSpikes = 0;
isLoadLfp = 0;
isLoadSpkc = 1;
isLoadDirect = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

assert(numel(blockNames) == numel(D.blockStartTimes));
blockName = strjoin(blockNames(gratingsTask3DIndices), '-');
fprintf('Analyzing block names: %s.\n', blockName);

%% plot spkcs
chsToPlot = 1:5:32;
yScale = 8;
yOffset = 1;
tInd = (1:20000) + 40000*2;

figure_tr_inch(4, 7);
subaxis(1, 1, 1, 'MT', 0.02, 'MB', 0.02, 'ML', 0.15, 'MR', 0.03);
hold on;
box off;
for i = 1:numel(chsToPlot)
    plot(D.spkcs(chsToPlot(i),tInd)*yScale + i*yOffset, 'LineWidth', 2, 'Color', 0.2 * ones(3, 1));
end
set(gca, 'YDir', 'reverse');
ylim([1-yOffset+0.5 numel(chsToPlot) + yOffset-0.5]);
% ylabel('Channel Number');
set(gca, 'XTick', []);
set(gca, 'XColor', 'none');
set(gca, 'YTick', 1:numel(chsToPlot));
set(gca, 'YTickLabel', chsToPlot);
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'LineWidth', 2);
set(gcf, 'Color', 'white');

% plotFileName = sprintf('%s-spkc_series.png', sessionName);
% export_fig(plotFileName, '-nocrop');

print(sprintf('%s/posterFigs/%s-spkc-sample.svg', processedDataRootDir, sessionName), '-dsvg');
