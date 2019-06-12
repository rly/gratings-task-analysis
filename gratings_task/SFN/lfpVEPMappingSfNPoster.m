clear;
sessionName = 'M20170201';
sessionInd = 3;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = 1:32;
muaChannelsToLoad = 1:32;
lfpChannelsToLoad = 1:32;
lfpChannels = lfpChannelsToLoad;
isZeroDistractors = 0;
numRandomizations = 2;

ref = 'CAR';
processedDataDir = sprintf('%s/%s/%s', processedDataRootDir, sessionName, 'LFP_VEPM');
blockName = 'vepm3';
areaName = 'PUL';
plotFileNamePrefix = sprintf('%s-ind%d-%s-ch%d-ch%d-%s', sessionName, sessionInd, areaName, channelsToLoad([1 end]), blockName);
Fs = 1000;
nChannels = 32;

v = 12;

%% save responses to mat file
saveFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDir, plotFileNamePrefix, ref, v);
load(saveFileName, 'D', 'R', 'responses', 't', 'periFlashWindowOffset', 'isNoisyChannel');

%%
nFlashes = size(responses, 3);


%% subtract out baseline
% units are standard deviations from baselined mean
averageResponse = mean(responses, 3); % average across flashes

% subtract mean baseline activity (-0.1, 0] seconds before flash
flashBaselineWindowOffset = [-0.1 0];

% TODO make so that baseline can have arbitrary end time
indexFlashTime = -round(periFlashWindowOffset(1) * Fs);
indexStartBaselineFlashTime = -round((periFlashWindowOffset(1) - flashBaselineWindowOffset(1)) * Fs) + 1;
for j = 1:nChannels
    averageResponse(j,:) = averageResponse(j,:) - mean(averageResponse(j,indexStartBaselineFlashTime:indexFlashTime));
end

%%
figure_tr_inch(4, 5.5);
ax1 = subaxis(1, 1, 1, 'ML', 0.07, 'MB', 0.14, 'MR', 0.06, 'MT', 0.03);
hold on;
cmap = getCoolWarmMap();
xBounds = [-0.05 0.251];

imagesc(t, 1:nChannels, averageResponse);

cBounds = max(abs(caxis)) * [-1 1] * 3/4;
caxis(cBounds); % symmetric cmap
colormap(cmap);
plot([0 0], [0.5 nChannels+0.5], '-', 'Color', 0.3*ones(3, 1)); 
set(gca, 'FontSize', 16);
set(gca, 'YDir', 'reverse');
xlabel('Time from Flash Onset (s)');
set(gca, 'XTick', -0.05:0.05:0.25);
set(gca, 'YTickLabel', '');
grid on;
xlim(xBounds);
ylim([0.5 nChannels+0.5]);
ax2 = axes('Position', get(ax1, 'Position'), 'Color', 'none', 'FontSize', 16);
set(ax2, 'XLim', get(ax1, 'XLim'), 'YLim', get(ax1, 'YLim'));
set(ax2, 'YTick', [], 'YColor', 'w', 'YAxisLocation', 'right');
set(ax2, 'XTick', get(ax1, 'XTick'), 'XTickLabel', [], 'XAxisLocation', 'bottom', 'TickDir', 'out');

plotFileName = sprintf('%s/%s-%s-lfpColorSfn-v%d.png', ...
        processedDataDir, plotFileNamePrefix, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

