clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
logRoot = 'gratings_task_eye_tracking_';

sessionInd = 3;

nLoc = 4;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis: Spike Cross-Correlations\n');
fprintf('Loading %s...\n', pl2FilePath);
tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
if sessionInd == 1 || sessionInd == 2 % temp
    isLoadDirect = 1;
else
    isLoadDirect = 0;
end
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

%% remove spike and event times not during task to save memory
D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);

%%  
fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
logDir = sprintf('%s/%s', dataDir, sessionName(2:end));
load('params.mat');

% process events and sort them into different conditions
UE = getUsefulEvents2(logDir, gratingsTask3DLogIndices, 4, D);

% unload these variables to workspace
% (cueOnset, cueOnsetByLoc, ...
%         arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
%         arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%         arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%         targetDim, targetDimByLoc, ...
%         targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
%         nTrialShortHold, nTrialLongHold, rt);
% struct2var(usefulEvents);
% stop

totalTimeOverall = sum(D.blockStopTimes(gratingsTask3DIndices) - D.blockStartTimes(gratingsTask3DIndices));
minFiringRateOverall = 0.2;

inRFLoc = 1;
exRFLoc = 3;
nTrialsInRF = numel(UE.arrayOnsetByLoc{inRFLoc});
nTrialsExRF = numel(UE.arrayOnsetByLoc{exRFLoc});

nUnits = numel(D.allSpikeStructs);

%% create binned binary spike variable
origSpikeFs = D.timestampFrequency;
spikeFs = 500; % 1000 Hz, period 1 ms; 200 Hz, period 5 ms

preArrayOnsetSpikeTimesWindow = [0.4 0];

preArrayOnsetSpikeTimesInRF = cell(nUnits, 1);
preArrayOnsetSpikeTimesExRF = cell(nUnits, 1);

preArrayOnsetSpikeTimesInRFBinary = false(nUnits, nTrialsInRF, round(sum(preArrayOnsetSpikeTimesWindow)/spikeFs));
preArrayOnsetSpikeTimesExRFBinary = false(nUnits, nTrialsExRF, round(sum(preArrayOnsetSpikeTimesWindow)/spikeFs));

preArrayOnsetSpikeCountInRF = nan(nUnits, nTrialsInRF);
preArrayOnsetSpikeCountExRF = nan(nUnits, nTrialsExRF);

for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;
    
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    
    preArrayOnsetSpikeTimesInRF{i} = createnonemptydatamatpt(spikeVar, UE.arrayOnsetByLoc{inRFLoc}, preArrayOnsetSpikeTimesWindow);
    preArrayOnsetSpikeTimesExRF{i} = createnonemptydatamatpt(spikeVar, UE.arrayOnsetByLoc{exRFLoc}, preArrayOnsetSpikeTimesWindow);
  
    for j = 1:nTrialsInRF
        preArrayOnsetSpikeTimesInRFBinary(i,j,(ceil(spikeFs*preArrayOnsetSpikeTimesInRF{i}(j).times))) = true;
        preArrayOnsetSpikeCountInRF(i,j) = numel(preArrayOnsetSpikeTimesInRF{i}(j).times);
    end
    
    for j = 1:nTrialsExRF
        preArrayOnsetSpikeTimesExRFBinary(i,j,(ceil(spikeFs*preArrayOnsetSpikeTimesExRF{i}(j).times))) = true;
        preArrayOnsetSpikeCountExRF(i,j) = numel(preArrayOnsetSpikeTimesExRF{i}(j).times);
    end
end

%%
minSpikesInd2 = 50;
inds(1) = 56;
for k =     2:nUnits
    if k == inds(1)
        continue;
    end
    
    spikeStruct = D.allSpikeStructs{k};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;
    
    if ~(strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking'))
        continue;
    end
    if sum(sum(preArrayOnsetSpikeTimesInRFBinary(k,:,:))) < minSpikesInd2
        continue;
    end
    inds(2) = k;
    
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, k, ...
            nUnits, round(k/nUnits*100));
    
maxLag = 50; % units of 1/spikeFs
x = (-maxLag:maxLag)/spikeFs;
preArrayOnsetXCorrInRFAll = nan(nTrialsInRF, numel(x));
preArrayOnsetXCorrExRFAll = nan(nTrialsExRF, numel(x));

for j = 1:nTrialsInRF
    preArrayOnsetXCorrInRFAll(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(1),j,:)), ...
            squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(2),j,:)), maxLag);
end
for j = 1:nTrialsExRF
    preArrayOnsetXCorrExRFAll(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(1),j,:)), ...
            squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(2),j,:)), maxLag);
end

spikeTimeSeriesAll = false(2, round(D.totalTicks/(origSpikeFs/spikeFs)));
spikeTimeSeriesAll(1,(ceil(spikeFs*D.allSpikeStructs{inds(1)}.ts))) = true;
spikeTimeSeriesAll(2,(ceil(spikeFs*D.allSpikeStructs{inds(2)}.ts))) = true;

xcorrAll = xcorr(spikeTimeSeriesAll(1,:), spikeTimeSeriesAll(2,:), maxLag);

figure_tr_inch(10, 6);
set(gcf, 'Color', 'white');
h1 = subaxis(1, 2, 1);
box off;
hold on;
plot(x, mean(preArrayOnsetXCorrInRFAll), 'r');
plot(x, mean(preArrayOnsetXCorrExRFAll), 'b');
plot(x, movmean(mean(preArrayOnsetXCorrInRFAll), 5), 'r', 'LineWidth', 3);
plot(x, movmean(mean(preArrayOnsetXCorrExRFAll), 5), 'b', 'LineWidth', 3);
origYLim = ylim();
plot([0 0], [-1e10 1e10], 'Color', 0.3*ones(3, 1));
ylim(origYLim);

subaxis(1, 2, 2);
box off;
hold on;
plot(x, xcorrAll - mean(xcorrAll), 'Color', [0 0.7 0.3], 'LineWidth', 3);
origYLim = ylim();
plot([0 0], [-1e10 1e10], 'Color', 0.3*ones(3, 1));
ylim(origYLim);
suptitle(sprintf('Xcorr Unit %d and %d', inds));
% pause;
% close;

drawnow;

end
stop


%%
% account for increased firing rate -- one shuffle
nShuffles = 50;
preArrayOnsetXCorrInRFAllShuffleAll = nan(nShuffles, numel(x));
preArrayOnsetXCorrExRFAllShuffleAll = nan(nShuffles, numel(x));

for i = 1:nShuffles
    preArrayOnsetXCorrInRFAllShuffle = nan(nTrialsInRF, numel(x));
    preArrayOnsetXCorrExRFAllShuffle = nan(nTrialsExRF, numel(x));
    permTrialsInRF = randperm(nTrialsInRF);
    permTrialsExRF = randperm(nTrialsExRF);
    for j = 1:nTrialsInRF
        preArrayOnsetXCorrInRFAllShuffle(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(1),j,:)), ...
                squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(2),permTrialsInRF(j),:)), maxLag);
    end
    for j = 1:nTrialsExRF
        preArrayOnsetXCorrExRFAllShuffle(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(1),j,:)), ...
                squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(2),permTrialsExRF(j),:)), maxLag);
    end
    preArrayOnsetXCorrInRFAllShuffleAll(i,:) = mean(preArrayOnsetXCorrInRFAllShuffle);
    preArrayOnsetXCorrExRFAllShuffleAll(i,:) = mean(preArrayOnsetXCorrExRFAllShuffle);
end
plot(h1, x, movmean(mean(preArrayOnsetXCorrInRFAllShuffleAll), 5), 'm', 'LineWidth', 3);
plot(h1, x, movmean(mean(preArrayOnsetXCorrExRFAllShuffleAll), 5), 'g', 'LineWidth', 3);

