clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
logRoot = 'gratings_task_eye_tracking_';

sessionInd = 2;

nLoc = 4;

% reference:
% ftp://nozdr.ru/biblio/kolxo3/B/BH/Gruen%20S.,%20Rotter%20S.%20(eds.)%20Analysis%20of%20Parallel%20Spike%20Trains%20(Springer,%202010)(ISBN%201441956743)(O)(462s)_BH_.pdf#page=212

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis: Spike Unitary Event Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
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

inRFLoc = 1;
exRFLoc = 3;
nTrialsInRF = numel(UE.arrayOnsetByLoc{inRFLoc});
nTrialsExRF = numel(UE.arrayOnsetByLoc{exRFLoc});

nUnits = numel(D.allSpikeStructs);

%% create binned binary spike variable
origSpikeFs = D.timestampFrequency;
spikeFs = 1000; % 1000 Hz, period 1 ms; 200 Hz, period 5 ms

preArrayOnsetSpikeTimesWindow = [0.5 0];

preArrayOnsetSpikeTimesInRF = cell(nUnits, 1);
preArrayOnsetSpikeTimesExRF = cell(nUnits, 1);

nTime = round(sum(preArrayOnsetSpikeTimesWindow)*spikeFs);
preArrayOnsetSpikeTimesInRFBinary = false(nUnits, nTrialsInRF, nTime);
preArrayOnsetSpikeTimesExRFBinary = false(nUnits, nTrialsExRF, nTime);

preArrayOnsetSpikeCountInRF = nan(nUnits, nTrialsInRF);
preArrayOnsetSpikeCountExRF = nan(nUnits, nTrialsExRF);

unitNames = cell(nUnits);

for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;
    
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    
    unitNames{i} = unitName;
        
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
isCell = false(nUnits, 1);
isFiringRateHigh = false(nUnits, 1);
minFiringRate = 2;
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    isCell(i) = strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking');
    meanFiringRate = sum(preArrayOnsetSpikeCountInRF(i,:)) / nTrialsInRF / sum(preArrayOnsetSpikeTimesWindow); % mean spikes/trial/sec
    isFiringRateHigh(i) = meanFiringRate > minFiringRate;
end

selection = isCell & isFiringRateHigh;
nCells = sum(selection);
fprintf('Looking at %d cells:\n', nCells);
fprintf('%s\n', unitNames{selection});
preArrayOnsetSpikeTimesInRFBinaryCells = preArrayOnsetSpikeTimesInRFBinary(selection,:,:);
preArrayOnsetSpikeTimesExRFBinaryCells = preArrayOnsetSpikeTimesExRFBinary(selection,:,:);


binWidth = 5; % 5 ms
binStartIndex = 1:binWidth:nTime - binWidth + 1;
nBins = numel(binStartIndex);

preArrayOnsetWindowedFiringInRFBinaryCells = false(nCells, nTrialsInRF, nBins);
preArrayOnsetWindowedFiringExRFBinaryCells = false(nCells, nTrialsExRF, nBins);

for i = 1:nCells
    for n = 1:nBins
        binIndex = binStartIndex(n):binStartIndex(n) + binWidth - 1;
        for j = 1:nTrialsInRF
            preArrayOnsetWindowedFiringInRFBinaryCells(i,j,n) = any(preArrayOnsetSpikeTimesInRFBinaryCells(i,j,binIndex));
        end
        for j = 1:nTrialsExRF
            preArrayOnsetWindowedFiringExRFBinaryCells(i,j,n) = any(preArrayOnsetSpikeTimesExRFBinaryCells(i,j,binIndex));
        end
    end
end

sumFiringRateEachTrialInRF = nan(nCells, nTrialsInRF);
for i = 1:nCells
    for j = 1:nTrialsInRF
        sumFiringRateEachTrialInRF(i,j) = sum(preArrayOnsetWindowedFiringInRFBinaryCells(i,j,:)) / nBins;
    end
    for j = 1:nTrialsExRF
        sumFiringRateEachTrialExRF(i,j) = sum(preArrayOnsetWindowedFiringExRFBinaryCells(i,j,:)) / nBins;
    end
end
probFiringInBinInRF = mean(sumFiringRateEachTrialInRF, 2); % mean over trials
probFiringInBinExRF = mean(sumFiringRateEachTrialExRF, 2); % mean over trials

% for each bin n, find unique firing patterns k
allCoincidencePatternsInRF = false(nBins * nTrialsInRF, nCells);
patternCount = 0;
% TODO vectorize/use reshape
for j = 1:nTrialsInRF
    for n = 1:nBins
        patternCount = patternCount + 1;
        allCoincidencePatternsInRF(patternCount,:) = preArrayOnsetWindowedFiringInRFBinaryCells(:,j,n);
    end
end

allCoincidencePatternsExRF = false(nBins * nTrialsExRF, nCells);
patternCount = 0;
% TODO vectorize/use reshape
for j = 1:nTrialsExRF
    for n = 1:nBins
        patternCount = patternCount + 1;
        allCoincidencePatternsExRF(patternCount,:) = preArrayOnsetWindowedFiringExRFBinaryCells(:,j,n);
    end
end

uniqueCoincidencePatternsInRF = unique(allCoincidencePatternsInRF, 'rows');
nUniquePatternsInRF = size(uniqueCoincidencePatternsInRF, 1);
actualCountUniquePatternsInRF = nan(nUniquePatternsInRF, 1);
fprintf('Found %d unique firing patterns InRF.\n', nUniquePatternsInRF);

uniqueCoincidencePatternsExRF = unique(allCoincidencePatternsExRF, 'rows');
nUniquePatternsExRF = size(uniqueCoincidencePatternsExRF, 1);
actualCountUniquePatternsExRF = nan(nUniquePatternsExRF, 1);
fprintf('Found %d unique firing patterns ExRF.\n', nUniquePatternsExRF);

% count instances of unique patterns over trials and time
for k = 1:nUniquePatternsInRF
    [~,matchIndex] = ismember(allCoincidencePatternsInRF, uniqueCoincidencePatternsInRF(k,:), 'rows');
    actualCountUniquePatternsInRF(k) = sum(matchIndex);
end
for k = 1:nUniquePatternsExRF
    [~,matchIndex] = ismember(allCoincidencePatternsExRF, uniqueCoincidencePatternsExRF(k,:), 'rows');
    actualCountUniquePatternsExRF(k) = sum(matchIndex);
end

uniqueCoincidencePatternsAsProbInRF = nan(size(uniqueCoincidencePatternsInRF));
uniqueCoincidencePatternsAsProbExRF = nan(size(uniqueCoincidencePatternsExRF));
for i = 1:nCells
    uniqueCoincidencePatternsAsProbInRF(:,i) = probFiringInBinInRF(i) * uniqueCoincidencePatternsInRF(:,i) + ...
            (1 - probFiringInBinInRF(i)) * ~uniqueCoincidencePatternsInRF(:,i);
    uniqueCoincidencePatternsAsProbExRF(:,i) = probFiringInBinExRF(i) * uniqueCoincidencePatternsExRF(:,i) + ...
            (1 - probFiringInBinExRF(i)) * ~uniqueCoincidencePatternsExRF(:,i);
end
% joint prob occurrence assuming statistical independence
expectedProbUniquePatternsInRF = prod(uniqueCoincidencePatternsAsProbInRF, 2);
expectedCountUniquePatternsInRF = expectedProbUniquePatternsInRF * nBins * nTrialsInRF;
expectedProbUniquePatternsExRF = prod(uniqueCoincidencePatternsAsProbExRF, 2);
expectedCountUniquePatternsExRF = expectedProbUniquePatternsExRF * nBins * nTrialsExRF;

oneSidedPActualCountUniquePatternsInRFByPoisson = poisscdf(actualCountUniquePatternsInRF, expectedCountUniquePatternsInRF, 'upper');
oneSidedPActualCountUniquePatternsInRFByPoisson(oneSidedPActualCountUniquePatternsInRFByPoisson > 0.5) = 1 - oneSidedPActualCountUniquePatternsInRFByPoisson(oneSidedPActualCountUniquePatternsInRFByPoisson > 0.5);
pActualCountUniquePatternsInRFByPoisson = oneSidedPActualCountUniquePatternsInRFByPoisson * 2;

pAlpha = 0.05 / nUniquePatternsInRF; % correct for multiple comparisons.........
significantCoincidentPatternsInRFByPoisson = uniqueCoincidencePatternsInRF(pActualCountUniquePatternsInRFByPoisson < pAlpha,:)
significantCoincidentPatternsInRFByPoissonActualCount = actualCountUniquePatternsInRF(pActualCountUniquePatternsInRFByPoisson < pAlpha)
significantCoincidentPatternsInRFByPoissonExpectedCount = expectedCountUniquePatternsInRF(pActualCountUniquePatternsInRFByPoisson < pAlpha)

oneSidedPActualCountUniquePatternsExRFByPoisson = poisscdf(actualCountUniquePatternsExRF, expectedCountUniquePatternsExRF, 'upper');
oneSidedPActualCountUniquePatternsExRFByPoisson(oneSidedPActualCountUniquePatternsExRFByPoisson > 0.5) = 1 - oneSidedPActualCountUniquePatternsExRFByPoisson(oneSidedPActualCountUniquePatternsExRFByPoisson > 0.5);
pActualCountUniquePatternsExRFByPoisson = oneSidedPActualCountUniquePatternsExRFByPoisson * 2;

pAlpha = 0.05 / nUniquePatternsExRF; % correct for multiple comparisons.........
significantCoincidentPatternsExRFByPoisson = uniqueCoincidencePatternsExRF(pActualCountUniquePatternsExRFByPoisson < pAlpha,:)
significantCoincidentPatternsExRFByPoissonActualCount = actualCountUniquePatternsExRF(pActualCountUniquePatternsExRFByPoisson < pAlpha)
significantCoincidentPatternsExRFByPoissonExpectedCount = expectedCountUniquePatternsExRF(pActualCountUniquePatternsExRFByPoisson < pAlpha)


%% randomly permute each cell's spikes on each trial
% probFiringInBin does not change
% this preserves number of spikes per trial but does not preserve ISIs --
% could circularly rotate spikes for that 
nPermutations = 500;

permCountUniquePatternsInRF = nan(nUniquePatternsInRF, nPermutations);
permCountUniquePatternsExRF = nan(nUniquePatternsExRF, nPermutations);

for m = 1:nPermutations
    if mod(m, 50) == 0
        fprintf('Permutation %d/%d = %d%%...\n', m, nPermutations, m / nPermutations * 100);
    end
    % permute
    preArrayOnsetWindowedFiringInRFBinaryCellsPerm = nan(size(preArrayOnsetWindowedFiringInRFBinaryCells));
    preArrayOnsetWindowedFiringExRFBinaryCellsPerm = nan(size(preArrayOnsetWindowedFiringInRFBinaryCells));
    for i = 1:nCells
        for j = 1:nTrialsInRF
            preArrayOnsetWindowedFiringInRFBinaryCellsPerm(i,j,:) = preArrayOnsetWindowedFiringInRFBinaryCells(i,j,randperm(nBins));
        end
        for j = 1:nTrialsExRF
            preArrayOnsetWindowedFiringExRFBinaryCellsPerm(i,j,:) = preArrayOnsetWindowedFiringExRFBinaryCells(i,j,randperm(nBins));
        end
    end
    
    % for each bin n, find unique firing patterns k
    allCoincidencePatternsInRFPerm = false(nBins * nTrialsInRF, nCells);
    patternCount = 0;
    % TODO vectorize/use reshape
    for j = 1:nTrialsInRF
        for n = 1:nBins
            patternCount = patternCount + 1;
            allCoincidencePatternsInRFPerm(patternCount,:) = preArrayOnsetWindowedFiringInRFBinaryCellsPerm(:,j,n);
        end
    end
    
    allCoincidencePatternsExRFPerm = false(nBins * nTrialsExRF, nCells);
    patternCount = 0;
    % TODO vectorize/use reshape
    for j = 1:nTrialsExRF
        for n = 1:nBins
            patternCount = patternCount + 1;
            allCoincidencePatternsExRFPerm(patternCount,:) = preArrayOnsetWindowedFiringExRFBinaryCellsPerm(:,j,n);
        end
    end
    
    % count instances of the actual patterns in the permutation
    for k = 1:nUniquePatternsInRF
        [~,matchIndex] = ismember(allCoincidencePatternsInRFPerm, uniqueCoincidencePatternsInRF(k,:), 'rows');
        permCountUniquePatternsInRF(k,m) = sum(matchIndex);
    end
    for k = 1:nUniquePatternsExRF
        [~,matchIndex] = ismember(allCoincidencePatternsExRFPerm, uniqueCoincidencePatternsExRF(k,:), 'rows');
        permCountUniquePatternsExRF(k,m) = sum(matchIndex);
    end
end

%% compute p values based on the permutation distribution
pActualCountUniquePatternsInRFByPerm = nan(nUniquePatternsInRF, 1);
for k = 1:nUniquePatternsInRF
    oneSidedP = sum(permCountUniquePatternsInRF(k,:) > actualCountUniquePatternsInRF(k)) / nPermutations;
    if oneSidedP > 0.5
        pActualCountUniquePatternsInRFByPerm(k) = (1 - oneSidedP) * 2;
    else
        pActualCountUniquePatternsInRFByPerm(k) = oneSidedP * 2;
    end
end

pAlpha = 0.05 / nUniquePatternsInRF; % correct for multiple comparisons.........
significantCoincidentPatternsInRFByPerm = uniqueCoincidencePatternsInRF(pActualCountUniquePatternsInRFByPerm < pAlpha,:)
significantCoincidentPatternsInRFByPermActualCount = actualCountUniquePatternsInRF(pActualCountUniquePatternsInRFByPerm < pAlpha)

pActualCountUniquePatternsExRFByPerm = nan(nUniquePatternsExRF, 1);
for k = 1:nUniquePatternsExRF
    oneSidedP = sum(permCountUniquePatternsExRF(k,:) > actualCountUniquePatternsExRF(k)) / nPermutations;
    if oneSidedP > 0.5
        pActualCountUniquePatternsExRFByPerm(k) = (1 - oneSidedP) * 2;
    else
        pActualCountUniquePatternsExRFByPerm(k) = oneSidedP * 2;
    end
end

pAlpha = 0.05 / nUniquePatternsExRF; % correct for multiple comparisons.........
significantCoincidentPatternsExRFByPerm = uniqueCoincidencePatternsExRF(pActualCountUniquePatternsExRFByPerm < pAlpha,:)
significantCoincidentPatternsExRFByPermActualCount = actualCountUniquePatternsExRF(pActualCountUniquePatternsExRFByPerm < pAlpha)

%% TODO do this over a moving window 

