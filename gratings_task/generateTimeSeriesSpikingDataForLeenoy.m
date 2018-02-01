
% 325ms fixation before pre-cue marker
% 25-125ms fixation between pre-cue marker and cue onset
% 100ms cue onset to cue offset
% 500-800ms cue offset to array onset
% 850-1050ms array onset to target dim for long hold trials
% 650-850ms array onset to target dim for long hold trials
% hold shape response window 280-730ms
% release shape response window 280-730ms
% min 280ms before saccade allowed

% evt5 = cue onset
% evt1,2,3,4 = cue offset
% evt6 = array onset and also release shape resp win start
% evt7 = target dim
% evt8 = juice

clear;

sessionInd = 3;

nLoc = 4;

logRoot = 'gratings_task_eye_tracking_';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
D = loadPL2(pl2FilePath, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 

processedDataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/%s', sessionName);
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
% for l = sessionRange
fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
logDir = sprintf('%s/%s', dataDir, sessionName(2:end));
load('params.mat');

% process events and sort them into different conditions
usefulEvents = getUsefulEvents2(logDir, gratingsTask3DLogIndices, 4, D);

% unload these variables to workspace
% (cueOnset, cueOnsetByLoc, ...
%         arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
%         arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%         arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%         targetDim, targetDimByLoc, ...
%         targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
%         nTrialShortHold, nTrialLongHold, rt);
struct2var(usefulEvents);

%%
nUnits = numel(D.allSpikeStructs);

binFs = 1000;

periCueOnsetWindow = [1 1];
periArrayOnsetWindow = [1 1];
periTargetDimWindow = [1 1];

nTimePeriCueOnset = round(sum(periCueOnsetWindow) * binFs);
nTimePeriArrayOnset = round(sum(periArrayOnsetWindow) * binFs);
nTimePeriTargetDim = round(sum(periTargetDimWindow) * binFs);

% for simplicity, hard-code inRFLoc and exRFLoc
inRFLoc = 3;
exRFLoc = 1;

periCueOnsetSpikeTimesInRFBinary = false(nUnits, numel(cueOnsetByLoc{inRFLoc}), nTimePeriCueOnset);
periCueOnsetSpikeTimesExRFBinary = false(nUnits, numel(cueOnsetByLoc{exRFLoc}), nTimePeriCueOnset);
periArrayOnsetSpikeTimesInRFBinary = false(nUnits, numel(arrayOnsetByLoc{inRFLoc}), nTimePeriArrayOnset);
periArrayOnsetSpikeTimesExRFBinary = false(nUnits, numel(arrayOnsetByLoc{exRFLoc}), nTimePeriArrayOnset);
periTargetDimSpikeTimesInRFBinary = false(nUnits, numel(targetDimByLoc{inRFLoc}), nTimePeriTargetDim);
periTargetDimSpikeTimesExRFBinary = false(nUnits, numel(targetDimByLoc{exRFLoc}), nTimePeriTargetDim);

% kernelSigma = 0.01;
% preCueOnsetSpdfWindowOffset = [-0.3 0];
% postCueOnsetSpdfWindowOffset = [0.025 0.225];
% preCueOnsetT = computeTForSpdf(periCueOnsetWindow(1), preCueOnsetSpdfWindowOffset, kernelSigma);
% postCueOnsetT = computeTForSpdf(periCueOnsetWindow(1), postCueOnsetSpdfWindowOffset, kernelSigma);
% inRFLocAll = nan(nUnits, 1);
% exRFLocAll = nan(nUnits, 1);

fprintf('-------------------------------------------------------------\n');
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;

    if ~isempty(spikeTimes)

%         cueResponseByLoc = nan(nLoc, 1);
%         % compute RF
%         for j = 1:nLoc
%             preCueOnsetSpikeTimesAtLocJ = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{j}, periCueOnsetWindow);
%             preCueOnsetSpdfByLoc = edpsth_notranspose(preCueOnsetSpikeTimesAtLocJ, kernelSigma, 'n', [], 0, preCueOnsetT);
%             postCueOnsetSpikeTimesAtLocJ = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{j}, periCueOnsetWindow);
%             postCueOnsetSpdfByLoc = edpsth_notranspose(postCueOnsetSpikeTimesAtLocJ, kernelSigma, 'n', [], 0, postCueOnsetT);
% 
%             % mean baseline-corrected post-cue response
%             cueResponseByLoc(j) = mean(postCueOnsetSpdfByLoc) - mean(preCueOnsetSpdfByLoc);
% 
%             % TODO test significance using time-SD and bootstrap
%         end
% 
%         [~,inRFLocAll(i)] = max(cueResponseByLoc);
%         assert(nLoc == 4); % next line based on nLoc == 4
%         exRFLocAll(i) = mod(inRFLocAll(i) + 1, 4) + 1; % opposite location

        periCueOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{inRFLoc}, periCueOnsetWindow);
        periCueOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{exRFLoc}, periCueOnsetWindow);
        
        periArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, periArrayOnsetWindow);
        periArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, periArrayOnsetWindow);
        
        periTargetDimSpikeTimesInRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{inRFLoc}, periTargetDimWindow);
        periTargetDimSpikeTimesExRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{exRFLoc}, periTargetDimWindow);
        
        for j = 1:numel(periCueOnsetSpikeTimesInRF)
        	periCueOnsetSpikeTimesInRFBinary(i,j,(ceil(binFs*periCueOnsetSpikeTimesInRF(j).times))) = true;
        end
        
        for j = 1:numel(periCueOnsetSpikeTimesExRF)
        	periCueOnsetSpikeTimesExRFBinary(i,j,(ceil(binFs*periCueOnsetSpikeTimesExRF(j).times))) = true;
        end
        
        for j = 1:numel(periArrayOnsetSpikeTimesInRF)
        	periArrayOnsetSpikeTimesInRFBinary(i,j,(ceil(binFs*periArrayOnsetSpikeTimesInRF(j).times))) = true;
        end
        
        for j = 1:numel(periArrayOnsetSpikeTimesExRF)
        	periArrayOnsetSpikeTimesExRFBinary(i,j,(ceil(binFs*periArrayOnsetSpikeTimesExRF(j).times))) = true;
        end
        
        for j = 1:numel(periTargetDimSpikeTimesInRF)
        	periTargetDimSpikeTimesInRFBinary(i,j,(ceil(binFs*periTargetDimSpikeTimesInRF(j).times))) = true;
        end
        
        for j = 1:numel(periTargetDimSpikeTimesExRF)
        	periTargetDimSpikeTimesExRFBinary(i,j,(ceil(binFs*periTargetDimSpikeTimesExRF(j).times))) = true;
        end
    end
end

%% save

saveFileName = sprintf('%s/%s-gratingsTaskBinaryTimePeriods.mat', processedDataDir, blockName);
save(saveFileName, 'periCueOnsetSpikeTimesInRFBinary', ...
        'periCueOnsetSpikeTimesExRFBinary', ...
        'periArrayOnsetSpikeTimesInRFBinary', ...
        'periArrayOnsetSpikeTimesExRFBinary', ...
        'periTargetDimSpikeTimesInRFBinary', ...
        'periTargetDimSpikeTimesExRFBinary');

