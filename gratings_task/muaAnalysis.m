function muaAnalysis(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad)
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

v = 9;
nLoc = 4;

%% input check
assert(numel(muaChannelsToLoad) == 1);

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channel to Load: %d\n', muaChannelsToLoad);
fprintf('Version: %d\n', v);

tic;
isLoadSpikes = 0;
isLoadMua = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 1;
spikeChannelsToLoad = [];
lfpChannelsToLoad = [];
spkcChannelsToLoad = [];
D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

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
minFiringRate = 1;

%%
nMUA = numel(D.allMUAStructs);
assert(nMUA == 1);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d MUAs...\n', nMUA);

for i = 1:nMUA
    muaStruct = D.allMUAStructs{i};
    unitName = muaStruct.name;
    spikeTimes = muaStruct.ts;
    nTrials = numel(UE.cueOnset);
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nMUA, round(i/nMUA*100));
    
    if firingRateOverall >= minFiringRateOverall
        tic;
        saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                processedDataDir, unitName, blockName, v);
        fprintf('\tWriting file %s...\n', saveFileName);
        
        computeEvokedSpiking(saveFileName, muaStruct, nLoc, UE);
        
        ES = load(saveFileName);
        if (any(ES.averageFiringRatesBySpdf.preEnterFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postEnterFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postEnterFixationLate.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.preCueBaseline.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayRelResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.preExitFixation.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.postExitFixation.byLoc >= minFiringRate))
            
            plotFileName = sprintf('%s/%s-%s-visual-v%d.png', processedDataDir, unitName, blockName, v);
        	fprintf('\tSaving figure to file %s...\n', plotFileName);
            
            quickSpdfAllEvents(saveFileName, blockName, ...
                    D, i, muaStruct, nLoc, nTrials, plotFileName);
            close;

            plotFileName = sprintf('%s/%s-%s-motor-v%d.png', processedDataDir, unitName, blockName, v);
        	fprintf('\tSaving figure to file %s...\n', plotFileName);
                
            quickSpdfAllMotorEvents(saveFileName, blockName, ...
                    D, i, muaStruct, nLoc, nTrials, plotFileName);
            close;
        else
            fprintf('\tTask-related firing rate < minimum task-related firing rate = %0.2f Hz - skipping.\n', ...
                    minFiringRate);
        end
        fprintf('\tTime elapsed: %0.2f s.\n', toc);
    else
        fprintf('\tOverall firing rate = %0.2f Hz < minimum firing rate = %0.2f Hz in these blocks - skipping.\n', ...
                firingRateOverall, minFiringRateOverall);
    end
end

end










