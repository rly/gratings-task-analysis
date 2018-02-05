
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
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
muaDataDirRoot = processedDataRootDir;

v = 9;
sessionIndAll = 1%[1 2 3 4 5 6 7 8];
for k = 1:numel(sessionIndAll)
    clearvars -except processedDataRootDir dataDirRoot sessionIndAll k v % cheap hack
    sessionInd = sessionIndAll(k);

nLoc = 4;

%% load recording information
recordingInfo = readRecordingInfo();
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
fprintf('Session index: %d\n', sessionInd);

tic;
isLoadSpikes = 1;
isLoadMua = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 1;
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
% for l = sessionRange
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

stop

%% process every cell in workspace
nUnits = numel(D.allSpikeStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d units...\n', nUnits);

for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    nTrials = numel(UE.cueOnset);
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    
    if firingRateOverall >= minFiringRateOverall
        tic;
        saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                processedDataDir, unitName, blockName, v);
        fprintf('\tWriting file %s...\n', saveFileName);
        
        computeEvokedSpiking(saveFileName, spikeStruct, nLoc, UE);
        
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
            
%             quickSpdfAllEvents(saveFileName, blockName, ...
%                     D, i, spikeStruct, nLoc, nTrials, plotFileName);
%             close;

            plotFileName = sprintf('%s/%s-%s-motor-v%d.png', processedDataDir, unitName, blockName, v);
        	fprintf('\tSaving figure to file %s...\n', plotFileName);
                
%             quickSpdfAllMotorEvents(saveFileName, blockName, ...
%                     D, i, spikeStruct, nLoc, nTrials, plotFileName);
%             close;
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
% stop

%% summarize across population of cells
% use only cells with time-locked response > 1 Hz

nUnits = numel(D.allSpikeStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d units...\n', nUnits);

minFiringRate = 1; % Hz
nIncludedUnits = 0;
isSignificantStats = false(nUnits, 10); % 5 periods > baseline, 5 periods info rate
isCell = false(nUnits, 1);
isBroadSpiking = false(nUnits, 1);
isNarrowSpiking = false(nUnits, 1);
statAlpha = 0.05/3; % reduce from 0.05 to account for multiple comparisons... 
% but also need to increase number of randomizations

for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    nTrials = numel(UE.cueOnset);
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;

    if firingRateOverall >= minFiringRateOverall
        saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
                processedDataDir, unitName, blockName, v);
            
        ES = load(saveFileName);
        
        if (any(ES.averageFiringRatesBySpdf.preCueBaseline.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.cueTargetDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.arrayHoldResponse.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimDelay.byLoc >= minFiringRate) || ...
                any(ES.averageFiringRatesBySpdf.targetDimResponse.byLoc >= minFiringRate))
            
            isCell(i) = strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking');
            if isCell(i)
                isBroadSpiking(i) = strcmp(spikeStruct.physClass, 'Broad-Spiking');
                isNarrowSpiking(i) = strcmp(spikeStruct.physClass, 'Narrow-Spiking');
            end
            
            isSignificantStats(i,:) = [...
                    min([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.cueResponseVsBaselineRankSumTestStatsByLoc.p]) ...
                    min([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.cueTargetDelayVsBaselineRankSumTestStatsByLoc.p]) ...
                    min([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.arrayHoldResponseVsBaselineRankSumTestStatsByLoc.p]) ...
                    min([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.targetDimDelayVsBaselineRankSumTestStatsByLoc.p]) ...
                    min([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) < statAlpha / numel([ES.targetDimResponseVsBaselineRankSumTestStatsByLoc.p]) ...
                    ES.cueResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.cueTargetDelayInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.arrayHoldResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.targetDimDelayInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha ...
                    ES.targetDimResponseInfoRateStruct.infoRatePValueByShuffleSpdf < statAlpha];
            nIncludedUnits = nIncludedUnits + 1;
        end
    end
end

nCells = sum(isBroadSpiking | isNarrowSpiking);

isSignificantAnyTaskMod = isCell & any(isSignificantStats(:,1:5), 2);
isSignificantAnySpatialSelectivity = isCell & any(isSignificantStats(:,6:10), 2);
isSignificantEvokedSelectivity = isCell & any(isSignificantStats(:,[6 8 10]), 2);
isSignificantDelaySelectivity = isCell & any(isSignificantStats(:,[7 9]), 2);
isSigSelectEvokedNotDelay = isCell & any(isSignificantStats(:,[6 8 10]), 2) & any(isSignificantStats(:,[7 9]), 2);
isSigSelectOnlyEvoked = isCell & any(isSignificantStats(:,[6 8 10]), 2) & ~any(isSignificantStats(:,[7 9]), 2);
isSigSelectOnlyDelay = isCell & ~any(isSignificantStats(:,[6 8 10]), 2) & any(isSignificantStats(:,[7 9]), 2);


fprintf('\n%s: \n', sessionName);
fprintf('%d/%d = %d%% units have spike rate during some task period >= %0.1f Hz.\n', ...
        nIncludedUnits, nUnits, ...
        round(nIncludedUnits/nUnits * 100), minFiringRate);
fprintf('%d/%d = %d%% units were classified as cells (BS or NS and rate >= %0.1f Hz).\n\n', ...
        nCells, nUnits, ...
        round(nCells/nUnits * 100), minFiringRate);
fprintf('%d/%d = %d%% cells were classified as broad-spiking and %d/%d = %d%% as narrow-spiking.\n', ...
        sum(isBroadSpiking(isCell)), nCells, ...
        round(sum(isBroadSpiking(isCell))/nCells * 100), ...
        sum(isNarrowSpiking(isCell)), nCells, ...
        round(sum(isNarrowSpiking(isCell))/nCells * 100));
fprintf('%d/%d = %d%% cells show significant task modulation compared to baseline.\n', ...
        sum(isSignificantAnyTaskMod), nCells, ...
        round(sum(isSignificantAnyTaskMod)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity during some task period.\n', ...
        sum(isSignificantAnySpatialSelectivity), nCells, ...
        round(sum(isSignificantAnySpatialSelectivity)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity in a task period involving a visual change.\n', ...
        sum(isSignificantEvokedSelectivity), nCells, ...
        round(sum(isSignificantEvokedSelectivity)/nCells * 100));
fprintf('%d/%d = %d%% cells show significant spatial selectivity in a task period involving sustained attention.\n\n', ...
        sum(isSignificantDelaySelectivity), nCells, ...
        round(sum(isSignificantDelaySelectivity)/nCells * 100));
fprintf('Of the %d cells that show significant spatial selectivity during some task period,\n', ...
        sum(isSignificantAnySpatialSelectivity));
fprintf('\t%d (%d%%) show significant spatial selectivity after a visual change AND during attention,\n', ...
        sum(isSigSelectEvokedNotDelay), ...
        round(sum(isSigSelectEvokedNotDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('\t%d (%d%%) show significant spatial selectivity ONLY after a visual change, and\n', ...
        sum(isSigSelectOnlyEvoked), ...
        round(sum(isSigSelectOnlyEvoked)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('\t%d (%d%%) show significant spatial selectivity ONLY during attention.\n\n', ...
        sum(isSigSelectOnlyDelay), ...
        round(sum(isSigSelectOnlyDelay)/sum(isSignificantAnySpatialSelectivity) * 100));
fprintf('Of the %d cells that show significant task modulation compared to baseline, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantAnyTaskMod), ...
        sum(isBroadSpiking(isSignificantAnyTaskMod)), ...
        round(sum(isBroadSpiking(isSignificantAnyTaskMod))/sum(isSignificantAnyTaskMod) * 100), ...
        sum(isNarrowSpiking(isSignificantAnyTaskMod)), ...
        round(sum(isNarrowSpiking(isSignificantAnyTaskMod))/sum(isSignificantAnyTaskMod) * 100));
fprintf('Of the %d cells that show significant spatial selectivity after a visual change, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantEvokedSelectivity), ...
        sum(isBroadSpiking(isSignificantEvokedSelectivity)), ...
        round(sum(isBroadSpiking(isSignificantEvokedSelectivity))/sum(isSignificantEvokedSelectivity) * 100), ...
        sum(isNarrowSpiking(isSignificantEvokedSelectivity)), ...
        round(sum(isNarrowSpiking(isSignificantEvokedSelectivity))/sum(isSignificantEvokedSelectivity) * 100));
fprintf('Of the %d cells that show significant spatial selectivity during attention, \n\t%d (%d%%) are broad-spiking and %d (%d%%) are narrow-spiking.\n', ...
        sum(isSignificantDelaySelectivity), ...
        sum(isBroadSpiking(isSignificantDelaySelectivity)), ...
        round(sum(isBroadSpiking(isSignificantDelaySelectivity))/sum(isSignificantDelaySelectivity) * 100), ...
        sum(isNarrowSpiking(isSignificantDelaySelectivity)), ...
        round(sum(isNarrowSpiking(isSignificantDelaySelectivity))/sum(isSignificantDelaySelectivity) * 100));

fprintf('\n');
fprintf('Cells with task modulation (%d): \n\t', sum(isSignificantAnyTaskMod));
cellfun(@(s) fprintf('%s, ', s), arrayfun(@(x) D.allSpikeStructs{x}.unitIDChar, find(isSignificantAnyTaskMod), 'UniformOutput', false));
fprintf('\n');
fprintf('Cells with evoked selectivity (%d): \n\t', sum(isSignificantEvokedSelectivity));
cellfun(@(s) fprintf('%s, ', s), arrayfun(@(x) D.allSpikeStructs{x}.unitIDChar, find(isSignificantEvokedSelectivity), 'UniformOutput', false));
fprintf('\n');
fprintf('Cells with delay selectivity (%d): \n\t', sum(isSignificantDelaySelectivity));
cellfun(@(s) fprintf('%s, ', s), arrayfun(@(x) D.allSpikeStructs{x}.unitIDChar, find(isSignificantDelaySelectivity), 'UniformOutput', false));
fprintf('\n');

end
stop

%%
nMUA = numel(D.allMUAStructs);
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
% stop














%% compute attention indices

nUnits = numel(D.allSpikeStructs);
kernelSigma = 0.01;

periCueOnsetWindow = [0.4 0.5];
preCueOnsetSpdfWindowOffset = [-0.3 0];
postCueOnsetSpdfWindowOffset = [0.025 0.225];

periArrayOnsetWindow = [0.7 0.5];
preArrayOnsetSpdfWindowOffset = [-0.4 0];
postArrayOnsetSpdfWindowOffset = [0 0.3];

periTargetDimWindow = [0.7 0.5];
preTargetDimSpdfWindowOffset = [-0.5 0];
postTargetDimSpdfWindowOffset = [-0.5 0.3];

preCueOnsetT = computeTForSpdf(periCueOnsetWindow(1), preCueOnsetSpdfWindowOffset, kernelSigma);
postCueOnsetT = computeTForSpdf(periCueOnsetWindow(1), postCueOnsetSpdfWindowOffset, kernelSigma);
preArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), preArrayOnsetSpdfWindowOffset, kernelSigma);
postArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), postArrayOnsetSpdfWindowOffset, kernelSigma);
preTargetDimT = computeTForSpdf(periTargetDimWindow(1), preTargetDimSpdfWindowOffset, kernelSigma);
postTargetDimT = computeTForSpdf(periTargetDimWindow(1), postTargetDimSpdfWindowOffset, kernelSigma);

wfAll = nan(nUnits, numel(D.allSpikeStructs{1}.meanWf));
isNSCell = nan(nUnits, 1);

preCueOnsetSpdfByLoc = nan(nUnits, nLoc, numel(preCueOnsetT));
postCueOnsetSpdfByLoc = nan(nUnits, nLoc, numel(postCueOnsetT));
cueResponseByLoc = nan(nUnits, nLoc);
inRFLocAll = nan(nUnits, 1);
exRFLocAll = nan(nUnits, 1);

preCueOnsetSpdfInRFAll = nan(nUnits, numel(preCueOnsetT));
preCueOnsetSpdfExRFAll = nan(nUnits, numel(preCueOnsetT));

preArrayOnsetSpdfInRFAll = nan(nUnits, numel(preArrayOnsetT));
preArrayOnsetSpdfExRFAll = nan(nUnits, numel(preArrayOnsetT));
preArrayOnsetSpdfAIAll = nan(nUnits, numel(preArrayOnsetT));

postArrayOnsetSpdfInRFAll = nan(nUnits, numel(postArrayOnsetT));
postArrayOnsetSpdfExRFAll = nan(nUnits, numel(postArrayOnsetT));
postArrayOnsetSpdfAIAll = nan(nUnits, numel(postArrayOnsetT));

preTargetDimSpdfInRFAll = nan(nUnits, numel(preTargetDimT));
preTargetDimSpdfExRFAll = nan(nUnits, numel(preTargetDimT));
preTargetDimSpdfAIAll = nan(nUnits, numel(preTargetDimT));

postTargetDimSpdfInRFAll = nan(nUnits, numel(postTargetDimT));
postTargetDimSpdfExRFAll = nan(nUnits, numel(postTargetDimT));
postTargetDimSpdfAIAll = nan(nUnits, numel(postTargetDimT));

cueTargetDelayAIAll = nan(nUnits, 1);
arrayResponseAIAll = nan(nUnits, 1);
targetDimDelayAIAll = nan(nUnits, 1);

minFiringRate = 1; % Hz
minAbsAI = 0.02; % TODO assess based on significance
unitCueTargetDelayCategory = ones(nUnits, 1);
unitTargetDimDelayCategory = ones(nUnits, 1);

fprintf('-------------------------------------------------------------\n');
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;

%     fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
%             nUnits, round(i/nUnits*100));
    if ~isempty(spikeTimes)
        % process only classified BS and NS units
        if strcmp(spikeStruct.physClass, 'Broad-Spiking') || strcmp(spikeStruct.physClass, 'Narrow-Spiking')
            wfAll(i,:) = spikeStruct.meanWf / max(spikeStruct.meanWf);
            isNSCell(i) = strcmp(spikeStruct.physClass, 'Narrow-Spiking');
            
            % compute RF
            for j = 1:nLoc
                preCueOnsetSpikeTimesAtLocJ = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{j}, periCueOnsetWindow);
                preCueOnsetSpdfByLoc(i,j,:) = edpsth_notranspose(preCueOnsetSpikeTimesAtLocJ, kernelSigma, 'n', [], 0, preCueOnsetT);
                postCueOnsetSpikeTimesAtLocJ = createnonemptydatamatpt(spikeVar, cueOnsetByLoc{j}, periCueOnsetWindow);
                postCueOnsetSpdfByLoc(i,j,:) = edpsth_notranspose(postCueOnsetSpikeTimesAtLocJ, kernelSigma, 'n', [], 0, postCueOnsetT);
                
                % mean baseline-corrected post-cue response
                cueResponseByLoc(i,j) = mean(postCueOnsetSpdfByLoc(i,j,:)) - mean(preCueOnsetSpdfByLoc(i,j,:));
                
                % TODO test significance using time-SD and bootstrap
            end
            
            [~,inRFLocAll(i)] = max(cueResponseByLoc(i,:));
            assert(nLoc == 4); % next line based on nLoc == 4
            exRFLocAll(i) = mod(inRFLocAll(i) + 1, 4) + 1; % opposite location
            
            cueResponseInRF = mean(cueResponseByLoc(i,inRFLocAll(i)));
            if cueResponseInRF < minFiringRate
                continue;
            end
            
            preCueOnsetSpdfInRFAll(i,:) = squeeze(preCueOnsetSpdfByLoc(i,inRFLocAll(i),:));
            preCueOnsetSpdfExRFAll(i,:) = squeeze(preCueOnsetSpdfByLoc(i,exRFLocAll(i),:));
            
            preArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLocAll(i)}, periArrayOnsetWindow);
            preArrayOnsetSpdfInRFAll(i,:) = edpsth_notranspose(preArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            preArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLocAll(i)}, periArrayOnsetWindow);
            preArrayOnsetSpdfExRFAll(i,:) = edpsth_notranspose(preArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            
            postArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLocAll(i)}, periArrayOnsetWindow);
            postArrayOnsetSpdfInRFAll(i,:) = edpsth_notranspose(postArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, postArrayOnsetT);
            postArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLocAll(i)}, periArrayOnsetWindow);
            postArrayOnsetSpdfExRFAll(i,:) = edpsth_notranspose(postArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, postArrayOnsetT);

            preTargetDimSpikeTimesInRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{inRFLocAll(i)}, periTargetDimWindow);
            preTargetDimSpdfInRFAll(i,:) = edpsth_notranspose(preTargetDimSpikeTimesInRF, kernelSigma, 'n', [], 0, preTargetDimT);
            preTargetDimSpikeTimesExRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{exRFLocAll(i)}, periTargetDimWindow);
            preTargetDimSpdfExRFAll(i,:) = edpsth_notranspose(preTargetDimSpikeTimesExRF, kernelSigma, 'n', [], 0, preTargetDimT);
            
            postTargetDimSpikeTimesInRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{inRFLocAll(i)}, periTargetDimWindow);
            postTargetDimSpdfInRFAll(i,:) = edpsth_notranspose(postTargetDimSpikeTimesInRF, kernelSigma, 'n', [], 0, postTargetDimT);
            postTargetDimSpikeTimesExRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{exRFLocAll(i)}, periTargetDimWindow);
            postTargetDimSpdfExRFAll(i,:) = edpsth_notranspose(postTargetDimSpikeTimesExRF, kernelSigma, 'n', [], 0, postTargetDimT);

            preArrayOnsetSpdfAIAll(i,:) = (preArrayOnsetSpdfInRFAll(i,:) - preArrayOnsetSpdfExRFAll(i,:)) ./ ...
                            (preArrayOnsetSpdfInRFAll(i,:) + preArrayOnsetSpdfExRFAll(i,:));
            postArrayOnsetSpdfAIAll(i,:) = (postArrayOnsetSpdfInRFAll(i,:) - postArrayOnsetSpdfExRFAll(i,:)) ./ ...
                            (postArrayOnsetSpdfInRFAll(i,:) + postArrayOnsetSpdfExRFAll(i,:));
            preTargetDimSpdfAIAll(i,:) = (preTargetDimSpdfInRFAll(i,:) - preTargetDimSpdfExRFAll(i,:)) ./ ...
                            (preTargetDimSpdfInRFAll(i,:) + preTargetDimSpdfExRFAll(i,:));
            postTargetDimSpdfAIAll(i,:) = (postTargetDimSpdfInRFAll(i,:) - postTargetDimSpdfExRFAll(i,:)) ./ ...
                            (postTargetDimSpdfInRFAll(i,:) + postTargetDimSpdfExRFAll(i,:));
            
            % attention index, using mean responses over time 
            % -- not equivalent to mean AI over time
            % TODO: do baseline correction here??
            cueTargetDelayAIAll(i) = (mean(preArrayOnsetSpdfInRFAll(i,:) - preArrayOnsetSpdfExRFAll(i,:))) / ...
                    (mean(preArrayOnsetSpdfInRFAll(i,:) + preArrayOnsetSpdfExRFAll(i,:)));
            arrayResponseAIAll(i) = (mean(postArrayOnsetSpdfInRFAll(i,:) - postArrayOnsetSpdfExRFAll(i,:))) / ...
                    (mean(postArrayOnsetSpdfInRFAll(i,:) + postArrayOnsetSpdfExRFAll(i,:)));
            targetDimDelayAIAll(i) = (mean(preTargetDimSpdfInRFAll(i,:) - preTargetDimSpdfExRFAll(i,:))) / ...
                    (mean(preTargetDimSpdfInRFAll(i,:) + preTargetDimSpdfExRFAll(i,:)));
            stop
            if any(isnan(cueTargetDelayAIAll(i)))
                continue;
            end

            if isNSCell(i)
                fprintf('Unit %d (%s) - NS -- FR(cueInRF): %0.2f Hz, AI(cue-target): %0.2f, AI(target-dim): %0.2f\n', ...
                        i, unitName, cueResponseInRF, cueTargetDelayAIAll(i), targetDimDelayAIAll(i));
                
                if cueTargetDelayAIAll(i) > minAbsAI
                    unitCueTargetDelayCategory(i) = 7;
                elseif cueTargetDelayAIAll(i) < -minAbsAI
                    unitCueTargetDelayCategory(i) = 5;
                else
                    unitCueTargetDelayCategory(i) = 6;
                end
                if targetDimDelayAIAll(i) > minAbsAI
                    unitTargetDimDelayCategory(i) = 7;
                elseif targetDimDelayAIAll(i) < -minAbsAI
                    unitTargetDimDelayCategory(i) = 5;
                else
                    unitTargetDimDelayCategory(i) = 6;
                end
                
            else
                fprintf('Unit %d (%s) - BS -- FR(cueInRF): %0.2f Hz, AI(cue-target): %0.2f, AI(target-dim): %0.2f\n', ...
                        i, unitName, cueResponseInRF, cueTargetDelayAIAll(i), targetDimDelayAIAll(i));
                
                if cueTargetDelayAIAll(i) > minAbsAI
                    unitCueTargetDelayCategory(i) = 4;
                elseif cueTargetDelayAIAll(i) < -minAbsAI
                    unitCueTargetDelayCategory(i) = 2;
                else
                    unitCueTargetDelayCategory(i) = 3;
                end
                if targetDimDelayAIAll(i) > minAbsAI
                    unitTargetDimDelayCategory(i) = 4;
                elseif targetDimDelayAIAll(i) < -minAbsAI
                    unitTargetDimDelayCategory(i) = 2;
                else
                    unitTargetDimDelayCategory(i) = 3;
                end
            end
        end
    end
end

%% draw channel map with NS/BS and attentional modulation
allChannelIDs = cellfun(@(x) x.channelID, D.allSpikeStructs);    
cueTargetDelayHistCounts = nan(5, 64);
targetDimDelayHistCounts = nan(5, 64);
for i = 1:7
    cueTargetDelayHistCounts(i,:) = histcounts(allChannelIDs(unitCueTargetDelayCategory == i), 1:65);
    targetDimDelayHistCounts(i,:) = histcounts(allChannelIDs(unitTargetDimDelayCategory == i), 1:65);
end

figure_tr_inch(12, 10);
subaxis(1, 3, 1);
bh = barh(cueTargetDelayHistCounts', 'stacked');
bh(1).FaceColor = 0.7*ones(3,1);
bh(2).FaceColor = [255 200 100]/255;
bh(3).FaceColor = [200 150 50]/255;
bh(4).FaceColor = [150 100 0]/255;
bh(5).FaceColor = [255 50 255]/255;
bh(6).FaceColor = [175 25 200]/255;
bh(7).FaceColor = [100 0 150]/255;

set(gca, 'YDir', 'reverse');
ylim([0 65]);
origXLim = xlim();
hold on;
plot(origXLim, [32.5 32.5], 'Color', 0.3*ones(3, 1));
xlim(origXLim);
xlabel('Number of Units');
ylabel('Channel Number (1 = top probe 1, 33 = top probe 2)');
title('Cue-Target Delay AI');
set(gca, 'FontSize', 12);

% TODO match cells from ax1 with ax2
subaxis(1, 3, 2);
bh = barh(targetDimDelayHistCounts', 'stacked');
bh(1).FaceColor = 0.7*ones(3,1);
bh(2).FaceColor = [255 200 100]/255;
bh(3).FaceColor = [200 150 50]/255;
bh(4).FaceColor = [150 100 0]/255;
bh(5).FaceColor = [255 50 255]/255;
bh(6).FaceColor = [175 25 200]/255;
bh(7).FaceColor = [100 0 150]/255;

set(gca, 'YDir', 'reverse');
ylim([0 65]);
origXLim = xlim();
hold on;
plot(origXLim, [32.5 32.5], 'Color', 0.3*ones(3, 1));
xlim(origXLim);
xlabel('Number of Units');
% ylabel('Channel Number (1 = top probe 1, 33 = top probe 2)');
title('Target-Dim Delay AI');
set(gca, 'FontSize', 12);

legendH = legend({'Axon? / Low FR', ...
        sprintf('BS AI<%0.2f', -minAbsAI), ...
        'BS', ...
        sprintf('BS AI>%0.2f', minAbsAI), ...
        sprintf('NS AI<%0.2f', -minAbsAI), ...
        'NS', ...
        sprintf('NS AI>%0.2f', minAbsAI)}, ...
        'Location', 'NorthEastOutside');
origLegendPosition = get(legendH, 'Position');

subaxis(1, 3, 3);
set(gca, 'Visible', 'off');
ax3Position = get(gca, 'Position');
set(legendH, 'Position', [ax3Position(1) origLegendPosition(2:4)]);

plotFileName = sprintf('%s/delayHistAIByCellTypeChannel.png', processedDataDir);
export_fig(plotFileName, '-nocrop');
