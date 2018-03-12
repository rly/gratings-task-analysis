function muaAnalysis(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad, isZeroDistractors)
% MUA gratings task analysis, one channel

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

v = 11;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task MUA Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('MUA Channel to Load: %d\n', muaChannelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

nLoc = 4;

%% input check
assert(numel(muaChannelsToLoad) == 1);

%% load recording information

if isZeroDistractors
    scriptName = 'MUA_GRATINGS_0D';
    taskName = 'GRATINGS_0D';
else
    scriptName = 'MUA_GRATINGS';
    taskName = 'GRATINGS';
end

[R, D, processedDataDir, blockName] = loadRecordingData(processedDataRootDir, ...
        dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad, ...
        taskName, scriptName, 1, 0);
sessionName = R.sessionName;

fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

% process events and sort them into different conditions
UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D);

totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
minFiringRateOverall = 0.2; % Hz
minFiringRate = 1; % Hz

nUnits = numel(D.allMUAStructs);
assert(nUnits == 1);
i = 1; % just one unit at a time in this script

%% compute evoked spiking and make SPDF plots for visual and motor events
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d MUAs...\n', nUnits);

muaStruct = D.allMUAStructs{i};
unitName = muaStruct.name;
spikeTimes = muaStruct.ts;
nTrials = numel(UE.cueOnset);
firingRateOverall = numel(spikeTimes) / totalTimeOverall;
fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
        nUnits, round(i/nUnits*100));

saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
    	processedDataDir, unitName, blockName, v);
if firingRateOverall >= minFiringRateOverall
    saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            processedDataDir, unitName, blockName, v);
    fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);

    computeEvokedSpiking(saveFileName, muaStruct, nLoc, UE);

    ES = load(saveFileName); % load once and pass around
    if isFiringRateGreaterThanMin(ES, minFiringRate)
        plotFileName = sprintf('%s/%s-%s-visual-v%d.png', processedDataDir, unitName, blockName, v);
        fprintf('\tSaving figure to file %s...\n', plotFileName);

        quickSpdfAllVisualEvents(ES, blockName, ...
                D, i, muaStruct, nLoc, nTrials, isZeroDistractors, plotFileName);
        close;

        plotFileName = sprintf('%s/%s-%s-motor-v%d.png', processedDataDir, unitName, blockName, v);
        fprintf('\tSaving figure to file %s...\n', plotFileName);

        quickSpdfAllMotorEvents(ES, blockName, ...
                D, i, muaStruct, nLoc, nTrials, isZeroDistractors, plotFileName);
        close;
    else
        fprintf('\tTask-related firing rate < minimum task-related firing rate = %0.2f Hz - skipping.\n', ...
                minFiringRate);
    end
else
    fprintf('\tOverall firing rate = %0.2f Hz < minimum firing rate = %0.2f Hz in these blocks - skipping.\n', ...
            firingRateOverall, minFiringRateOverall);
    if exists(saveFileName, 'file')
        fprintf('\t%s exists... deleting.\n', saveFileName);
        delete(saveFileName);
    end 
end

fprintf('Time elapsed: %0.2f s.\n', toc);

