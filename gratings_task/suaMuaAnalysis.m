function suaMuaAnalysis(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ...
        isZeroDistractors, numRandomizations, isLoadSortedSua, isLoadMua)
% SUA & MUA gratings task analysis, one channel

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

v = 15;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task SUA/MUA Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('Channels to Load: %d\n', channelsToLoad);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('SUA/MUA data root dir: %s\n', suaMuaDataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

nLoc = 4;

%% input check
assert(numel(channelsToLoad) == 1);

%% load recording information
if isZeroDistractors
    scriptName = 'SUA_MUA_GRATINGS_0D';
    taskName = 'GRATINGS_0D';
else
    scriptName = 'SUA_MUA_GRATINGS';
    taskName = 'GRATINGS';
end

isLoadLfp = 0;
isLoadMetaDataOnly = 0;
minSuaSepQuality = 3;
paramsStruct = var2struct(processedDataRootDir, ...
        dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ...
        taskName, scriptName, isLoadSortedSua, isLoadMua, isLoadLfp, isLoadMetaDataOnly, minSuaSepQuality);

[R, D, processedDataDir, blockName] = loadRecordingData2(paramsStruct);
sessionName = R.sessionName;

fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

% process events and sort them into different conditions
UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D, blockName);

totalTimeOverall = sum(D.blockStopTimes(R.blockIndices) - D.blockStartTimes(R.blockIndices));
minFiringRateOverall = 0.2; % Hz

nUnits = numel(D.allUnitStructs);
fprintf('-------------------------------------------------------------\n');
fprintf('Processing %d Units...\n', nUnits);

%% compute evoked spiking and make SPDF plots for visual and motor events
for i = 1:nUnits
    unitStruct = D.allUnitStructs{i};
    unitName = unitStruct.name;
    spikeTimes = unitStruct.ts;
    firingRateOverall = numel(spikeTimes) / totalTimeOverall;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));

    saveFileName = sprintf('%s/%s-%s-evokedSpiking-v%d.mat', ...
            processedDataDir, unitName, blockName, v);
    if firingRateOverall >= minFiringRateOverall
        fprintf('\tOverall firing rate = %0.2f Hz > minimum firing rate = %0.2f Hz in these blocks.\n', ...
                firingRateOverall, minFiringRateOverall);
        fprintf('\tComputing evoked spiking and writing file %s...\n', saveFileName);

        computeEvokedSpiking(saveFileName, unitStruct, nLoc, UE, numRandomizations);
    else
        fprintf('\tOverall firing rate = %0.2f Hz < minimum firing rate = %0.2f Hz in these blocks - skipping.\n', ...
                firingRateOverall, minFiringRateOverall);
        if exist(saveFileName, 'file')
            fprintf('\t%s exists... deleting.\n', saveFileName);
            delete(saveFileName);
        end 
    end
end

fprintf('Time elapsed: %0.2f s.\n\n', toc);

