readDataRemotely;
sessionInds = 38:57;
v = 12;

processedDataRootDirLocal = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/Ferdy_LFP_VEPM/';
mkdir(processedDataRootDirLocal);

recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    blockName = strjoin(R.blockNames(R.vepmIndices), '-');
    processedDataDir = sprintf('%s/%s/LFP_VEPM/', processedDataRootDir, sessionName);
    
    fileNameTemplate = sprintf('%s/%s-ind%d-PUL-ch%d-ch%d-%s-CAR-lfpColor-v%d.png', ...
            processedDataDir, sessionName, sessionInd, R.lfpChannelsToLoad([1 end]), blockName, v);
    fileMatches = dir(fileNameTemplate);
    assert(numel(fileMatches) == 1);
    sourceFileName = sprintf('%s/%s', processedDataDir, fileMatches(1).name);
    destFileName = sprintf('%s/%s', processedDataRootDirLocal, fileMatches(1).name);
    fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
    copyfile(sourceFileName, destFileName);
    
    fileNameTemplate = sprintf('%s/%s-ind%d-PUL-ch%d-ch%d-%s-RAW-lfpColor-v%d.png', ...
            processedDataDir, sessionName, sessionInd, R.lfpChannelsToLoad([1 end]), blockName, v);
    fileMatches = dir(fileNameTemplate);
    assert(numel(fileMatches) == 1);
    sourceFileName = sprintf('%s/%s', processedDataDir, fileMatches(1).name);
    destFileName = sprintf('%s/%s', processedDataRootDirLocal, fileMatches(1).name);
    fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
    copyfile(sourceFileName, destFileName);
end