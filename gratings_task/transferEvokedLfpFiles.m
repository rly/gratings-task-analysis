readDataRemotely;
sessionInds = 1:37;
v = 13;

processedDataRootDirLocal = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/LFP_GRATINGS_ALL/';
mkdir(processedDataRootDirLocal);

recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    blockName = strjoin(R.blockNames(R.gratingsTask3DIndices), '-');
    processedDataDir = sprintf('%s/%s/LFP_GRATINGS/', processedDataRootDir, sessionName);
    fileNameTemplate = sprintf('%s/%s-ch%d-ch%d-%s-evokedLfps-v%d*', processedDataDir, sessionName, R.lfpChannelsToLoad([1 end]), blockName, v);
    fileMatches = dir(fileNameTemplate);
    assert(numel(fileMatches) == 1);
    sourceFileName = sprintf('%s/%s', processedDataDir, fileMatches(1).name);
    destFileName = sprintf('%s/%s', processedDataRootDirLocal, fileMatches(1).name);
    fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
    copyfile(sourceFileName, destFileName);
end