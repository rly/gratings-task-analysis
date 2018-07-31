readDataRemotely;
sessionInds = [1:16 19 20 23];
ref = 'RAW';
v = 12;

processedDataRootDirLocal = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';

recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    lfpChannelsToLoad = R.lfpChannelsToLoad;
    blockName = strjoin(R.blockNames(R.vepmIndices), '-');
    processedDataDir = sprintf('%s/%s/LFP_VEPM/', processedDataRootDir, sessionName);
    fileNamePrefix = sprintf('%s-ind%d-%s-ch%d-ch%d-%s', sessionName, sessionInd, areaName, lfpChannelsToLoad([1 end]), blockName);
    sourceFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDir, fileNamePrefix, ref, v);

    processedDataDirLocal = sprintf('%s/%s/LFP_VEPM/', processedDataRootDirLocal, sessionName);
    mkdir(processedDataDirLocal);
    destFileName = sprintf('%s/%s-%s-responses-v%d.mat', processedDataDirLocal, fileNamePrefix, ref, v);
    
    fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
    copyfile(sourceFileName, destFileName)
end