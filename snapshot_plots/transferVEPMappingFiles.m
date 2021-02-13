%% copy LFP VEPM files
readDataRemotely;
sessionInds = [1:16 19 20 23];
sessionInds = 39:57
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

%% copy MUA VEPM files
readDataRemotely;
sessionInds = [1:16 19 20 23];
v = 11;

processedDataRootDirLocal = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';

recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

%% session loop
for i = 9:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    muaChannelsToLoad = R.muaChannelsToLoad;
    blockName = strjoin(R.blockNames(R.vepmIndices), '-');
    processedDataDir = sprintf('%s/%s/MUA_VEPM/', processedDataRootDir, sessionName);
    for j = 1:numel(muaChannelsToLoad)
        c = muaChannelsToLoad(j);
        sourceFileName = sprintf('%s/%s_%s_%dM-%s-vepm-v%d.mat', processedDataDir, sessionName, areaName, c, blockName, v);

        processedDataDirLocal = sprintf('%s/%s/MUA_VEPM/', processedDataRootDirLocal, sessionName);
        if ~exist(processedDataDirLocal, 'dir')
            mkdir(processedDataDirLocal);
        end
        destFileName = sprintf('%s/%s_%s_%dM-%s-vepm-v%d.mat', processedDataDirLocal, sessionName, areaName, c, blockName, v);

        fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
        if exist(sourceFileName, 'file')
            copyfile(sourceFileName, destFileName)
        else
            warning('File %s does not exist', sourceFileName);
        end
    end
end