readDataRemotely;
sessionInds = 1:37;
v = 12;

processedDataRootDirLocal = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/PUL_MUA_GRATINGS_ALL/';
mkdir(processedDataRootDirLocal);

recordingInfo = readRecordingInfo(recordingInfoFileName);
nSessions = numel(sessionInds);

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    areaName = R.areaName;
    pulChannels = [R.dPulChannels R.vPulChannels];
    blockName = strjoin(R.blockNames(R.gratingsTask3DIndices), '-');
    processedDataDir = sprintf('%s/%s/MUA_GRATINGS/', processedDataRootDir, sessionName);
    fileNameTemplate = sprintf('%s/%s_%s_*-%s-*-v%d*', processedDataDir, sessionName, areaName, blockName, v);
    fileMatches = dir(fileNameTemplate);
    for j = 1:numel(fileMatches)
        tokens = regexp(fileMatches(j).name, '.+_.+_(\d)+\w-.+-.+-v\d+.+', 'tokens');
        assert(numel(tokens) == 1);
        channelID = str2double(tokens{1}{1});
        
        if find(channelID == pulChannels)
            sourceFileName = sprintf('%s/%s', processedDataDir, fileMatches(j).name);
            destFileName = sprintf('%s/%s', processedDataRootDirLocal, fileMatches(j).name);
            fprintf('Copying file %s to %s...\n', sourceFileName, destFileName);
            copyfile(sourceFileName, destFileName);
        end
    end
end