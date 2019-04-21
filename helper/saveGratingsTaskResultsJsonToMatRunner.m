function saveGratingsTaskResultsJsonToMatRunner(dataDirRoot, recordingInfoFileName, sessionInds, isZeroDistractors)

fprintf('\n-------------------------------------------------------\n');
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Session indices: %d\n', sessionInds);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('------------------------\n');

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
for sessionInd = sessionInds
    tic;
    % note this is not efficient because when there are multiple 
    % sessionInds for the same experimental session, those sessionInds will
    % have the same log file and this will generate the same mat file
    % multiple times
    R = recordingInfo(sessionInd);
    if ~isZeroDistractors
        R.blockIndices = R.gratingsTask3DIndices;
        R.gratingsTaskLogIndices = R.gratingsTask3DLogIndices;
    else
        R.blockIndices = R.gratingsTask0DIndices;
        R.gratingsTaskLogIndices = R.gratingsTask0DLogIndices;
        if isnan(R.blockIndices)
            error('No Block Indices defined for Gratings Task 0D');
        end
    end
    blockName = strjoin(R.blockNames(R.blockIndices), '-');
    sessionName = R.sessionName;
    
    fprintf('------------------------\n');
    fprintf('Processing session index %d, %s...\n', sessionInd, sessionName);

    dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
    gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

    %% process events and sort them into different conditions
    saveFileName = sprintf('%s/gratingsTaskResultsJson-%s.mat', gratingsTaskLogDir, blockName);
    saveGratingsTaskResultsJsonToMat(gratingsTaskLogDir, R.gratingsTaskLogIndices, saveFileName);
    toc
end