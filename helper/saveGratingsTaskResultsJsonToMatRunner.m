function saveGratingsTaskResultsJsonToMatRunner(dataDirRoot, recordingInfoFileName, sessionInd, isZeroDistractors)
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task MUA Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('------------------------\n');

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
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

fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

%% process events and sort them into different conditions
saveFileName = sprintf('%s/gratingsTaskResultsJson-%s.mat', gratingsTaskLogDir, blockName);
readGratingsTaskResultsJson(gratingsTaskLogDir, R.gratingsTaskLogIndices, saveFileName);