function [R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, taskName, scriptName, isLoadMua, isLoadLfp, rfMappingNewInfoFileName, rfMappingNewMode)
% loads MUA data and eyetracking/lever data into D struct and recording
% metadata into R struct

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
recordingInfo = rmfield(recordingInfo, 'lfpChannelsToLoad'); % remove b/c we're using passed value
recordingInfo = rmfield(recordingInfo, 'muaChannelsToLoad'); % remove b/c we're using passed value
R = recordingInfo(sessionInd);
sessionName = R.sessionName;
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, R.pl2FileName);
fprintf('Loading %s...\n', pl2FilePath);

%% load recording data
isLoadSpikes = 0;
isLoadSpkc = 0;
isLoadDirect = 1;

R.spikeChannelsToLoad = NaN;
R.muaChannelsToLoad = channelsToLoad;
R.lfpChannelsToLoad = channelsToLoad;
R.spkcChannelsToLoad = NaN;

D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, R.areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        R.spikeChannelPrefix, R.spikeChannelsToLoad, R.muaChannelsToLoad, R.lfpChannelsToLoad, R.spkcChannelsToLoad, R.directChannelsToLoad); 

processedDataDirPre = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDirPre, 'dir') == 0
    mkdir(processedDataDirPre);
end
processedDataDir = sprintf('%s/%s', processedDataDirPre, scriptName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

%% get block indices
assert(numel(R.blockNames) == numel(D.blockStartTimes));
if strcmp(taskName, 'Gratings')
    blockIndices = R.gratingsTask3DIndices;
elseif strcmp(taskName, 'VEPM')
    blockIndices = R.vepmIndices;
elseif strcmp(taskName, 'RFM_OLD')
    blockIndices = R.rfmOldIndices;
    if isnan(blockIndices)
        error('No Block Indices defined for RF Mapping Old Task');
    end
elseif strcmp(taskName, 'RFM_NEW')
    rfMappingNewInfo = readRFMappingNewInfo(rfMappingNewInfoFileName);
    matchSession = cellfun(@(x) strcmp(x, R.sessionName), {rfMappingNewInfo.sessionName});
    matchMode = [rfMappingNewInfo.mode] == rfMappingNewMode;
    blockIndices = [rfMappingNewInfo(matchSession & matchMode).blockInd];
    R.rfmResultsRootDir = sprintf('%s/%s/%s', dataDirRoot, sessionName, sessionName(2:end));
    R.rfmResultsFileNames = {rfMappingNewInfo(matchSession & matchMode).resultsFileName};
else
    error('Unknown task name: %s\n', taskName);
end

blockName = strjoin(R.blockNames(blockIndices), '-');
fprintf('Analyzing task name: %s, block names: %s.\n', taskName, blockName);

%% remove spike and event times not during task to save memory
D = trimSpikeTimesAndEvents(D, blockIndices);
if isLoadLfp
    D = adjustSpikeTimesLfpsAndEvents(D, blockIndices);
end