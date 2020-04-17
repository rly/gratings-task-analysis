function [R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, ...
        sessionInd, channelsToLoad, taskName, scriptName, isLoadSortedSua, isLoadMua, isLoadLfp, isLoadMetaDataOnly, ...
        rfMappingNewInfoFileName, rfMappingNewMode, isLoadAllSpikes)
% loads MUA data and eyetracking/lever data into D struct and recording
% metadata into R struct

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
if ~isempty(channelsToLoad)
    recordingInfo = rmfield(recordingInfo, 'lfpChannelsToLoad'); % remove b/c we're using passed value
    recordingInfo = rmfield(recordingInfo, 'muaChannelsToLoad'); % remove b/c we're using passed value
    recordingInfo = rmfield(recordingInfo, 'spikeChannelsToLoad'); % remove b/c we're using passed value
end
R = recordingInfo(sessionInd);
sessionName = R.sessionName;
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, R.pl2FileName);

%% load recording data
isLoadSpkc = 0;

% if strcmp(taskName, 'GRATINGS') || strcmp(taskName, 'GRATINGS_0D')
    isLoadDirect = 1;
% else
%     isLoadDirect = 0;
% end

if ~isempty(channelsToLoad)
    R.spikeChannelsToLoad = channelsToLoad;
    R.muaChannelsToLoad = channelsToLoad;
    R.lfpChannelsToLoad = channelsToLoad;
end

processedDataDirPre = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDirPre, 'dir') == 0
%    mkdir(processedDataDirPre);
end
processedDataDir = sprintf('%s/%s', processedDataDirPre, scriptName);
if exist(processedDataDir, 'dir') == 0
%     mkdir(processedDataDir);
end

tic;
if isLoadMetaDataOnly
    R.metaDataFilePath = sprintf('%s/%s-sessionInd%d-sua%d-mua%d-gratings-metadata.mat', ...
            processedDataDir, sessionName, sessionInd, isLoadSortedSua, isLoadMua);
    fprintf('Loading metadata %s...\n', R.metaDataFilePath);
    MD = load(R.metaDataFilePath);
    D = MD.MD;
else
    fprintf('Loading data %s...\n', pl2FilePath);
    D = loadPL2(pl2FilePath, suaMuaDataDirRoot, sessionName, R.areaName, isLoadSortedSua, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
            R.spikeChannelPrefix, R.spikeChannelsToLoad, R.muaChannelsToLoad, R.lfpChannelsToLoad, R.spkcChannelsToLoad, R.directChannelsToLoad); 
end
fprintf('... done (%0.2f s).\n', toc);

%% get block indices
fprintf('%d block names specified, %d entries in block start times.\n', numel(R.blockNames), numel(D.blockStartTimes));
assert(numel(R.blockNames) == numel(D.blockStartTimes));
if strcmp(taskName, 'GRATINGS')
    R.blockIndices = R.gratingsTask3DIndices;
    R.gratingsTaskLogIndices = R.gratingsTask3DLogIndices;
elseif strcmp(taskName, 'GRATINGS_0D')
    R.blockIndices = R.gratingsTask0DIndices;
    R.gratingsTaskLogIndices = R.gratingsTask0DLogIndices;
    if isnan(R.blockIndices)
        error('No Block Indices defined for Gratings Task 0D');
    end
elseif strcmp(taskName, 'VEPM')
    R.blockIndices = R.vepmIndices;
elseif strcmp(taskName, 'AEPM')
    R.blockIndices = R.aepmIndices;
elseif strcmp(taskName, 'RFM_OLD')
    R.blockIndices = R.rfmOldIndices;
    if isnan(R.blockIndices)
        error('No Block Indices defined for RF Mapping Old Task');
    end
elseif strcmp(taskName, 'RFM_EIGHTHS')
    R.blockIndices = R.rfmEighthsIndices;
    if isnan(R.blockIndices)
        error('No Block Indices defined for RF Mapping Eighths Task');
    end
elseif strcmp(taskName, 'RFM_NEW')
    rfMappingNewInfo = readRFMappingNewInfo(rfMappingNewInfoFileName);
    matchSession = cellfun(@(x) strcmp(x, R.sessionName), {rfMappingNewInfo.sessionName});
    matchMode = [rfMappingNewInfo.mode] == rfMappingNewMode;
    R.blockIndices = [rfMappingNewInfo(matchSession & matchMode).blockInd];
    R.rfmResultsRootDir = sprintf('%s/%s/%s', dataDirRoot, sessionName, sessionName(2:end));
    R.rfmResultsFileNames = {rfMappingNewInfo(matchSession & matchMode).resultsFileName};
else
    error('Unknown task name: %s\n', taskName);
end

blockName = strjoin(R.blockNames(R.blockIndices), '-');
fprintf('Analyzing task name: %s, block names: %s.\n', taskName, blockName);

%% remove spike and event times not during task to save memory
% meta data already had spike times adjusted
if ~isLoadMetaDataOnly
    D = trimSpikeTimesAndEvents(D, R.blockIndices);
    if isLoadLfp
        D = adjustSpikeTimesLfpsAndEvents(D, R.blockIndices);
    end
end
% %% remove spike and event times not during task to save memory
% % meta data already had spike times adjusted
% if ~isLoadMetaDataOnly 
%     if isLoadAllSpikes
%         DallSpikes = D;
%         [D, spikesToKeep] = trimSpikeTimesAndEvents(D, R.blockIndices, isLoadAllSpikes);
%         for uniti = 1:numel(DallSpikes.allUnitStructs)
%             D.allUnitStructs{uniti}.tsAll = DallSpikes.allUnitStructs{uniti}.ts;
%             D.spikesToKeep = spikesToKeep;
%         end
%     else
%         [D, spikesToKeep] = trimSpikeTimesAndEvents(D, R.blockIndices, isLoadAllSpikes);
%     end
%     if isLoadLfp
%         D = adjustSpikeTimesLfpsAndEvents(D, R.blockIndices);
%     end
% end