function MD = loadRecordingDataIntoSpikeMetaData(paramsStruct)
% read SUA/MUA data from PL2 file, trim the spike times based on the blocks
% of interest, and save only the cell array of unit structs and the block
% start and stop times
%
% do this once to save time loading the full POL2 file for the sua/mua
% analysis extract summary script and other summary scripts
%% check params struct for all required fields
% the fields in paramsStruct will be unpacked into the workspace
requiredFields = {'processedDataRootDir', 'dataDirRoot', 'suaMuaDataDirRoot', 'recordingInfoFileName', ...
        'sessionInd', 'channelsToLoad', 'isLoadSortedSua', 'isLoadMua', 'minSuaSepQuality'};
% certain fields are used only for loading certain data
fn = fieldnames(paramsStruct);
for i = 1:numel(requiredFields)
    if ~ismember(requiredFields{i}, fn)
        error(['loadRecordingDataIntoSpikeMetaData input struct is missing field "' requiredFields{i} '"'])
    end
end
struct2varfn(paramsStruct)

taskName = 'GRATINGS';
scriptName = 'SUA_MUA_GRATINGS';

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
isLoadDirect = 1;
isLoadLfp = 0;

assert(numel(channelsToLoad) == 32);
if ~isempty(channelsToLoad)
    R.spikeChannelsToLoad = channelsToLoad;
    R.muaChannelsToLoad = channelsToLoad;
    R.lfpChannelsToLoad = channelsToLoad;
end

processedDataDirPre = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDirPre, 'dir') == 0
    mkdir(processedDataDirPre);
end
processedDataDir = sprintf('%s/%s', processedDataDirPre, scriptName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end

tic;
fprintf('Loading data %s...\n', pl2FilePath);
D = loadPL2(pl2FilePath, suaMuaDataDirRoot, sessionName, R.areaName, isLoadSortedSua, isLoadMua, isLoadLfp, ...
            isLoadSpkc, isLoadDirect, R.spikeChannelPrefix, R.spikeChannelsToLoad, R.muaChannelsToLoad, ...
            R.lfpChannelsToLoad, R.spkcChannelsToLoad, R.directChannelsToLoad, minSuaSepQuality);  
fprintf('... done (%0.2f s).\n', toc);

%% get block indices
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
else
    error('Unknown task name: %s\n', taskName);
end

blockName = strjoin(R.blockNames(R.blockIndices), '-');
fprintf('Analyzing task name: %s, block names: %s.\n', taskName, blockName);

%% remove spike and event times not during task to save memory
D = trimSpikeTimesAndEvents(D, R.blockIndices);
if isLoadLfp
    D = adjustSpikeTimesLfpsAndEvents(D, R.blockIndices);
end

%% save meta data into smaller file
R.metaDataFilePath = sprintf('%s/%s-sessionInd%d-sua%d-minsq%d-mua%d-gratings-metadata.mat', ...
        processedDataDir, sessionName, sessionInd, isLoadSortedSua, minSuaSepQuality, isLoadMua);
fprintf('Writing metadata %s...\n', R.metaDataFilePath);
MD = createMetaDataFile(D, R.metaDataFilePath);
