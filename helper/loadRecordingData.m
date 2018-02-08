function [R, D, processedDataDir, blockName] = loadRecordingData(...
        processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, ...
        sessionInd, muaChannelsToLoad)

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
recordingInfo = rmfield(recordingInfo, 'muaChannelsToLoad'); % remove b/c we're using passed value
R = recordingInfo(sessionInd);
sessionName = R.sessionName;
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, R.pl2FileName);
fprintf('Loading %s...\n', pl2FilePath);

%% load recording data
isLoadSpikes = 0;
isLoadMua = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 1;

R.spikeChannelsToLoad = NaN;
R.muaChannelsToLoad = muaChannelsToLoad;
R.lfpChannelsToLoad = NaN;
R.spkcChannelsToLoad = NaN;

D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, R.areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        R.spikeChannelPrefix, R.spikeChannelsToLoad, R.muaChannelsToLoad, R.lfpChannelsToLoad, R.spkcChannelsToLoad, R.directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

assert(numel(R.blockNames) == numel(D.blockStartTimes));
blockName = strjoin(R.blockNames(R.gratingsTask3DIndices), '-');
fprintf('Analyzing block names: %s.\n', blockName);

%% remove spike and event times not during task to save memory
D = trimSpikeTimesAndEvents(D, R.gratingsTask3DIndices);
