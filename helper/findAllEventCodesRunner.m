function findAllEventCodesRunner(dataDirRoot, recordingInfoFileName, sessionInd, blockIndices)
% LFP RF Mapping, all channels on a probe
% can't really do one channel at a time because of Common Average
% Referencing

%% setup and load data
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Find all event codes\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('Block indices: %s\n', blockIndices);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('------------------------\n');

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
R = recordingInfo(sessionInd);
sessionName = R.sessionName;
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, R.pl2FileName);
fprintf('Loading %s...\n', pl2FilePath);

%% load recording data
muaDataDirRoot = '';

isLoadLfp = 0;
isLoadMua = 0;
isLoadSpikes = 0;
isLoadSpkc = 0;
isLoadDirect = 0;

R.spikeChannelsToLoad = NaN;
R.muaChannelsToLoad = NaN;
R.lfpChannelsToLoad = NaN;
R.spkcChannelsToLoad = NaN;

D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, R.areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        R.spikeChannelPrefix, R.spikeChannelsToLoad, R.muaChannelsToLoad, R.lfpChannelsToLoad, R.spkcChannelsToLoad, R.directChannelsToLoad); 
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(R.blockNames(blockIndices), '-');
fprintf('Analyzing block names: %s.\n', blockName);

D = trimSpikeTimesAndEvents(D, blockIndices);

%% find event codes
[eventTimes,eventCodesPort1,eventCodesPort2] = findAllEventCodes(D.events);
fprintf('Event Codes: \n');
for i = 1:numel(eventTimes)
    fprintf('\t%0.5f\t\t%d\t\t%d\n', eventTimes(i), eventCodesPort1(i), eventCodesPort2(i));
end
fprintf('\n');

UE1 = unique(eventCodesPort1);
UE2 = unique(eventCodesPort2);
fprintf('Unique codes port 1 (N=%d):\n', numel(UE1));
for i = 1:numel(UE1)
    fprintf('%d, ', UE1(i));
end
fprintf('\n\n');
fprintf('Unique codes port 2 (N=%d):\n', numel(UE2));
for i = 1:numel(UE2)
    fprintf('%d, ', UE2(i));
end
fprintf('\n\n');

