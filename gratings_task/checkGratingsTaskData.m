

clear;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'Z:/ryanly/McCartney/originals';
muaDataDirRoot = processedDataRootDir;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';

v = 9;
sessionIndAll = 1:37;
for k = 1:numel(sessionIndAll)
    clearvars -except processedDataRootDir dataDirRoot sessionIndAll k v % cheap hack
    sessionInd = sessionIndAll(k);

nLoc = 4;

%% load recording information
recordingInfo = readRecordingInfo(recordingInfoFileName);
struct2var(recordingInfo(sessionInd));
pl2FilePath = sprintf('%s/%s/%s', dataDirRoot, sessionName, pl2FileName);

%% load recording data
fprintf('\n-------------------------------------------------------\n');
fprintf('Check Data for Gratings Task Analysis\n');
fprintf('Loading %s...\n', pl2FilePath);
fprintf('Session index: %d\n', sessionInd);

tic;
isLoadSpikes = 0;
isLoadMua = 0;
isLoadLfp = 0;
isLoadSpkc = 0;
isLoadDirect = 1;
D = loadPL2(pl2FilePath, muaDataDirRoot, sessionName, areaName, isLoadSpikes, isLoadMua, isLoadLfp, isLoadSpkc, isLoadDirect, ...
        spikeChannelPrefix, spikeChannelsToLoad, muaChannelsToLoad, lfpChannelsToLoad, spkcChannelsToLoad, directChannelsToLoad); 

processedDataDir = sprintf('%s/%s', processedDataRootDir, sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

assert(numel(blockNames) == numel(D.blockStartTimes));
blockName = strjoin(blockNames(gratingsTask3DIndices), '-');
fprintf('Analyzing block names: %s.\n', blockName);

%% remove spike and event times not during task to save memory
D = trimSpikeTimesAndEvents(D, gratingsTask3DIndices);

%%  
% for l = sessionRange
fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
logDir = sprintf('%s/%s', dataDir, sessionName(2:end));
load('params.mat');

% process events and sort them into different conditions
UE = getUsefulEvents2(logDir, gratingsTask3DLogIndices, 4, D);

end
