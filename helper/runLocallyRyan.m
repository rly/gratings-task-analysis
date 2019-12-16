clear;
sessionName = 'M20170127';
sessionInd = 1;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
muaDataDirRoot = ['C:/Users/Ryan/Documents/MATLAB/gratings-task-data/' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% channelsToLoad = 1:32;
% muaChannelsToLoad = 1:32;
% lfpChannelsToLoad = 1:32;
% lfpChannels = lfpChannelsToLoad;
isZeroDistractors = 0;
numRandomizations = 2;
isLoadSortedSua = 1;
isLoadMua = 0;

% clear;
% sessionName = 'F20171009';
% sessionInd = 40;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
% dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
% muaDataDirRoot = ['C:/Users/Ryan/Documents/MATLAB/gratings-task-data/' sessionName];
% suaMuaDataDirRoot = muaDataDirRoot;
% recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% % channelsToLoad = 1:32;
% % muaChannelsToLoad = 1:32;
% % lfpChannelsToLoad = 1:32;
% % lfpChannels = lfpChannelsToLoad;
% isZeroDistractors = 0;
% numRandomizations = 2;
% isLoadSortedSua = 1;
% isLoadMua = 0;
% 
% saveGratingsTaskResultsJsonToMatRunner(dataDirRoot, recordingInfoFileName, sessionInd, isZeroDistractors)
% 
for chan = 1:32
    suaMuaAnalysis(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, chan, ...
            isZeroDistractors, numRandomizations, isLoadSortedSua, isLoadMua)
    suaMuaAnalysisPlots(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, chan, ...
            isZeroDistractors, isLoadSortedSua, isLoadMua)
end
% 
% channelsToLoad = 1:32;
% loadRecordingDataIntoSpikeMetaData(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, ...
%         sessionInd, channelsToLoad, isLoadSortedSua, isLoadMua)
% suaMuaAnalysisExtractSummaryData(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, ...
%         sessionInd, channelsToLoad, isLoadSortedSua, isLoadMua)
% 
% 
% suaMuaAnalysis(processedDataRootDir, dataDirRoot, suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, ...
%         isZeroDistractors, numRandomizations, isLoadSortedSua, isLoadMua)