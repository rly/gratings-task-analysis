% clear;
% sessionName = 'M20170311';
% sessionInd = 7;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
% dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
% muaDataDirRoot = ['C:/Users/Ryan/Documents/MATLAB/gratings-task-data/' sessionName];
% suaMuaDataDirRoot = muaDataDirRoot;
% recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% channelsToLoad = 1:32;
% muaChannelsToLoad = 1:32;
% lfpChannelsToLoad = 1:32;
% lfpChannels = lfpChannelsToLoad;
% isZeroDistractors = 0;
% numRandomizations = 2;
% isLoadSortedSua = 1;
% isLoadMua = 0;

% 
% clear;
% sessionName = 'M20170311';
% sessionInd = 8;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
% dataDirRoot = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-data/';
% muaDataDirRoot = ['C:/Users/Ryan/Documents/MATLAB/gratings-task-data/' sessionName];
% suaMuaDataDirRoot = muaDataDirRoot;
% recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% channelsToLoad = 33:64;
% muaChannelsToLoad = 33:64;
% lfpChannelsToLoad = 33:64;
% lfpChannels = lfpChannelsToLoad;
% isZeroDistractors = 0;
% numRandomizations = 2;
% isLoadSortedSua = 1;
% isLoadMua = 0;

% clear;
% sessionName = 'M20170529';
% sessionInd = 24;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
% dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
% muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
% suaMuaDataDirRoot = muaDataDirRoot;
% recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% channelsToLoad = 1:32;
% muaChannelsToLoad = 1:32;
% lfpChannelsToLoad = 1:32;
% lfpChannels = lfpChannelsToLoad;
% isZeroDistractors = 0;
% numRandomizations = 2;
% isLoadSortedSua = 1;
% isLoadMua = 0;

% clear;
% sessionName = 'M20170201';
% sessionInd = 3;
% processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
% dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
% muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
% suaMuaDataDirRoot = muaDataDirRoot;
% recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
% channelsToLoad = 1:32;
% muaChannelsToLoad = 1:32;
% lfpChannelsToLoad = 1:32;
% lfpChannels = lfpChannelsToLoad;
% isZeroDistractors = 0;
% numRandomizations = 2;
% isLoadSortedSua = 1;
% isLoadMua = 0;

clear;
sessionName = 'M20170127';
sessionInd = 1;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = 1:32;
muaChannelsToLoad = 1:32;
lfpChannelsToLoad = 1:32;
lfpChannels = lfpChannelsToLoad;
isZeroDistractors = 0;
numRandomizations = 2;
isLoadSortedSua = 1;
isLoadMua = 0;

clear;
sessionName = 'M20170130';
sessionInd = 2;
processedDataRootDir = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = '/Users/labmanager/Documents/MATLAB/rlyRawData/';
muaDataDirRoot = ['/Users/labmanager/Documents/MATLAB/rlyRawData/' sessionName];
suaMuaDataDirRoot = muaDataDirRoot;
recordingInfoFileName = '/Users/labmanager/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
channelsToLoad = 1:32;
muaChannelsToLoad = channelsToLoad;
lfpChannelsToLoad = 1:32;
lfpChannels = lfpChannelsToLoad;
isZeroDistractors = 0;
numRandomizations = 2;
isLoadSortedSua = 1;
isLoadMua = 0;

%sessionInds = [26];
%suaMuaAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds)
%
%for chani = 1:length(channelsToLoad)
%    suaMuaAnalysisExtractSummaryData(processedDataRootDir, dataDirRoot, ...
%        suaMuaDataDirRoot, recordingInfoFileName, sessionInd, channelsToLoad, isLoadSortedSua, isLoadMua)
%end