clear;
sessionName = 'M20170311';
sessionInd = 7;
processedDataRootDir = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/';
dataDirRoot = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\';
muaDataDirRoot = ['C:\Users\Ryan\Documents\MATLAB\gratings-task-data\' sessionName];
recordingInfoFileName = 'C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/recordingInfo2.csv';
isZeroDistractors = 0;
numRandomizations = 2;

for muaChannelsToLoad = 1:32
    muaAnalysis(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, muaChannelsToLoad, isZeroDistractors, numRandomizations);
end