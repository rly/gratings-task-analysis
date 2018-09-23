function [gain,dataInfo] = getPL2AnalogGain(fileName)
% function to remind self how to find analog gain of PL2 file

dataInfo = PL2GetFileIndex(fileName);
gain = dataInfo.AnalogChannels{1}.TotalGain;