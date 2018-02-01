% clear;
% 
% lfpsFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170127\20170127-all_merged-sort1all.pl2';
% animalName = 'M';
% sessionName = '20170127';
% areaName = 'PUL';
% blockNames = {'vepm1', 'vepm2', 'rfm1', 'g1', 'g2', 'g3', 'g4', 'g5', 'rest1', 'rfm2', 'vepm3'};
% blockInds = [1 2];

%%
clear;

pl2FileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170327\20170327-mcc-pul-d1_34.00mm_d2_34.00mm_allExceptRest_merged_noWB-sort1.pl2';
sessionName = 'M20170327';
areaName = 'PUL';
blockNames = {'vepm1', 'vepm2', 'aepm1', 'rfm1', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'rfm2', 'rfm3', 'vepm3', 'vepm4', 'rfm4', 'rest1'};
% blockInds = [4 15 16 18];
blockInds = [1 2];

%%
lfpVEPMapping(pl2FileName, sessionName, areaName, blockNames, blockInds);

