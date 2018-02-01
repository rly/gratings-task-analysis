clear;

spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\20170127\20170127-all_merged-sort1.pl2';
animalName = 'M';
sessionName = '20170127';
areaName = 'PUL';
blockNames = {'vepm1', 'vepm2', 'rfm1', 'g1', 'g2', 'g3', 'g4', 'g5', 'rest1', 'rfm2', 'vepm3'};

%%

fprintf('\n-------------------------------------------------------\n');
fprintf('RF Mapping Analysis - Spiking\n');
fprintf('Loading %s...\n', spikesFileName);
tic;

isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
D = loadPL2(spikesFileName, sessionName, animalName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 

processedDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/%s%s', animalName, sessionName);
if exist(processedDir, 'dir') == 0
    mkdir(processedDir);
end
fprintf('... done (%0.2f s).\n', toc);

%%
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);

%%
for j = 1:nUnits
    spikeStruct = D.allSpikeStructs{j};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
%     fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, j, ...
%             nUnits, round(j/nUnits*100));
        
    fsClassThresh = 0.35/1000;
    if startsWith(spikeStruct.inflectionPattern, 'tp')
        if spikeStruct.troughToPeakTime <= fsClassThresh
            fsClass = 'Narrow';
        else
            fsClass = 'Broad';
        end
    else
        fsClass = '';
    end
    fprintf('Unit #%d: %s - %s %s\n', j, spikeStruct.name, spikeStruct.inflectionPattern, fsClass);

end