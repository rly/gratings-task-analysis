function saveLfpForChethan(processedDataRootDir, dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, lfpChannels, isZeroDistractors)
% 325ms fixation before pre-cue marker
% 25-125ms fixation between pre-cue marker and cue onset
% 100ms cue onset to cue offset
% 500-800ms cue offset to array onset
% 850-1050ms array onset to target dim for long hold trials
% 650-850ms array onset to target dim for long hold trials
% hold shape response window 280-730ms
% release shape response window 280-730ms
% min 280ms before saccade allowed

% evt5 = cue onset
% evt1,2,3,4 = cue offset
% evt6 = array onset and also release shape resp win start
% evt7 = target dim
% evt8 = juice

v = 12;
tic;

fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task LFP Analysis\n');
fprintf('Session index: %d\n', sessionInd);
fprintf('LFP Channel to Load: %d\n', lfpChannels);
fprintf('Recording info file name: %s\n', recordingInfoFileName);
fprintf('Processed data root dir: %s\n', processedDataRootDir);
fprintf('Data root dir: %s\n', dataDirRoot);
fprintf('MUA data root dir: %s\n', muaDataDirRoot);
fprintf('Is zero distractors: %s\n', isZeroDistractors);
fprintf('Version: %d\n', v);
fprintf('------------------------\n');

nLoc = 4;

%% input check
% assert(numel(lfpChannelsToLoad) == 1);
assert(numel(lfpChannels) > 1);
assert(numel(lfpChannels) <= 32);

%% load recording information
if isZeroDistractors
    scriptName = 'LFP_GRATINGS_0D';
    taskName = 'GRATINGS_0D';
else
    scriptName = 'LFP_GRATINGS';
    taskName = 'GRATINGS';
end

[R, D, processedDataDir, blockName] = loadRecordingData(processedDataRootDir, ...
        dataDirRoot, muaDataDirRoot, recordingInfoFileName, sessionInd, lfpChannels, ...
        taskName, scriptName, 1, 1);
sessionName = R.sessionName;
areaName = R.areaName;

fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('%s/%s/', dataDirRoot, sessionName);
gratingsTaskLogDir = sprintf('%s/%s', dataDir, sessionName(2:end));

% process events and sort them into different conditions
UE = getUsefulEvents2(gratingsTaskLogDir, R.gratingsTaskLogIndices, 4, D, blockName);

fileNamePrefix = sprintf('%s-ch%d-ch%d-%s', sessionName, lfpChannels([1 end]), blockName);

%% preprocess LFPs
tic;
Fs = D.lfpFs;
nChannels = D.nLfpCh;

D.adjLfpsClean = interpolateLfpOverSpikeTimes(D.adjLfps, lfpChannels, Fs, D.allMUAStructs);

hiCutoffFreq = 200;
[channelDataNorm,isTrialOutlier,isNoisyChannel,meanBeforeCAR] = preprocessLfpsGratingsTask(D.adjLfpsClean, Fs, D.lfpNames, UE, ...
        processedDataDir, fileNamePrefix, hiCutoffFreq, v); % TODO USE isTrialOutlier
D.adjLfps = [];
D.adjLfpsClean = [];
D.adjDirects = [];

% LAST CHANNEL IS MEAN BEFORE COMMON AVERAGE REFERENCING
channelDataNorm(size(channelDataNorm, 1)+1,:) = meanBeforeCAR;

channelNames = cell(nChannels + 1, 1);
for i = 1:nChannels
    channelNames{i} = sprintf('FP%03d', lfpChannels(i));
end
channelNames{nChannels + 1} = 'Common Average';

fprintf('Took %0.1f minutes.\n', toc/60);

%% save for chethan/feng collaboration
saveFileName = sprintf('%s/%s-continuousLfps-v%d.mat', ...
        processedDataDir, fileNamePrefix, v);
save(saveFileName, 'channelDataNorm', 'channelNames', 'Fs', 'hiCutoffFreq', ...
        'isNoisyChannel', 'isTrialOutlier', 'UE', 'v', 'sessionName', 'lfpChannels', 'D');
fileInfo = dir(saveFileName);
fprintf('%d MB file created. Took %0.1f minutes.\n', round(fileInfo.bytes/1024^2), toc/60);
