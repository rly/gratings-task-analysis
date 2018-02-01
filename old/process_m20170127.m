% running code for M 2017-01-27

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

% 1105 g3-g6, ch65-69 above brain
% 1104 g3-g6, ch65-66 above brain
% 1102 g2-g9, ch65-68 above brain
% 1101 g2-g3, ch65-67 above brain
clear;

sessionRange = 1:1;

sessions = {...
        'M20170127', 5:7
%         'M20170327', 2:5%[2:5 7:10]
%         'M20170529', 2:4
        };
numSessions = size(sessions, 1);

numFreqInS400 = 92;
numFreqInS250 = 46;
SPreArrayInRFAll = nan(numSessions, numFreqInS400);
SPreArrayExRFAll = nan(numSessions, numFreqInS400);
SPreDimInRFAll = nan(numSessions, numFreqInS400);
SPreDimExRFAll = nan(numSessions, numFreqInS400);
SPreCueInRFAll = nan(numSessions, numFreqInS250);

nChannels = 32;
startChannel = 65;

supChannelOrig = 6;
deepChannelOrig = 27;

inRFLoc = 3;
exRFLoc = 1;

nLoc = 4;

logRoot = 'gratings_task_eye_tracking_';
% cohCueOnsetAll = cell(numSessions, nLoc, nChannels);
% cohArrayOnsetAll = cell(numSessions, nLoc, nChannels);
% cohArrayOnsetHoldAll = cell(numSessions, nLoc, nChannels);
% cohTargetDimAll = cell(numSessions, nLoc, nChannels);
% 
% cohLineCueOnsetAll = cell(numSessions, nLoc, nChannels);
% cohLineArrayOnsetAll = cell(numSessions, nLoc, nChannels);
% cohLineArrayOnsetHoldAll = cell(numSessions, nLoc, nChannels);
% cohLineTargetDimAll = cell(numSessions, nLoc, nChannels);

%%
spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170127\20170127-all_merged-sort2.pl2';
sessionName = 'M20170127';
areaName = 'PUL';
blockNames = {'vepm1', 'vepm2', 'rfm1', 'g1', 'g2', 'g3', 'g4', 'g5', 'rest1', 'rfm2', 'vepm3'};
blockInds = 6:8;

%%
% spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170327\20170327-mcc-pul-d1_34.00mm_d2_34.00mm_allExceptRest_merged_noWB-sort1.pl2';
% sessionName = 'M20170327';
% areaName = 'PUL';
% blockNames = {'vepm1', 'vepm2', 'aepm1', 'rfm1', 'g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'rfm2', 'rfm3', 'vepm3', 'vepm4', 'rfm4', 'rest1'};
% blockInds = 6:9;%6:13;

%%
% spikesFileName = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-data\M20170529\d1-33.00mm-d2-32.00mm-all_merged_noWB-sort.pl2';
% sessionName = 'M20170529';
% areaName = 'PUL';
% blockNames = {'vepm4', 'vepm5', 'aepm1', 'aepm2', 'aepm3', 'rfm1', 'rfm2', 'rfm3', 'rfm4', 'g1', 'g2', 'g3', 'g4', 'rfm5', 'rest1'};
% blockInds = 11:13;

%%
fprintf('\n-------------------------------------------------------\n');
fprintf('Gratings Task Analysis\n');
fprintf('Loading %s...\n', spikesFileName);
tic;
isLoadSpikes = 1;
isLoadLfp = 0;
isLoadSpkc = 0;
D = loadPL2(spikesFileName, sessionName, areaName, isLoadSpikes, isLoadLfp, isLoadSpkc); 

processedDataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data/%s', sessionName);
if exist(processedDataDir, 'dir') == 0
    mkdir(processedDataDir);
end
fprintf('... done (%0.2f s).\n', toc);

blockName = strjoin(blockNames(blockInds), '-');

%% remove spike and event times not during RFM task to save memory
D = trimSpikeTimesAndEvents(D, blockInds);

%%
% for l = sessionRange
l = 1;
    sessionName = sessions{l,1};
    logIndices = sessions{l,2};
    fprintf('Processing %s...\n', sessionName);

dataDir = sprintf('C:/Users/Ryan/Documents/MATLAB/gratings-task-data/%s/', ...
        sessionName);
logDir = sprintf('%s/%s', dataDir, sessionName(2:end));
load('params.mat')

% process events and sort them into different conditions
usefulEvents = getUsefulEvents2(logDir, logIndices, 4, D);

% unload these variables to workspace
% (cueOnset, cueOnsetByLoc, ...
%         arrayOnset, arrayOnsetRel, arrayOnsetHold, arrayOnsetByLoc, ...
%         arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%         arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%         targetDim, targetDimByLoc, ...
%         targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
%         nTrialShortHold, nTrialLongHold, rt);
struct2var(usefulEvents);
stop

%% process every cell in workspace
nUnits = numel(D.allSpikeStructs);
fprintf('Processing %d units...\n', nUnits);

tic;
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    if ~isempty(spikeTimes)
        process_spikes_m20170127(processedDataDir, blockName, spikeStruct, nLoc, ...
                arrayOnsetRel, arrayOnsetHold, ...
                arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
                targetDim, targetDimByLoc, ...
                targetDimShortHoldByLoc, targetDimLongHoldByLoc, ...
                cueOnset, cueOnsetByLoc);
    % sfc with pulvinar
%     process_sfc_m20161105(sessionName, eval(spikeVars{i}), spikeVars{i}, ...
%             adjLfp, lfpVarName, nLoc, params, ...
%             cueOnset, cueOnsetByLoc, ...
%             arrayOnset, arrayOnsetByLoc, ...
%             arrayOnsetRel, arrayOnsetHold, ...
%             arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%             arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%             targetDim, targetDimByLoc, ...
%             targetDimShortHoldDurByLoc, targetDimLongHoldDurByLoc);
    % local sfc
%     lfpVarName = sprintf('FP%s', spikeVars{i}(4:6));
%     adjLfp = padNaNsToAdjustLfpOffset(eval(lfpVarName), lfpTs, lfpInd, Fs);
%         process_sfc_m20170127(processedDataDir, blockName, spikeStruct, ...
%                 D.adjLfps(spikeStruct.channelID,:)', ...
%                 sprintf('FP%03d', spikeStruct.channelID), nLoc, params, ...
%                 cueOnset, cueOnsetByLoc, ...
%                 arrayOnset, arrayOnsetByLoc, ...
%                 arrayOnsetRel, arrayOnsetHold, ...
%                 arrayOnsetRelByLoc, arrayOnsetHoldByLoc, ...
%                 arrayOnsetShortHoldByLoc, arrayOnsetLongHoldByLoc, ...
%                 targetDim, targetDimByLoc, ...
%                 targetDimShortHoldByLoc, targetDimLongHoldByLoc);
    end
%     close all;
end
fprintf('... done (%0.2f s).\n', toc);
stop

%% for grant
nUnits = numel(D.allSpikeStructs);
kernelSigma = 0.01;
fsClassThresh = 0.35/1000;

periArrayOnsetWindow = [0.7 0.5];
preArrayOnsetSpdfWindowOffset = [-0.5 0];
postArrayOnsetSpdfWindowOffset = [0 0.3];

periTargetDimWindow = [0.7 0.2];
targetDimSpdfWindowOffset = [-0.5 0];

nsCueTargetDelayActivity = nan(nUnits, 1);
bsCueTargetDelayActivity = nan(nUnits, 1);
nsArrayResponseActivity = nan(nUnits, 1);
bsArrayResponseActivity = nan(nUnits, 1);
nsTargetDimDelayActivity = nan(nUnits, 1);
bsTargetDimDelayActivity = nan(nUnits, 1);

unitCategory = ones(nUnits, 1);

nsPreArraySpdfInRFAll = [];
nsPreArraySpdfExRFAll = [];
bsPreArraySpdfInRFAll = [];
bsPreArraySpdfExRFAll = [];
nsPreArraySpdfAIAll = [];
bsPreArraySpdfAIAll = [];
nsPostArraySpdfInRFAll = [];
nsPostArraySpdfExRFAll = [];
bsPostArraySpdfInRFAll = [];
bsPostArraySpdfExRFAll = [];
nsPostArraySpdfAIAll = [];
bsPostArraySpdfAIAll = [];
nsWfAll = [];
bsWfAll = [];
for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;

    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    if ~isempty(spikeTimes)
        if all(spikeStruct.peakSmoothedAmps > 0) && ...
                (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
                (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
                spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
                spikeStruct.peakSmoothedAmps(1) < 1/4 * max(spikeStruct.peakSmoothedAmps))) && ...
                any(spikeStruct.peakSmoothedAmps > -1/4 * spikeStruct.troughSmoothedAmps(1))
            preArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), preArrayOnsetSpdfWindowOffset, kernelSigma);
            preArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfInRF = edpsth_notranspose(preArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            preArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfExRF = edpsth_notranspose(preArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            
            postArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), postArrayOnsetSpdfWindowOffset, kernelSigma);
            postArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, periArrayOnsetWindow);
            postArrayOnsetSpdfInRF = edpsth_notranspose(postArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, postArrayOnsetT);
            postArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, periArrayOnsetWindow);
            postArrayOnsetSpdfExRF = edpsth_notranspose(postArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, postArrayOnsetT);

            targetDimT = computeTForSpdf(periTargetDimWindow(1), targetDimSpdfWindowOffset, kernelSigma);
            targetDimSpikeTimesInRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{inRFLoc}, periTargetDimWindow);
            targetDimSpdfInRF = edpsth_notranspose(targetDimSpikeTimesInRF, kernelSigma, 'n', [], 0, targetDimT);
            targetDimSpikeTimesExRF = createnonemptydatamatpt(spikeVar, targetDimByLoc{exRFLoc}, periTargetDimWindow);
            targetDimSpdfExRF = edpsth_notranspose(targetDimSpikeTimesExRF, kernelSigma, 'n', [], 0, targetDimT);

            % attention index, mean over time
            cueTargetDelayAI = (mean(preArrayOnsetSpdfInRF) - mean(preArrayOnsetSpdfExRF)) / (mean(preArrayOnsetSpdfInRF) + mean(preArrayOnsetSpdfExRF));
            arrayResponseAI = (mean(postArrayOnsetSpdfInRF) - mean(postArrayOnsetSpdfExRF)) / (mean(postArrayOnsetSpdfInRF) + mean(postArrayOnsetSpdfExRF));
            targetDimDelayAI = (mean(targetDimSpdfInRF) - mean(targetDimSpdfExRF)) / (mean(targetDimSpdfInRF) + mean(targetDimSpdfExRF));

            if any(isnan(preArrayOnsetSpdfInRF)) || any(isnan(preArrayOnsetSpdfExRF))
                continue;
            end
            nWfTime = numel(D.allSpikeStructs{i}.meanWf);
            spikeFs = 40000;
            t = (1:nWfTime)/(spikeFs/1000);
%             tShift = spikeStruct.firstTroughIndex/(spikeFs/1000);
            [~, troughIndex] = min(spikeStruct.meanWf);
            tShift = troughIndex/(spikeFs/1000);
            minFR = 2;
            if mean(preArrayOnsetSpdfInRF) >= minFR || mean(preArrayOnsetSpdfExRF) >= minFR
                if spikeStruct.troughToPeakTime <= fsClassThresh
                    nsWfAll = [nsWfAll; spikeStruct.meanWf / max(spikeStruct.meanWf)];
                    nsPreArraySpdfInRFAll = [nsPreArraySpdfInRFAll; preArrayOnsetSpdfInRF];
                    nsPreArraySpdfExRFAll = [nsPreArraySpdfExRFAll; preArrayOnsetSpdfExRF];
                    nsPreArraySpdfAIAll = [nsPreArraySpdfAIAll; (preArrayOnsetSpdfInRF - preArrayOnsetSpdfExRF) ./ (preArrayOnsetSpdfInRF + preArrayOnsetSpdfExRF)];
                    nsPostArraySpdfInRFAll = [nsPostArraySpdfInRFAll; postArrayOnsetSpdfInRF];
                    nsPostArraySpdfExRFAll = [nsPostArraySpdfExRFAll; postArrayOnsetSpdfExRF];
                    nsPostArraySpdfAIAll = [nsPostArraySpdfAIAll; (postArrayOnsetSpdfInRF - postArrayOnsetSpdfExRF) ./ (postArrayOnsetSpdfInRF + postArrayOnsetSpdfExRF)];
                    % attention index
                    nsCueTargetDelayActivity(i) = cueTargetDelayAI;
                    nsArrayResponseActivity(i) = arrayResponseAI;
                    nsTargetDimDelayActivity(i) = targetDimDelayAI;
                    fprintf('Unit %d - NS - AI: %0.1f\n', i, cueTargetDelayAI);
                    unitCategory(i) = 6;
                    if cueTargetDelayAI > 0%0.05
                        unitCategory(i) = 7;
                    end
                    if cueTargetDelayAI < 0%-0.05
                        unitCategory(i) = 5;
                    end
                else
                    bsWfAll = [bsWfAll; spikeStruct.meanWf / max(spikeStruct.meanWf)];
                    bsPreArraySpdfInRFAll = [bsPreArraySpdfInRFAll; preArrayOnsetSpdfInRF];
                    bsPreArraySpdfExRFAll = [bsPreArraySpdfExRFAll; preArrayOnsetSpdfExRF];
                    bsPreArraySpdfAIAll = [bsPreArraySpdfAIAll; (preArrayOnsetSpdfInRF - preArrayOnsetSpdfExRF) ./ (preArrayOnsetSpdfInRF + preArrayOnsetSpdfExRF)];
                    bsPostArraySpdfInRFAll = [bsPostArraySpdfInRFAll; postArrayOnsetSpdfInRF];
                    bsPostArraySpdfExRFAll = [bsPostArraySpdfExRFAll; postArrayOnsetSpdfExRF];
                    bsPostArraySpdfAIAll = [bsPostArraySpdfAIAll; (postArrayOnsetSpdfInRF - postArrayOnsetSpdfExRF) ./ (postArrayOnsetSpdfInRF + postArrayOnsetSpdfExRF)];
                    % attention index
                    bsCueTargetDelayActivity(i) = cueTargetDelayAI;
                    bsArrayResponseActivity(i) = arrayResponseAI;
                    bsTargetDimDelayActivity(i) = targetDimDelayAI;
                    fprintf('Unit %d - BS - AI: %0.1f\n', i, cueTargetDelayAI);
                    unitCategory(i) = 3;
                    if cueTargetDelayAI > 0%0.05
                        unitCategory(i) = 4;
                    end
                    if cueTargetDelayAI < 0%-0.05
                        unitCategory(i) = 2;
                    end
                end
            end
        end
    end
end

%%
allChannelIDs = cellfun(@(x) x.channelID, D.allSpikeStructs);    
bsNsHistCounts = nan(5, 64);
for i = 1:7
    bsNsHistCounts(i,:) = histcounts(allChannelIDs(unitCategory == i), 1:65);
end
figure_tr_inch(6, 10);
bh = barh(bsNsHistCounts', 'stacked');
bh(1).FaceColor = 0.7*ones(3,1);
bh(2).FaceColor = [255 200 100]/255;
bh(3).FaceColor = [200 150 50]/255;
bh(4).FaceColor = [150 100 0]/255;
bh(5).FaceColor = [255 50 255]/255;
bh(6).FaceColor = [175 25 200]/255;
bh(7).FaceColor = [100 0 150]/255;

legend({'Axon? / Low FR', 'BS AI<0.05', 'BS', 'BS AI>0.05', 'NS AI<0.05', 'NS', 'NS AI>0.05'}, 'Location', 'NorthEastOutside');
set(gca, 'YDir', 'reverse');
ylim([0 65]);
origXLim = xlim();
hold on;
plot(origXLim, [32.5 32.5], 'Color', 0.3*ones(3, 1));
xlim(origXLim);
xlabel('Number of Units');
ylabel('Channel Number (1 = top probe 1, 33 = top probe 2)');
title('Delay Period AI by Cell Type, Channel Number');
set(gca, 'FontSize', 12);

plotFileName = sprintf('%s/delayPeriodHistAIByCellTypeChannel.png', processedDataDir);
export_fig(plotFileName, '-nocrop');


%%
save('for-grant-fig9-temp4.mat', 'bsNsHistCounts');

%%
a1 = load('for-grant-fig9-temp1.mat');
a2 = load('for-grant-fig9-temp2.mat');

%%
a1 = load('for-grant-fig9-temp3.mat');
a2 = load('for-grant-fig9-temp4.mat');

%%
% 20170327 Vprobe2 contact 32 corresponds to Vprobe1 contact 31
allCounts = zeros(7, 33);
allCounts(:,2:33) = a2.bsNsHistCounts(:,1:32); % offset 1
allCounts(:,1:32) = allCounts(:,1:32) + a2.bsNsHistCounts(:,33:64);
allCounts(:,2:33) = allCounts(:,2:33) + a1.bsNsHistCounts(:,1:32); % offset 1
% channels 30-33 PI

figure_tr_inch(2.5, 2);
subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.13, 'MB', 0.15, 'MR', 0.33);
set(gcf, 'Color', 'white');
bh = barh(allCounts(2:7,:)', 'stacked');
cwMap = getCoolWarmMap();
cols = cwMap(round(linspace(1, size(cwMap, 1), 7)),:);
bh(1).FaceColor = cols(3,:);%[255 200 100]/255;
bh(2).FaceColor = cols(2,:);%[200 150 50]/255;
bh(3).FaceColor = cols(1,:);%[150 100 0]/255;
bh(4).FaceColor = cols(5,:);%[255 50 255]/255;
bh(5).FaceColor = cols(6,:);%[175 25 200]/255;
bh(6).FaceColor = cols(7,:);%[100 0 150]/255;
box off;

hold on;
plot([-0.6 4], [29.5 29.5], 'Color', 0.2*ones(3, 1));
plot([-0.6 4], [22.5 22.5], 'Color', 0.2*ones(3, 1));
plot([-0.6 4], [0.5 0.5], 'Color', 0.2*ones(3, 1));

lh = legend({'BS AI<0.05', 'BS', 'BS AI>0.05', 'NS AI<0.05', 'NS', 'NS AI>0.05'}, 'Location', 'SouthEast', 'FontSize', 8);
currPos = get(lh, 'Position');
set(lh, 'Position', [0.58 0.23 currPos(3) currPos(4)]);
set(gca, 'YDir', 'reverse');
ylim([0.5 33.5]);
xlim([-0.6 4]);
hold on;
xlabel('Number of Units');
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'FontSize', 10);
set(gca, 'Color', 'none');
set(gca, 'YColor', 'none');
yl = ylabel({'Depth (0.15mm bins)'}, 'Color', 'k');
set(yl, 'Visible', 'on');
ylPos = get(yl, 'Position');
set(yl, 'Position', ylPos + [-0.3 0 0]);
text(-0.1, 31.25, 'PI', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(-0.1, 26, 'bSC', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(-0.1, 11, 'PLd/', 'FontSize', 9, 'HorizontalAlignment', 'right');
text(-0.1, 14, 'PM', 'FontSize', 9, 'HorizontalAlignment', 'right');

plotFileName = sprintf('grant-fig-delaySpiking-byDepth.png');
export_fig(plotFileName, '-nocrop');

%%
nsCueTargetDelayActivity(isnan(nsCueTargetDelayActivity)) = [];
bsCueTargetDelayActivity(isnan(bsCueTargetDelayActivity)) = [];
nsArrayResponseActivity(isnan(nsArrayResponseActivity)) = [];
bsArrayResponseActivity(isnan(bsArrayResponseActivity)) = [];
nsTargetDimDelayActivity(isnan(nsTargetDimDelayActivity)) = [];
bsTargetDimDelayActivity(isnan(bsTargetDimDelayActivity)) = [];

%%
outlierMaxSDs = 2;
nsOutlier = abs((nsCueTargetDelayActivity - mean(nsCueTargetDelayActivity)) / std(nsCueTargetDelayActivity)) > outlierMaxSDs | ...
        abs((nsArrayResponseActivity - mean(nsArrayResponseActivity)) / std(nsArrayResponseActivity)) > outlierMaxSDs;
bsOutlier = abs((bsCueTargetDelayActivity - mean(bsCueTargetDelayActivity)) / std(bsCueTargetDelayActivity)) > outlierMaxSDs | ...
        abs((bsArrayResponseActivity - mean(bsArrayResponseActivity)) / std(bsArrayResponseActivity)) > outlierMaxSDs;
fprintf('%d NS outliers\n', sum(nsOutlier));
fprintf('%d BS outliers\n', sum(bsOutlier));

nsCueTargetDelayActivity(nsOutlier) = [];
bsCueTargetDelayActivity(bsOutlier) = [];
nsArrayResponseActivity(nsOutlier) = [];
bsArrayResponseActivity(bsOutlier) = [];
nsPreArraySpdfInRFAll(nsOutlier,:) = [];
nsPreArraySpdfExRFAll(nsOutlier,:) = [];
nsPreArraySpdfAIAll(nsOutlier,:) = [];
bsPreArraySpdfInRFAll(bsOutlier,:) = [];
bsPreArraySpdfExRFAll(bsOutlier,:) = [];
bsPreArraySpdfAIAll(bsOutlier,:) = [];
nsPostArraySpdfInRFAll(nsOutlier,:) = [];
nsPostArraySpdfExRFAll(nsOutlier,:) = [];
nsPostArraySpdfAIAll(nsOutlier,:) = [];
bsPostArraySpdfInRFAll(bsOutlier,:) = [];
bsPostArraySpdfExRFAll(bsOutlier,:) = [];
bsPostArraySpdfAIAll(bsOutlier,:) = [];

%%
figure_tr_inch(4.5, 2);
set(gcf, 'Color', 'white');

sa1 = subaxis(1, 2, 1, 'SH', 0.13, 'ML', 0.1, 'MT', 0.05, 'MB', 0.25);
% hold on;
% box off;
% plot(sa1, t - tShift, nsWfAll, 'Color', [234 130 127]/255, 'LineWidth', 2);
% plot(sa1, t - tShift, mean(bsWfAll), 'Color', [108 175 222]/255, 'LineWidth', 2);
% ylim([-3 1]);
% xlim([-0.2 0.6]);
% xlabel('Time from Trough (ms)');
% ylabel('Normalized Amplitude');
% title('Spike Waveforms');
hold on;
box off;
arrayOnsetPlotT = preArrayOnsetT - periArrayOnsetWindow(1);
bsInRFSE = std(bsPreArraySpdfInRFAll)/sqrt(size(bsPreArraySpdfInRFAll, 1));
bsExRFSE = std(bsPreArraySpdfExRFAll)/sqrt(size(bsPreArraySpdfExRFAll, 1));
jbfill(arrayOnsetPlotT, mean(bsPreArraySpdfInRFAll) + bsInRFSE, mean(bsPreArraySpdfInRFAll) - bsInRFSE, [85 90 152]/255, [85 90 152]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(bsPreArraySpdfExRFAll) + bsExRFSE, mean(bsPreArraySpdfExRFAll) - bsExRFSE, [108 175 222]/255, [108 175 222]/255, 1, 0.3);
hold on;
plot(sa1, arrayOnsetPlotT, mean(bsPreArraySpdfInRFAll), 'Color', [85 90 152]/255, 'LineWidth', 2);
plot(sa1, arrayOnsetPlotT, mean(bsPreArraySpdfExRFAll), 'Color', [108 175 222]/255, 'LineWidth', 2);
arrayOnsetPlotT = postArrayOnsetT - periArrayOnsetWindow(1);
bsInRFSE = std(bsPostArraySpdfInRFAll)/sqrt(size(bsPostArraySpdfInRFAll, 1));
bsExRFSE = std(bsPostArraySpdfExRFAll)/sqrt(size(bsPostArraySpdfExRFAll, 1));
jbfill(arrayOnsetPlotT, mean(bsPostArraySpdfInRFAll) + bsInRFSE, mean(bsPostArraySpdfInRFAll) - bsInRFSE, [85 90 152]/255, [85 90 152]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(bsPostArraySpdfExRFAll) + bsExRFSE, mean(bsPostArraySpdfExRFAll) - bsExRFSE, [108 175 222]/255, [108 175 222]/255, 1, 0.3);
hold on;
plot(sa1, arrayOnsetPlotT, mean(bsPostArraySpdfInRFAll), 'Color', [85 90 152]/255, 'LineWidth', 2);
plot(sa1, arrayOnsetPlotT, mean(bsPostArraySpdfExRFAll), 'Color', [108 175 222]/255, 'LineWidth', 2);
plot(sa1, [0 0], [0 100], 'Color', 0.3*ones(3, 1));
ylim([0 25]);
xlim([-0.4 0.3]);
xlabel('Time from Array Onset (s)');
ylabel('Firing Rate (Hz)');
textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
text(0.02, 1.0, 'BS: Attend-RF', 'Color', [85 90 152]/255, textParams{:});
text(0.02, 0.9, 'BS: Attend-Away', 'Color', [108 175 222]/255, textParams{:});

sa2 = subaxis(1, 2, 2);
hold on;
box off;
arrayOnsetPlotT = preArrayOnsetT - periArrayOnsetWindow(1);
nsInRFSE = std(nsPreArraySpdfInRFAll)/sqrt(size(nsPreArraySpdfInRFAll, 1));
nsExRFSE = std(nsPreArraySpdfExRFAll)/sqrt(size(nsPreArraySpdfExRFAll, 1));
jbfill(arrayOnsetPlotT, mean(nsPreArraySpdfInRFAll) + nsInRFSE, mean(nsPreArraySpdfInRFAll) - nsInRFSE, [136 85 83]/255, [136 85 83]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(nsPreArraySpdfExRFAll) + nsExRFSE, mean(nsPreArraySpdfExRFAll) - nsExRFSE, [234 130 127]/255, [234 130 127]/255, 1, 0.3);
hold on;
plot(sa2, arrayOnsetPlotT, mean(nsPreArraySpdfInRFAll), 'Color', [136 85 83]/255, 'LineWidth', 2);
plot(sa2, arrayOnsetPlotT, mean(nsPreArraySpdfExRFAll), 'Color', [234 130 127]/255, 'LineWidth', 2);
arrayOnsetPlotT = postArrayOnsetT - periArrayOnsetWindow(1);
nsInRFSE = std(nsPostArraySpdfInRFAll)/sqrt(size(nsPostArraySpdfInRFAll, 1));
nsExRFSE = std(nsPostArraySpdfExRFAll)/sqrt(size(nsPostArraySpdfExRFAll, 1));
jbfill(arrayOnsetPlotT, mean(nsPostArraySpdfInRFAll) + nsInRFSE, mean(nsPostArraySpdfInRFAll) - nsInRFSE, [136 85 83]/255, [136 85 83]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(nsPostArraySpdfExRFAll) + nsExRFSE, mean(nsPostArraySpdfExRFAll) - nsExRFSE, [234 130 127]/255, [234 130 127]/255, 1, 0.3);
hold on;
plot(sa2, arrayOnsetPlotT, mean(nsPostArraySpdfInRFAll), 'Color', [136 85 83]/255, 'LineWidth', 2);
plot(sa2, arrayOnsetPlotT, mean(nsPostArraySpdfExRFAll), 'Color', [234 130 127]/255, 'LineWidth', 2);
plot(sa2, [0 0], [0 100], 'Color', 0.3*ones(3, 1));
ylim([0 25]);
xlim([-0.4 0.3]);
xlabel('Time from Array Onset (s)');
ylabel('Firing Rate (Hz)');
textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
text(0.02, 1.0, 'NS: Attend-RF', 'Color', [136 85 83]/255, textParams{:});
text(0.02, 0.9, 'NS: Attend-Away', 'Color', [234 130 127]/255, textParams{:});

plotFileName = sprintf('%s/grantFig-attentionEffectOnFR.png', processedDataDir);
export_fig(plotFileName, '-nocrop');

figure_tr_inch(2.5, 2);
set(gcf, 'Color', 'white');
subaxis(1, 1, 1, 'ML', 0.22, 'MT', 0.05, 'MB', 0.25);
hold on;
box off;
arrayOnsetPlotT = preArrayOnsetT - periArrayOnsetWindow(1);
bsAISE = std(bsPreArraySpdfAIAll)/sqrt(size(bsPreArraySpdfAIAll, 1));
nsAISE = std(nsPreArraySpdfAIAll)/sqrt(size(nsPreArraySpdfAIAll, 1));
jbfill(arrayOnsetPlotT, mean(bsPreArraySpdfAIAll) + bsAISE, mean(bsPreArraySpdfAIAll) - bsAISE, [108 175 222]/255, [108 175 222]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(nsPreArraySpdfAIAll) + nsAISE, mean(nsPreArraySpdfAIAll) - nsAISE, [234 130 127]/255, [234 130 127]/255, 1, 0.3);
hold on;
plot(arrayOnsetPlotT, mean(nsPreArraySpdfAIAll), 'Color', [234 130 127]/255, 'LineWidth', 2);
plot(arrayOnsetPlotT, mean(bsPreArraySpdfAIAll), 'Color', [108 175 222]/255, 'LineWidth', 2);
arrayOnsetPlotT = postArrayOnsetT - periArrayOnsetWindow(1);
bsAISE = std(bsPostArraySpdfAIAll)/sqrt(size(bsPostArraySpdfAIAll, 1));
nsAISE = std(nsPostArraySpdfAIAll)/sqrt(size(nsPostArraySpdfAIAll, 1));
jbfill(arrayOnsetPlotT, mean(bsPostArraySpdfAIAll) + bsAISE, mean(bsPostArraySpdfAIAll) - bsAISE, [108 175 222]/255, [108 175 222]/255, 1, 0.3);
hold on;
jbfill(arrayOnsetPlotT, mean(nsPostArraySpdfAIAll) + nsAISE, mean(nsPostArraySpdfAIAll) - nsAISE, [234 130 127]/255, [234 130 127]/255, 1, 0.3);
hold on;
plot(arrayOnsetPlotT, mean(nsPostArraySpdfAIAll), 'Color', [234 130 127]/255, 'LineWidth', 2);
plot(arrayOnsetPlotT, mean(bsPostArraySpdfAIAll), 'Color', [108 175 222]/255, 'LineWidth', 2);
plot([0 0], [-100 100], 'Color', 0.3*ones(3, 1));
plot(arrayOnsetPlotT, zeros(size(arrayOnsetPlotT)), 'Color', 0.3*ones(3, 1));
ylim([-0.1 0.4]);
xlim([-0.4 0.3]);
xlabel('Time from Array Onset (s)');
ylabel('Attention Index');
textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
text(0.02, 0.95, sprintf('BS N=%d', size(bsPreArraySpdfInRFAll, 1)), 'Color', [108 175 222]/255, textParams{:});
text(0.02, 0.85, sprintf('NS N=%d', size(nsPreArraySpdfInRFAll, 1)), 'Color', [234 130 127]/255, textParams{:});

plotFileName = sprintf('%s/grantFig-attentionEffectOnAI.png', processedDataDir);
export_fig(plotFileName, '-nocrop');



%%
figure_tr_inch(3, 2);
set(gcf, 'Color', 'white');
sa1 = subaxis(1, 1, 1, 'ML', 0.15, 'MT', 0.05, 'MB', 0.25);
hold on;
box off;
histogram(bsCueTargetDelayActivity, -0.5:0.1:0.5);
histogram(nsCueTargetDelayActivity, -0.5:0.1:0.5);
origYLim = ylim();
plot([0 0], origYLim, 'Color', 0.3*ones(3, 1));
ylim(origYLim);
legend({'BS', 'NS'});
xlabel('Attention Index');
ylabel('Number of Units');

plotFileName = sprintf('%s/grantFig-attentionEffectOnAIHistogram.png', processedDataDir);
export_fig(plotFileName, '-nocrop');

%%
figure_tr_inch(3, 2);
set(gcf, 'Color', 'white');
sa1 = subaxis(1, 1, 1, 'ML', 0.15, 'MT', 0.05, 'MB', 0.25);
hold on;
box off;
histogram(bsArrayResponseActivity, -0.5:0.1:0.5);
histogram(nsArrayResponseActivity, -0.5:0.1:0.5);
origYLim = ylim();
plot([0 0], origYLim, 'Color', 0.3*ones(3, 1));
ylim(origYLim);
legend({'BS', 'NS'});
xlabel('Attention Index');
ylabel('Number of Units');

plotFileName = sprintf('%s/grantFig-attentionEffectOnAIHistogramArray.png', processedDataDir);
export_fig(plotFileName, '-nocrop');

%%
figure_tr_inch(3, 2);
set(gcf, 'Color', 'white');
sa1 = subaxis(1, 1, 1, 'ML', 0.2, 'MT', 0.05, 'MB', 0.25);
hold on;
box off;
means = [mean(bsCueTargetDelayActivity) mean(nsCueTargetDelayActivity); ...
        mean(bsArrayResponseActivity) mean(nsArrayResponseActivity)];
ses = [std(bsCueTargetDelayActivity)/sqrt(size(bsCueTargetDelayActivity, 1)) ...
        std(nsCueTargetDelayActivity)/sqrt(size(nsCueTargetDelayActivity, 1)); ...
        std(bsArrayResponseActivity)/sqrt(size(bsArrayResponseActivity, 1)) ...
        std(nsArrayResponseActivity)/sqrt(size(nsArrayResponseActivity, 1))];
% b = bar(means);
b = barwitherr(ses, means);
b(1).FaceColor = [108 175 222]/255;
b(2).FaceColor = [234 130 127]/255;

text(0.02, 0.95, sprintf('BS N=%d', size(bsCueTargetDelayActivity, 1)), 'Color', [108 175 222]/255, textParams{:});
text(0.02, 0.85, sprintf('NS N=%d', size(nsCueTargetDelayActivity, 1)), 'Color', [234 130 127]/255, textParams{:});
% histogram(bsCueTargetDelayActivity, -0.5:0.1:0.5);
% histogram(nsCueTargetDelayActivity, -0.5:0.1:0.5);
% origYLim = ylim();
% plot([0 0], origYLim, 'Color', 0.3*ones(3, 1));
% ylim(origYLim);
ylabel('Attention Index');
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel', {'Delay', 'Array'});
ylim([-0.05 0.15]);


plotFileName = sprintf('%s/grantFig-attentionEffectOnAIBar.png', processedDataDir);
export_fig(plotFileName, '-nocrop');

%% for grant: noise correlations

isCell = zeros(nUnits, 1);
minFR = 5;

kernelSigma = 0.01;
fsClassThresh = 0.35/1000;

periArrayOnsetWindow = [0.7 0.5];
preArrayOnsetSpdfWindowOffset = [-0.5 0];
postArrayOnsetSpdfWindowOffset = [0 0.3];
preArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), preArrayOnsetSpdfWindowOffset, kernelSigma);
postArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), postArrayOnsetSpdfWindowOffset, kernelSigma);

for i = 1:nUnits
    spikeStruct = D.allSpikeStructs{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;

    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nUnits, round(i/nUnits*100));
    
    if ~isempty(spikeTimes)
        if all(spikeStruct.peakSmoothedAmps > 0) && ...
                (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
                (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
                spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
                spikeStruct.peakSmoothedAmps(1) < 1/4 * max(spikeStruct.peakSmoothedAmps))) && ...
                any(spikeStruct.peakSmoothedAmps > -1/4 * spikeStruct.troughSmoothedAmps(1))
            
            preArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfInRF = edpsth_notranspose(preArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            preArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfExRF = edpsth_notranspose(preArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            
            postArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, periArrayOnsetWindow);
            postArrayOnsetSpdfInRF = edpsth_notranspose(postArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, postArrayOnsetT);
            postArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, periArrayOnsetWindow);
            postArrayOnsetSpdfExRF = edpsth_notranspose(postArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, postArrayOnsetT);
            
            if mean(preArrayOnsetSpdfInRF) >= minFR && mean(preArrayOnsetSpdfExRF) >= minFR && ...
                    mean(postArrayOnsetSpdfInRF) >= minFR && mean(postArrayOnsetSpdfExRF) >= minFR
                isCell(i) = 1;
            end
        end
    end
end

%% create binned binary spike variable
spikeStructCells = D.allSpikeStructs;%D.allSpikeStructs(isCell == 1);
nCells = numel(spikeStructCells);
nUnits = numel(D.allSpikeStructs);
origSpikeFs = D.allSpikeStructs{1}.Fs;
spikeFs = 200; % 1000 Hz, period 1 ms; 200 Hz, period 5 ms

preArrayOnsetSpikeTimesWindow = [0.4 0];
postArrayOnsetSpikeTimesWindow = [0 0.3];

preArrayOnsetSpikeTimesInRF = cell(nCells, 1);
preArrayOnsetSpikeTimesExRF = cell(nCells, 1);
postArrayOnsetSpikeTimesInRF = cell(nCells, 1);
postArrayOnsetSpikeTimesExRF = cell(nCells, 1);

preArrayOnsetSpikeTimesInRFBinary = false(nCells, numel(arrayOnsetByLoc{inRFLoc}), round(sum(preArrayOnsetSpikeTimesWindow)/spikeFs));
preArrayOnsetSpikeTimesExRFBinary = false(nCells, numel(arrayOnsetByLoc{exRFLoc}), round(sum(preArrayOnsetSpikeTimesWindow)/spikeFs));
postArrayOnsetSpikeTimesInRFBinary = false(nCells, numel(arrayOnsetByLoc{inRFLoc}), round(sum(postArrayOnsetSpikeTimesWindow)/spikeFs));
postArrayOnsetSpikeTimesExRFBinary = false(nCells, numel(arrayOnsetByLoc{exRFLoc}), round(sum(postArrayOnsetSpikeTimesWindow)/spikeFs));

preArrayOnsetSpikeCountInRF = nan(nCells, numel(arrayOnsetByLoc{inRFLoc}));
postArrayOnsetSpikeCountInRF = nan(nCells, numel(arrayOnsetByLoc{inRFLoc}));
preArrayOnsetSpikeCountExRF = nan(nCells, numel(arrayOnsetByLoc{exRFLoc}));
postArrayOnsetSpikeCountExRF = nan(nCells, numel(arrayOnsetByLoc{exRFLoc}));

for i = 1:nCells
    spikeStruct = spikeStructCells{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;
    
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nCells, round(i/nCells*100));
    
    preArrayOnsetSpikeTimesInRF{i} = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, preArrayOnsetSpikeTimesWindow);
    preArrayOnsetSpikeTimesExRF{i} = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, preArrayOnsetSpikeTimesWindow);

    postArrayOnsetSpikeTimesInRF{i} = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{inRFLoc}, postArrayOnsetSpikeTimesWindow);
    postArrayOnsetSpikeTimesExRF{i} = createnonemptydatamatpt(spikeVar, arrayOnsetByLoc{exRFLoc}, postArrayOnsetSpikeTimesWindow);
    
    for j = 1:numel(arrayOnsetByLoc{inRFLoc})
        preArrayOnsetSpikeTimesInRFBinary(i,j,(ceil(spikeFs*preArrayOnsetSpikeTimesInRF{i}(j).times))) = true;
        postArrayOnsetSpikeTimesInRFBinary(i,j,(ceil(spikeFs*postArrayOnsetSpikeTimesInRF{i}(j).times))) = true;
        preArrayOnsetSpikeCountInRF(i,j) = numel(preArrayOnsetSpikeTimesInRF{i}(j).times);
        postArrayOnsetSpikeCountInRF(i,j) = numel(postArrayOnsetSpikeTimesInRF{i}(j).times);
    end
    
    for j = 1:numel(arrayOnsetByLoc{exRFLoc})
        preArrayOnsetSpikeTimesExRFBinary(i,j,(ceil(spikeFs*preArrayOnsetSpikeTimesExRF{i}(j).times))) = true;
        postArrayOnsetSpikeTimesExRFBinary(i,j,(ceil(spikeFs*postArrayOnsetSpikeTimesExRF{i}(j).times))) = true;
        preArrayOnsetSpikeCountExRF(i,j) = numel(preArrayOnsetSpikeTimesExRF{i}(j).times);
        postArrayOnsetSpikeCountExRF(i,j) = numel(postArrayOnsetSpikeTimesExRF{i}(j).times);
    end
end

noSpikingCells = all(preArrayOnsetSpikeCountInRF == 0, 2) | all(preArrayOnsetSpikeCountExRF == 0, 2) | ...
        all(postArrayOnsetSpikeCountInRF == 0, 2) | all(postArrayOnsetSpikeCountExRF == 0, 2) ;
% 
% preArrayOnsetSpikeCountInRF(noSpikingCells,:) = [];
% preArrayOnsetSpikeCountExRF(noSpikingCells,:) = [];
% postArrayOnsetSpikeCountInRF(noSpikingCells,:) = [];
% postArrayOnsetSpikeCountExRF(noSpikingCells,:) = [];

fprintf('%d remaining cells\n', size(preArrayOnsetSpikeCountInRF, 1));

%% option 1: correlate spike count responses in window
% preArrayOnsetSpikeCountCorrInRF = nan(nCells, nCells);
% preArrayOnsetSpikeCountCorrExRF = nan(nCells, nCells);

preArrayOnsetSpikeCountCorrInRF = corr(preArrayOnsetSpikeCountInRF');
triuPreArrayOnsetSpikeCountCorrInRF = preArrayOnsetSpikeCountCorrInRF(~tril(ones(size(preArrayOnsetSpikeCountCorrInRF))));
preArrayOnsetSpikeCountCorrInRFAll = triuPreArrayOnsetSpikeCountCorrInRF(:);
meanPreArrayOnsetSpikeCountCorrInRF = mean(preArrayOnsetSpikeCountCorrInRFAll)

preArrayOnsetSpikeCountCorrExRF = corr(preArrayOnsetSpikeCountExRF');
triuPreArrayOnsetSpikeCountCorrExRF = preArrayOnsetSpikeCountCorrExRF(~tril(ones(size(preArrayOnsetSpikeCountCorrExRF))));
preArrayOnsetSpikeCountCorrExRFAll = triuPreArrayOnsetSpikeCountCorrExRF(:);
meanPreArrayOnsetSpikeCountCorrExRF = mean(preArrayOnsetSpikeCountCorrExRFAll)

postArrayOnsetSpikeCountCorrInRF = corr(postArrayOnsetSpikeCountInRF');
triuPostArrayOnsetSpikeCountCorrInRF = postArrayOnsetSpikeCountCorrInRF(~tril(ones(size(postArrayOnsetSpikeCountCorrInRF))));
postArrayOnsetSpikeCountCorrInRFAll = triuPostArrayOnsetSpikeCountCorrInRF(:);
meanPostArrayOnsetSpikeCountCorrInRF = mean(postArrayOnsetSpikeCountCorrInRFAll)

postArrayOnsetSpikeCountCorrExRF = corr(postArrayOnsetSpikeCountExRF');
triuPostArrayOnsetSpikeCountCorrExRF = postArrayOnsetSpikeCountCorrExRF(~tril(ones(size(postArrayOnsetSpikeCountCorrExRF))));
postArrayOnsetSpikeCountCorrExRFAll = triuPostArrayOnsetSpikeCountCorrExRF(:);
meanPostArrayOnsetSpikeCountCorrExRF = mean(postArrayOnsetSpikeCountCorrExRFAll)

figure;
subaxis(2, 2, 1);
histogram(preArrayOnsetSpikeCountCorrInRFAll, -1:0.1:1);
subaxis(2, 2, 2);
histogram(preArrayOnsetSpikeCountCorrExRFAll, -1:0.1:1);
subaxis(2, 2, 3);
histogram(postArrayOnsetSpikeCountCorrInRFAll, -1:0.1:1);
subaxis(2, 2, 4);
histogram(postArrayOnsetSpikeCountCorrExRFAll, -1:0.1:1);

%%
% crossCorrInRF = nan(nCells, nCells);
% crossCorrExRF = nan(nCells, nCells);
% 
% for i = 1:nCells
%     spikeStruct = spikeStructCells{i};
%     unitName = spikeStruct.name;
%     spikeTimes = spikeStruct.ts;
%     spikeVar = spikeStruct.ts;
% 
%     fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
%             nCells, round(i/nCells*100));
%     
%     for j = (i+1):nCells
%         preArrayOnsetInRFXcorrAll = nan(numel(preArrayOnsetSpikeTimesInRF), 1);
%         postArrayOnsetInRFXcorrAll = nan(numel(postArrayOnsetSpikeTimesInRF), 1);
%         preArrayOnsetExRFXcorrAll = nan(numel(preArrayOnsetSpikeTimesExRF), 1);
%         postArrayOnsetExRFXcorrAll = nan(numel(postArrayOnsetSpikeTimesExRF), 1);
%         for k = 1:numel(preArrayOnsetSpikeTimesInRF)
% %             trialWiseXcorr(k) = xcorr(preArrayOnsetSpikeTimesInRFBinary(i,k,:), preArrayOnsetSpikeTimesInRFBinary(j,k,:), 
%         end
%         
%         crossCorrInRF(j,i) = crossCorrInRF(i,j);
%         crossCorrExRF(j,i) = crossCorrExRF(i,j);
%     end
% end

%% for grant

for k = 1%1:130
inds = [83 82]; % [41 112]
spikeStructCells = D.allSpikeStructs(inds,:); %D.allSpikeStructs(isCell == 1);
nCells = numel(spikeStructCells);

maxLag = 80;
x = (-maxLag:maxLag)/spikeFs;
preArrayOnsetXCorrInRFAll = nan(numel(arrayOnsetByLoc{inRFLoc}), numel(x));
preArrayOnsetXCorrExRFAll = nan(numel(arrayOnsetByLoc{exRFLoc}), numel(x));

for j = 1:numel(arrayOnsetByLoc{inRFLoc})
    preArrayOnsetXCorrInRFAll(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(1),j,:)), ...
            squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(2),j,:)), maxLag);
end
for j = 1:numel(arrayOnsetByLoc{exRFLoc})
    preArrayOnsetXCorrExRFAll(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(1),j,:)), ...
            squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(2),j,:)), maxLag);
end

spikeTimeSeriesAll = false(2, round(D.totalTicks/(origSpikeFs/spikeFs)));
for i = 1:numel(inds)
    spikeTimeSeriesAll(i,(ceil(spikeFs*D.allSpikeStructs{inds(i)}.ts))) = true;
end
xcorrAll = xcorr(spikeTimeSeriesAll(1,:), spikeTimeSeriesAll(2,:), maxLag);

figure_tr_inch(10, 6);
set(gcf, 'Color', 'white');
h1 = subaxis(1, 2, 1);
box off;
hold on;
plot(x, mean(preArrayOnsetXCorrInRFAll), 'r');
plot(x, mean(preArrayOnsetXCorrExRFAll), 'b');
plot(x, movmean(mean(preArrayOnsetXCorrInRFAll), 5), 'r', 'LineWidth', 3);
plot(x, movmean(mean(preArrayOnsetXCorrExRFAll), 5), 'b', 'LineWidth', 3);

subaxis(1, 2, 2);
box off;
hold on;
plot(x, xcorrAll - mean(xcorrAll), 'Color', [0 0.7 0.3], 'LineWidth', 3);

drawnow;
end

% account for increased firing rate -- one shuffle
nShuffles = 50;
preArrayOnsetXCorrInRFAllShuffleAll = nan(nShuffles, numel(x));
preArrayOnsetXCorrExRFAllShuffleAll = nan(nShuffles, numel(x));

for i = 1:nShuffles
    preArrayOnsetXCorrInRFAllShuffle = nan(numel(arrayOnsetByLoc{inRFLoc}), numel(x));
    preArrayOnsetXCorrExRFAllShuffle = nan(numel(arrayOnsetByLoc{exRFLoc}), numel(x));
    permTrialsInRF = randperm(numel(arrayOnsetByLoc{inRFLoc}));
    permTrialsExRF = randperm(numel(arrayOnsetByLoc{exRFLoc}));
    for j = 1:numel(arrayOnsetByLoc{inRFLoc})
        preArrayOnsetXCorrInRFAllShuffle(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(1),j,:)), ...
                squeeze(preArrayOnsetSpikeTimesInRFBinary(inds(2),permTrialsInRF(j),:)), maxLag);
    end
    for j = 1:numel(arrayOnsetByLoc{exRFLoc})
        preArrayOnsetXCorrExRFAllShuffle(j,:) = xcorr(squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(1),j,:)), ...
                squeeze(preArrayOnsetSpikeTimesExRFBinary(inds(2),permTrialsExRF(j),:)), maxLag);
    end
    preArrayOnsetXCorrInRFAllShuffleAll(i,:) = mean(preArrayOnsetXCorrInRFAllShuffle);
    preArrayOnsetXCorrExRFAllShuffleAll(i,:) = mean(preArrayOnsetXCorrExRFAllShuffle);
end
plot(h1, x, movmean(mean(preArrayOnsetXCorrInRFAllShuffleAll), 5), 'm', 'LineWidth', 3);
plot(h1, x, movmean(mean(preArrayOnsetXCorrExRFAllShuffleAll), 5), 'g', 'LineWidth', 3);




%% for grant spike field coherence

spikeStructCells = D.allSpikeStructs;
nCells = numel(spikeStructCells);
nUnits = numel(D.allSpikeStructs);
origSpikeFs = D.allSpikeStructs{1}.Fs;

preArrayOnsetHoldSpikeTimesWindow = [0.3 0];

preArrayOnsetHoldSpikeTimesByLoc = cell(nCells, nLoc);

for i = 1:nCells
    spikeStruct = spikeStructCells{i};
    unitName = spikeStruct.name;
    spikeTimes = spikeStruct.ts;
    spikeVar = spikeStruct.ts;
    
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nCells, round(i/nCells*100));
    
    for k = [inRFLoc exRFLoc]
        preArrayOnsetHoldSpikeTimesByLoc{i,k} = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{k}, preArrayOnsetHoldSpikeTimesWindow);
    end
end

ref = 'CAR';
isVprobe1 = 1;

if strcmp(ref, 'BIP')
    if isVprobe1
        channelsToProcess = 28:31;
    else
        channelsToProcess = 61:63;
    end
else
    if isVprobe1
        channelsToProcess = 28:32;
    else
        channelsToProcess = 61:64;
    end
end

Fs = D.lfpFs;
if isVprobe1
    car = mean(D.adjLfps(1:32,:))';
else
    car = mean(D.adjLfps(33:64,:))';
end
averagePI = mean(D.adjLfps(channelsToProcess,:))' - car;

arrayOnsetLfpWindow = [0.3 0];
arrayOnsetLfpT = -arrayOnsetLfpWindow(1) * Fs : arrayOnsetLfpWindow(2) * Fs - 1;

lfpAroundArrayOnsetHoldByLocAll = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundArrayOnsetHoldByLocAll{i} = nan(numel(arrayOnsetLfpT), numel(arrayOnsetHoldByLoc{i}));
    if strcmp(ref, 'CAR')
        adjLfp = averagePI;
    elseif strcmp(ref, 'BIP')
%             adjLfp = D.adjLfps(channelsToProcess(j)+1,:)' - D.adjLfps(channelsToProcess(j),:)'; % bipolar LFP
    else
%             adjLfp = D.adjLfps(channelsToProcess(j),:)';
    end
    lfpAroundArrayOnsetHoldByLocAll{i}(:,:) = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
end

load('params.mat');
params.fpass = [5 100];
params.tapers = [2 3];
params.pad = 1;

kernelSigma = 0.01;
periArrayOnsetWindow = [0.7 0.5];
preArrayOnsetSpdfWindowOffset = [-0.3 0];
preArrayOnsetT = computeTForSpdf(periArrayOnsetWindow(1), preArrayOnsetSpdfWindowOffset, kernelSigma);
sfcLineArrayOnsetHoldAll = cell(nCells, nLoc);

if isVprobe1
    pmChannelRange = [1 20];
else
    pmChannelRange = [33 53];
end

for i = 1:nCells
    spikeStruct = spikeStructCells{i};
    unitName = spikeStruct.name;
    spikeVar = spikeStruct.ts;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nCells, round(i/nCells*100));
    if ~isempty(spikeTimes) && spikeStruct.channelID >= pmChannelRange(1) && spikeStruct.channelID <= pmChannelRange(2)
        if all(spikeStruct.peakSmoothedAmps > 0) && ...
                (startsWith(spikeStruct.inflectionPattern, 'tp') || ...
                (startsWith(spikeStruct.inflectionPattern, 'ptp') && ...
                spikeStruct.peakSmoothedAmps(1) < -1/2 * spikeStruct.troughSmoothedAmps(1) && ...
                spikeStruct.peakSmoothedAmps(1) < 1/4 * max(spikeStruct.peakSmoothedAmps))) && ...
                any(spikeStruct.peakSmoothedAmps > -1/4 * spikeStruct.troughSmoothedAmps(1))
            preArrayOnsetSpikeTimesInRF = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{inRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfInRF = edpsth_notranspose(preArrayOnsetSpikeTimesInRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            preArrayOnsetSpikeTimesExRF = createnonemptydatamatpt(spikeVar, arrayOnsetHoldByLoc{exRFLoc}, periArrayOnsetWindow);
            preArrayOnsetSpdfExRF = edpsth_notranspose(preArrayOnsetSpikeTimesExRF, kernelSigma, 'n', [], 0, preArrayOnsetT);
            
            minFR = 2;
            if mean(preArrayOnsetSpdfInRF) >= minFR && mean(preArrayOnsetSpdfExRF) >= minFR
                for k = [inRFLoc exRFLoc]
                    [sfcLineArrayOnsetHoldAll{i,k},phi,S12,~,~,fArr3,~,~,~] = coherencycpt(...
                            lfpAroundArrayOnsetHoldByLocAll{k},...
                            preArrayOnsetHoldSpikeTimesByLoc{i,k}, ...
                            params);
                end
            end
        end
    end
end

%%
fsClassThresh = 0.35/1000;

allNSSFCDiff = [];
allBSSFCDiff = [];
allNSSFCInRF = [];
allNSSFCExRF = [];
allBSSFCInRF = [];
allBSSFCExRF = [];

figure;
hold on;
for i = 1:nCells
    spikeStruct = spikeStructCells{i};
    unitName = spikeStruct.name;
    fprintf('Processing %s (%d/%d = %d%%)... \n', unitName, i, ...
            nCells, round(i/nCells*100));
    if ~isempty(sfcLineArrayOnsetHoldAll{i,inRFLoc}) && ~any(isnan(sfcLineArrayOnsetHoldAll{i,inRFLoc})) && ~any(isnan(sfcLineArrayOnsetHoldAll{i,exRFLoc}))
        fprintf('\tPlotting %s (%d/%d = %d%%)... \n', unitName, i, ...
                    nCells, round(i/nCells*100));
        if spikeStruct.troughToPeakTime <= fsClassThresh
%             plot(fArr3, sfcLineArrayOnsetHoldAll{i,inRFLoc} - sfcLineArrayOnsetHoldAll{i,exRFLoc}, 'r');
            allNSSFCDiff = [allNSSFCDiff; (sfcLineArrayOnsetHoldAll{i,inRFLoc} - sfcLineArrayOnsetHoldAll{i,exRFLoc})'];
            allNSSFCInRF = [allNSSFCInRF; sfcLineArrayOnsetHoldAll{i,inRFLoc}'];
            allNSSFCExRF = [allNSSFCExRF; sfcLineArrayOnsetHoldAll{i,exRFLoc}'];
        else
%             plot(fArr3, sfcLineArrayOnsetHoldAll{i,inRFLoc} - sfcLineArrayOnsetHoldAll{i,exRFLoc}, 'b');
            allBSSFCDiff = [allBSSFCDiff; (sfcLineArrayOnsetHoldAll{i,inRFLoc} - sfcLineArrayOnsetHoldAll{i,exRFLoc})'];
            allBSSFCInRF = [allBSSFCInRF; sfcLineArrayOnsetHoldAll{i,inRFLoc}'];
            allBSSFCExRF = [allBSSFCExRF; sfcLineArrayOnsetHoldAll{i,exRFLoc}'];
        end
    end
end

plot(fArr3, mean(allNSSFCDiff), 'm', 'LineWidth', 3);
plot(fArr3, mean(allBSSFCDiff), 'c', 'LineWidth', 3);
plot([0 100], [0 0], 'Color', 0.3*ones(3, 1));

jbfill(fArr3, mean(allNSSFCDiff) + std(allNSSFCDiff)/sqrt(size(allNSSFCDiff, 1)), mean(allNSSFCDiff) - std(allNSSFCDiff)/sqrt(size(allNSSFCDiff, 1)), [1 0 1], [1 0 1], 1, 0.3);
jbfill(fArr3, mean(allBSSFCDiff) + std(allBSSFCDiff)/sqrt(size(allBSSFCDiff, 1)), mean(allBSSFCDiff) - std(allBSSFCDiff)/sqrt(size(allBSSFCDiff, 1)), [0 1 1], [0 1 1], 1, 0.3);

save('for-grant-fig10-temp1-2tapers.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF', 'allBSSFCInRF', 'allBSSFCExRF');

%%
% a1 = load('for-grant-fig10-temp1.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF');
% a2 = load('for-grant-fig10-temp2.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF');
% a3 = load('for-grant-fig10-temp3.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF');

a1 = load('for-grant-fig10-temp1-2tapers.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF', 'allBSSFCInRF', 'allBSSFCExRF');
a2 = load('for-grant-fig10-temp2-2tapers.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF', 'allBSSFCInRF', 'allBSSFCExRF');
a3 = load('for-grant-fig10-temp3-2tapers.mat', 'allNSSFCDiff', 'allBSSFCDiff', 'allNSSFCInRF', 'allNSSFCExRF', 'allBSSFCInRF', 'allBSSFCExRF');

allNSSFCDiffAll = [a1.allNSSFCDiff; a2.allNSSFCDiff; a3.allNSSFCDiff];
allBSSFCDiffAll = [a1.allBSSFCDiff; a2.allBSSFCDiff; a3.allBSSFCDiff];
allNSSFCInRFAll = [a1.allNSSFCInRF; a2.allNSSFCInRF; a3.allNSSFCInRF];
allNSSFCExRFAll = [a1.allNSSFCExRF; a2.allNSSFCExRF; a3.allNSSFCExRF];
allBSSFCInRFAll = [a1.allBSSFCInRF; a2.allBSSFCInRF; a3.allBSSFCInRF];
allBSSFCExRFAll = [a1.allBSSFCExRF; a2.allBSSFCExRF; a3.allBSSFCExRF];

figure;
hold on;
plot(fArr3, mean(allNSSFCDiffAll), 'm', 'LineWidth', 3);
plot(fArr3, mean(allBSSFCDiffAll), 'c', 'LineWidth', 3);
plot([0 100], [0 0], 'Color', 0.3*ones(3, 1));

jbfill(fArr3, mean(allNSSFCDiffAll) + std(allNSSFCDiffAll)/sqrt(size(allNSSFCDiffAll, 1)), mean(allNSSFCDiffAll) - std(allNSSFCDiffAll)/sqrt(size(allNSSFCDiffAll, 1)), [1 0 1], [1 0 1], 1, 0.3);
jbfill(fArr3, mean(allBSSFCDiffAll) + std(allBSSFCDiffAll)/sqrt(size(allBSSFCDiffAll, 1)), mean(allBSSFCDiffAll) - std(allBSSFCDiffAll)/sqrt(size(allBSSFCDiffAll, 1)), [0 1 1], [0 1 1], 1, 0.3);

%%
figure_tr_inch(3, 2);
subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.22, 'MB', 0.26, 'MR', 0.04);
set(gcf, 'Color', 'white');
box off;
hold on;
jbfill(fArr3, mean(allNSSFCInRFAll) + std(allNSSFCInRFAll)/sqrt(size(allNSSFCInRFAll, 1)), mean(allNSSFCInRFAll) - std(allNSSFCInRFAll)/sqrt(size(allNSSFCInRFAll, 1)), [1 0 0], [1 0 0], 1, 0.3);
jbfill(fArr3, mean(allNSSFCExRFAll) + std(allNSSFCExRFAll)/sqrt(size(allNSSFCExRFAll, 1)), mean(allNSSFCExRFAll) - std(allNSSFCExRFAll)/sqrt(size(allNSSFCExRFAll, 1)), [0 0 1], [0 0 1], 1, 0.3);
hold on;
plot(fArr3, mean(allNSSFCInRFAll), 'r', 'LineWidth', 3);
plot(fArr3, mean(allNSSFCExRFAll), 'b', 'LineWidth', 3);
xlim([5 50]);
ylim([0.02 0.072]);
ylabel('Coherence');
xlabel('Frequency (Hz)');
textParams = {'Units', 'normalized'};
text(0.64, 0.95, 'Attend-RF', 'Color', [1 0 0], textParams{:});
text(0.64, 0.85, 'Attend-Away', 'Color', [0 0 1], textParams{:});
set(gca, 'FontSize', 12);

plotFileName = sprintf('for-grant-fig10-2tapers-sfc-NS-to50Hz.png');
export_fig(plotFileName, '-nocrop');

%%
figure_tr_inch(3, 2);
subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.22, 'MB', 0.26, 'MR', 0.04);
set(gcf, 'Color', 'white');
box off;
hold on;
jbfill(fArr3, mean(allBSSFCInRFAll) + std(allBSSFCInRFAll)/sqrt(size(allBSSFCInRFAll, 1)), mean(allBSSFCInRFAll) - std(allBSSFCInRFAll)/sqrt(size(allBSSFCInRFAll, 1)), [1 0 0], [1 0 0], 1, 0.3);
jbfill(fArr3, mean(allBSSFCExRFAll) + std(allBSSFCExRFAll)/sqrt(size(allBSSFCExRFAll, 1)), mean(allBSSFCExRFAll) - std(allBSSFCExRFAll)/sqrt(size(allBSSFCExRFAll, 1)), [0 0 1], [0 0 1], 1, 0.3);
hold on;
plot(fArr3, mean(allBSSFCInRFAll), 'r', 'LineWidth', 3);
plot(fArr3, mean(allBSSFCExRFAll), 'b', 'LineWidth', 3);
xlim([5 50]);
ylim([0.02 0.072]);
ylabel('Coherence');
xlabel('Frequency (Hz)');
textParams = {'Units', 'normalized'};
text(0.64, 0.95, 'Attend-RF', 'Color', [1 0 0], textParams{:});
text(0.64, 0.85, 'Attend-Away', 'Color', [0 0 1], textParams{:});
set(gca, 'FontSize', 12);

plotFileName = sprintf('for-grant-fig10-2tapers-sfc-BS-to50Hz.png');
export_fig(plotFileName, '-nocrop');

xlim([5 100]);

plotFileName = sprintf('for-grant-fig10-2tapers-sfc-BS-to100Hz.png');
export_fig(plotFileName, '-nocrop');

%%
allSFCInRFAll = [allBSSFCInRFAll; allNSSFCInRFAll];
allSFCExRFAll = [allBSSFCExRFAll; allNSSFCExRFAll];

figure_tr_inch(3, 2);
subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.22, 'MB', 0.26, 'MR', 0.04);
set(gcf, 'Color', 'white');
box off;
hold on;
jbfill(fArr3, mean(allSFCInRFAll) + std(allSFCInRFAll)/sqrt(size(allSFCInRFAll, 1)), mean(allSFCInRFAll) - std(allSFCInRFAll)/sqrt(size(allSFCInRFAll, 1)), [1 0 0], [1 0 0], 1, 0.3);
jbfill(fArr3, mean(allSFCExRFAll) + std(allSFCExRFAll)/sqrt(size(allSFCExRFAll, 1)), mean(allSFCExRFAll) - std(allSFCExRFAll)/sqrt(size(allSFCExRFAll, 1)), [0 0 1], [0 0 1], 1, 0.3);
hold on;
plot(fArr3, mean(allSFCInRFAll), 'r', 'LineWidth', 3);
plot(fArr3, mean(allSFCExRFAll), 'b', 'LineWidth', 3);
xlim([5 50]);
ylim([0.025 0.06]);
ylabel('Coherence');
xlabel('Frequency (Hz)');
textParams = {'Units', 'normalized'};
text(0.64, 0.95, 'Attend-RF', 'Color', [1 0 0], textParams{:});
text(0.64, 0.85, 'Attend-Away', 'Color', [0 0 1], textParams{:});
set(gca, 'FontSize', 12);

plotFileName = sprintf('for-grant-fig10-2tapers-sfc-all-to50Hz.png');
export_fig(plotFileName, '-nocrop');

xlim([5 100]);

plotFileName = sprintf('for-grant-fig10-2tapers-sfc-all-to100Hz.png');
export_fig(plotFileName, '-nocrop');

%%
isProbe1 = 0;
ref = 'CAR';
fprintf('REFERENCING SCHEME: %s\n', ref);

if strcmp(ref, 'BIP')
    if isProbe1
        channelsToProcess = 1:31;
    else
        channelsToProcess = 33:63;
    end
else
    if isProbe1
        channelsToProcess = 1:32;
    else
        channelsToProcess = 33:64;
    end
end
Fs = D.lfpFs;

car = mean(D.adjLfps(channelsToProcess,:))';

cueOnsetLfpWindow = [0.5 1];
cueOnsetLfpXLim = [-0.3 0.7];
cueOnsetLfpT = -cueOnsetLfpWindow(1) * Fs : cueOnsetLfpWindow(2) * Fs - 1;
lfpAroundCueOnsetAll = nan(numel(channelsToProcess), numel(cueOnsetLfpT), numel(cueOnset));

arrayOnsetLfpWindow = [1 0.5];
arrayOnsetLfpXLim = [-0.5 0.3];
arrayOnsetLfpT = -arrayOnsetLfpWindow(1) * Fs : arrayOnsetLfpWindow(2) * Fs - 1;
lfpAroundArrayOnsetAll = nan(numel(channelsToProcess), numel(arrayOnsetLfpT), numel(arrayOnset));

targetDimLfpWindow = [1 0.5];
targetDimLfpXLim = [-0.5 0.3];
targetDimLfpT = -targetDimLfpWindow(1) * Fs : targetDimLfpWindow(2) * Fs - 1;
lfpAroundTargetDimAll = nan(numel(channelsToProcess), numel(targetDimLfpT), numel(targetDim));

for j = 1:numel(channelsToProcess)
    if strcmp(ref, 'CAR')
        adjLfp = D.adjLfps(channelsToProcess(j),:)' - car;
    elseif strcmp(ref, 'BIP')
        adjLfp = D.adjLfps(channelsToProcess(j)+1,:)' - D.adjLfps(channelsToProcess(j),:)'; % bipolar LFP
    else
        adjLfp = D.adjLfps(channelsToProcess(j),:)';
    end
    lfpAroundCueOnsetAll(j,:,:) = createdatamatc(adjLfp, cueOnset, Fs, cueOnsetLfpWindow);
    lfpAroundArrayOnsetAll(j,:,:) = createdatamatc(adjLfp, arrayOnset, Fs, arrayOnsetLfpWindow);
    lfpAroundTargetDimAll(j,:,:) = createdatamatc(adjLfp, targetDim, Fs, targetDimLfpWindow);
end
assert(~any(isnan(lfpAroundCueOnsetAll(:))));
assert(~any(isnan(lfpAroundArrayOnsetAll(:))));
assert(~any(isnan(lfpAroundTargetDimAll(:))));

lfpAroundCueOnsetByLocAll = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLocAll = cell(nLoc, 1);
lfpAroundArrayOnsetShortHoldByLocAll = cell(nLoc, 1);
lfpAroundArrayOnsetLongHoldByLocAll = cell(nLoc, 1);
lfpAroundArrayOnsetRelByLocAll = cell(nLoc, 1);
lfpAroundTargetDimByLocAll = cell(nLoc, 1);
for i = 1:nLoc
    lfpAroundCueOnsetByLocAll{i} = nan(numel(channelsToProcess), numel(cueOnsetLfpT), numel(cueOnsetByLoc{i}));
    lfpAroundArrayOnsetHoldByLocAll{i} = nan(numel(channelsToProcess), numel(arrayOnsetLfpT), numel(arrayOnsetHoldByLoc{i}));
    lfpAroundArrayOnsetShortHoldByLocAll{i} = nan(numel(channelsToProcess), numel(arrayOnsetLfpT), numel(arrayOnsetShortHoldByLoc{i}));
    lfpAroundArrayOnsetLongHoldByLocAll{i} = nan(numel(channelsToProcess), numel(arrayOnsetLfpT), numel(arrayOnsetLongHoldByLoc{i}));
    lfpAroundArrayOnsetRelByLocAll{i} = nan(numel(channelsToProcess), numel(arrayOnsetLfpT), numel(arrayOnsetRelByLoc{i}));
    lfpAroundTargetDimByLocAll{i} = nan(numel(channelsToProcess), numel(targetDimLfpT), numel(targetDimByLoc{i}));
    for j = 1:numel(channelsToProcess)
        fprintf('Processing channel %d...\n', channelsToProcess(j));
        if strcmp(ref, 'CAR')
            adjLfp = D.adjLfps(channelsToProcess(j),:)' - car;
        elseif strcmp(ref, 'BIP')
            adjLfp = D.adjLfps(channelsToProcess(j)+1,:)' - D.adjLfps(channelsToProcess(j),:)'; % bipolar LFP
        else
            adjLfp = D.adjLfps(channelsToProcess(j),:)';
        end
        lfpAroundCueOnsetByLocAll{i}(j,:,:) = createdatamatc(adjLfp, cueOnsetByLoc{i}, Fs, cueOnsetLfpWindow);
        lfpAroundArrayOnsetHoldByLocAll{i}(j,:,:) = createdatamatc(adjLfp, arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
        lfpAroundArrayOnsetShortHoldByLocAll{i}(j,:,:) = createdatamatc(adjLfp, arrayOnsetShortHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    %     lfpAroundArrayOnsetShortHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetShortHoldByLoc{i});
        lfpAroundArrayOnsetLongHoldByLocAll{i}(j,:,:) = createdatamatc(adjLfp, arrayOnsetLongHoldByLoc{i}, Fs, arrayOnsetLfpWindow);
    %     lfpAroundArrayOnsetLongHoldByLocFilt{i} = filtfilt(bFirLowPass, 1, lfpAroundArrayOnsetLongHoldByLoc{i});
        lfpAroundArrayOnsetRelByLocAll{i}(j,:,:) = createdatamatc(adjLfp, arrayOnsetRelByLoc{i}, Fs, arrayOnsetLfpWindow);
        lfpAroundTargetDimByLocAll{i}(j,:,:) = createdatamatc(adjLfp, targetDimByLoc{i}, Fs, targetDimLfpWindow);
    end
end

%%
load('params.mat');
params.fpass = [50 100];
params.tapers = [2 3];
params.pad = 1;
cBounds = [-0.15 0.15];

% cohgramSlidingWindow = [0.3 0.05];

cohCueOnsetAll = cell(nLoc, 1);
cohArrayOnsetHoldAll = cell(nLoc, 1);
cohTargetDimAll = cell(nLoc, 1);

cohLinePreCueWindowOffset = [-0.3 0];
preCueLfpSelection = cueOnsetLfpT >= cohLinePreCueWindowOffset(1)*Fs & ...
        cueOnsetLfpT <= cohLinePreCueWindowOffset(2)*Fs;
cohLineCueOnsetAll = cell(numel(channelsToProcess), numel(channelsToProcess), nLoc);

cohLinePreArrayWindowOffset = [-0.3 0];
preArrayLfpSelection = arrayOnsetLfpT >= cohLinePreArrayWindowOffset(1)*Fs & ...
        arrayOnsetLfpT <= cohLinePreArrayWindowOffset(2)*Fs;
cohLineArrayOnsetHoldAll = cell(numel(channelsToProcess), numel(channelsToProcess), nLoc);

for i = 1:numel(channelsToProcess)
    for j = i+1:numel(channelsToProcess) % TODO make loop more efficient
        if i == j
            continue;
        end
        ch1 = channelsToProcess(i);
        ch2 = channelsToProcess(j);
        
        fprintf('Processing coherence between channel %d and channel %d\n', ch1, ch2);
        for k = [inRFLoc exRFLoc]
%             [cohCueOnsetAll{k},phi,S12,~,~,tCue,fCue,~,~,~] = cohgramc(...
%                     squeeze(lfpAroundCueOnsetByLocAll{k}(i,:,:)), ...
%                     squeeze(lfpAroundCueOnsetByLocAll{k}(j,:,:)), ...
%                     cohgramSlidingWindow, params);
%             [cohArrayOnsetHoldAll{k},phi,S12,~,~,tArr,fArr,~,~,~] = cohgramc(...
%                     squeeze(lfpAroundArrayOnsetHoldByLocAll{k}(i,:,:)), ...
%                     squeeze(lfpAroundArrayOnsetHoldByLocAll{k}(j,:,:)), ...
%                     cohgramSlidingWindow, params);
%             [cohTargetDimAll{k},phi,S12,~,~,tDim,fDim,~,~,~] = cohgramc(...
%                     squeeze(lfpAroundTargetDimByLocAll{k}(i,:,:)), ...
%                     squeeze(lfpAroundTargetDimByLocAll{k}(j,:,:)), ...
%                     cohgramSlidingWindow, params);

            [cohLineCueOnsetAll{i,j,k},phi,S12,~,~,fCue2,~,~,~] = coherencyc(...
                    squeeze(lfpAroundCueOnsetByLocAll{k}(i,preCueLfpSelection,:)),...
                    squeeze(lfpAroundCueOnsetByLocAll{k}(j,preCueLfpSelection,:)),...
                    params);
                
            [cohLineArrayOnsetHoldAll{i,j,k},phi,S12,~,~,fArr2,~,~,~] = coherencyc(...
                    squeeze(lfpAroundArrayOnsetHoldByLocAll{k}(i,preArrayLfpSelection,:)),...
                    squeeze(lfpAroundArrayOnsetHoldByLocAll{k}(j,preArrayLfpSelection,:)),...
                    params);
                
            cohLineCueOnsetAll{j,i,k} = cohLineCueOnsetAll{i,j,k};
            cohLineArrayOnsetHoldAll{j,i,k} = cohLineArrayOnsetHoldAll{i,j,k};

%             figure;
%             hold on;
%             plot_matrix(cohCueOnsetAll{k}, tCue - cueOnsetLfpWindow(1), fCue, 'n');
%             plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%             ylim(params.fpass);
%             caxis([0 1]);
%             title(sprintf('Coherence between Channel %d and %d at P%d - Cue Onset', ch1, ch2, k));
%             plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-cueOnset-P%d.png', processedDataDir, ch1, ch2, blockName, k);
%             export_fig(plotFileName, '-nocrop');
% 
%             figure;
%             hold on;
%             plot_matrix(cohArrayOnsetHoldAll{k}, tArr - arrayOnsetLfpWindow(1), fArr, 'n');
%             plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%             ylim(params.fpass);
%             caxis([0 1]);
%             title(sprintf('Coherence between Channel %d and %d at P%d - Array Onset Hold', ch1, ch2, k));
%             plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-arrayOnset-P%d.png', processedDataDir, ch1, ch2, blockName, k);
%             export_fig(plotFileName, '-nocrop');
% 
%             figure;
%             hold on;
%             plot_matrix(cohTargetDimAll{k}, tDim - targetDimLfpWindow(1), f, 'n');
%             plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%             ylim(params.fpass);
%             caxis([0 1]);
%             title(sprintf('Coherence between Channel %d and %d at P%d - Target Dim', ch1, ch2, k));
%             plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-targetDim-P%d.png', processedDataDir, ch1, ch2, blockName, k);
%             export_fig(plotFileName, '-nocrop');
        end
        
%         figure;
%         hold on;
%         plot_matrix(cohCueOnsetAll{inRFLoc} - cohCueOnsetAll{exRFLoc}, tCue - cueOnsetLfpWindow(1), fCue, 'n');
%         plot([0 0], [0 100], '-', 'Color', 0.3*ones(3, 1));
%         xlim(cueOnsetLfpXLim);
%         ylim(params.fpass);
%         caxis(cBounds);
%         title(sprintf('Difference Coherence between Channel %d - %d - Cue Onset', ch1, ch2));
%         plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-cueOnset-diff.png', processedDataDir, ch1, ch2, blockName);
%         export_fig(plotFileName, '-nocrop');

%         figure;
%         hold on;
%         plot_matrix(cohArrayOnsetHoldAll{inRFLoc} - cohArrayOnsetHoldAll{exRFLoc}, tArr - arrayOnsetLfpWindow(1), fArr, 'n');
%         plot([0 0], [0 100], '-', 'Color', 0.3*ones(3, 1));
%         xlim(arrayOnsetLfpXLim);
%         ylim(params.fpass);
%         caxis(cBounds);
%         title(sprintf('Difference Coherence between Channel %d and %d - Array Onset Hold', ch1, ch2));
%         plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-arrayOnset-diff.png', processedDataDir, ch1, ch2, blockName);
%         export_fig(plotFileName, '-nocrop');

%         figure;
%         hold on;
%         plot_matrix(cohTargetDimAll{inRFLoc} - cohTargetDimAll{exRFLoc}, tDim - targetDimLfpWindow(1), fDim, 'n');
%         plot([0 0], [0 100], '-', 'Color', 0.3*ones(3, 1));
%         xlim(targetDimLfpXLim);
%         ylim(params.fpass);
%         caxis(cBounds);
%         title(sprintf('Difference Coherence between Channel %d and %d - Target Dim', ch1, ch2));
%         plotFileName = sprintf('%s/coh-FP%03d-FP%03d-%s-targetDim-diff.png', processedDataDir, ch1, ch2, blockName);
%         export_fig(plotFileName, '-nocrop');

%         close all;
    end
end

if params.fpass(1) == 5
    if channelsToProcess(1) == 1
        saveFileName = sprintf('%s/all-pairs-coh-%s-preCueOnset-Vprobe1.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

        saveFileName = sprintf('%s/all-pairs-coh-%s-preArrayOnset-Vprobe1.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');
    else
        saveFileName = sprintf('%s/all-pairs-coh-%s-preCueOnset-Vprobe2.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

        saveFileName = sprintf('%s/all-pairs-coh-%s-preArrayOnset-Vprobe2.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');
    end
elseif params.fpass(1) == 50
    if channelsToProcess(1) == 1
        saveFileName = sprintf('%s/all-pairs-coh-%s-preCueOnset-highGamma-Vprobe1.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

        saveFileName = sprintf('%s/all-pairs-coh-%s-preArrayOnset-highGamma-Vprobe1.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');
    else
        saveFileName = sprintf('%s/all-pairs-coh-%s-preCueOnset-highGamma-Vprobe2.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

        saveFileName = sprintf('%s/all-pairs-coh-%s-preArrayOnset-highGamma-Vprobe2.mat', processedDataDir, lower(ref));
        save(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');
    end
end

%%
saveFileName = sprintf('%s/all-pairs-coh-bip-preCueOnset-Vprobe1.mat', processedDataDir);
load(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

saveFileName = sprintf('%s/all-pairs-coh-bip-preArrayOnset-Vprobe1.mat', processedDataDir);
load(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');

%%
saveFileName = sprintf('%s/all-pairs-coh-bip-preCueOnset-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

saveFileName = sprintf('%s/all-pairs-coh-bip-preArrayOnset-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');

%%
saveFileName = sprintf('%s/all-pairs-coh-car-preCueOnset-Vprobe1.mat', processedDataDir);
load(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

saveFileName = sprintf('%s/all-pairs-coh-car-preArrayOnset-Vprobe1.mat', processedDataDir);
load(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');

%%
saveFileName = sprintf('%s/all-pairs-coh-car-preCueOnset-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

saveFileName = sprintf('%s/all-pairs-coh-car-preArrayOnset-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');

%%
saveFileName = sprintf('%s/all-pairs-coh-car-preCueOnset-highGamma-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineCueOnsetAll', 'fCue2', 'channelsToProcess', 'ref');

saveFileName = sprintf('%s/all-pairs-coh-car-preArrayOnset-highGamma-Vprobe2.mat', processedDataDir);
load(saveFileName, 'cohLineArrayOnsetHoldAll', 'fArr2', 'channelsToProcess', 'ref');

%%
% freqRanges = {[5 8], [8 15], [15 30], [15 20], [20 25], [25 30], [30 50], [40 50]};
freqRanges = {[50 100]};
for i = 1:numel(freqRanges)
    freqRange = freqRanges{i};
    fCueInds = fCue2 >= freqRange(1) & fCue2 < freqRange(2);
    fArrInds = fArr2 >= freqRange(1) & fArr2 < freqRange(2);

    % subset
    cohLineCueOnsetInRF = squeeze(cohLineCueOnsetAll(:,:,inRFLoc));
    cohLineCueOnsetExRF = squeeze(cohLineCueOnsetAll(:,:,exRFLoc));
    cohLineArrayOnsetHoldInRF = squeeze(cohLineArrayOnsetHoldAll(:,:,inRFLoc));
    cohLineArrayOnsetHoldExRF = squeeze(cohLineArrayOnsetHoldAll(:,:,exRFLoc));
    for j = 1:size(cohLineArrayOnsetHoldInRF, 1)
        cohLineCueOnsetInRF{j,j} = nan(numel(fCue2), 1);
        cohLineCueOnsetExRF{j,j} = nan(numel(fCue2), 1);
        cohLineArrayOnsetHoldInRF{j,j} = nan(numel(fArr2), 1);
        cohLineArrayOnsetHoldExRF{j,j} = nan(numel(fArr2), 1);
    end
    cueInRFCohMatrixF = cellfun(@(x) mean(x(fCueInds)), cohLineCueOnsetInRF);
    cueExRFCohMatrixF = cellfun(@(x) mean(x(fCueInds)), cohLineCueOnsetExRF);
    arrayInRFCohMatrixF = cellfun(@(x) mean(x(fArrInds)), cohLineArrayOnsetHoldInRF);
    arrayExRFCohMatrixF = cellfun(@(x) mean(x(fArrInds)), cohLineArrayOnsetHoldExRF);
    
    baselineF = (cueInRFCohMatrixF + cueInRFCohMatrixF)/2;

    if strcmp(ref, 'BIP')
        cBounds = [0 0.15];
    else
        cBounds = [min([min(arrayInRFCohMatrixF(:)) min(arrayExRFCohMatrixF(:))]) max([max(arrayInRFCohMatrixF(:)) max(arrayExRFCohMatrixF(:))])];
    end
    
    diffMaskThresh = 0.02;
    
    figure_tr_inch(14, 9);
    subaxis(2, 3, 1);
    imagesc(baselineF');
    caxis(cBounds);
    colorbar;
    axis square;
    title('Baseline');
    
    subaxis(2, 3, 2);
    imagesc(arrayInRFCohMatrixF');
    caxis(cBounds);
    colorbar;
    axis square;
    title('Delay InRF');

    subaxis(2, 3, 3);
    imagesc(arrayExRFCohMatrixF');
    caxis(cBounds);
    colorbar;
    axis square;
    title('Delay ExRF');
    
    h4 = subaxis(2, 3, 4);
    diffCohMatrixF = arrayInRFCohMatrixF' - baselineF';
    diffCohMatrixF(abs(diffCohMatrixF) < diffMaskThresh) = 0;
    imagesc(diffCohMatrixF);
    caxis([-0.1 0.1]);
    colorbar;
    colormap(h4, getCoolWarmMap());
    axis square;
    title('Delay InRF - Baseline');

    h5 = subaxis(2, 3, 5);
    diffCohMatrixF = arrayInRFCohMatrixF' - arrayExRFCohMatrixF';
    diffCohMatrixF(abs(diffCohMatrixF) < diffMaskThresh) = 0;
    imagesc(diffCohMatrixF);
    caxis([-0.1 0.1]);
    colorbar;
    colormap(h5, getCoolWarmMap());
    axis square;
    title('Delay InRF - Delay ExRF');

    suptitle(sprintf('Across-Contact Coherence (%0.1f-%0.1f Hz)', min(fArr2(fArrInds)), max(fArr2(fArrInds))));

    plotFileName = sprintf('%s/coh-mat-all-FP%03d-FP%03d-%s-%s-%0.1f-%0.1f.png', ...
            processedDataDir, channelsToProcess([1 end]), lower(ref), blockName, min(fArr2(fArrInds)), max(fArr2(fArrInds)));
    export_fig(plotFileName, '-nocrop');
end

%%
    freqRange = [15 30];
    fArrInds = fArr2 >= freqRange(1) & fArr2 < freqRange(2);

    % subset
    cohLineArrayOnsetHoldInRF = squeeze(cohLineArrayOnsetHoldAll(:,:,inRFLoc));
    cohLineArrayOnsetHoldExRF = squeeze(cohLineArrayOnsetHoldAll(:,:,exRFLoc));
    for i = 1:size(cohLineArrayOnsetHoldInRF, 1)
        cohLineArrayOnsetHoldInRF{i,i} = nan(numel(fArr2), 1);
        cohLineArrayOnsetHoldExRF{i,i} = nan(numel(fArr2), 1);
    end
    arrayInRFCohMatrixF = cellfun(@(x) mean(x(fArrInds)), cohLineArrayOnsetHoldInRF)';
    arrayExRFCohMatrixF = cellfun(@(x) mean(x(fArrInds)), cohLineArrayOnsetHoldExRF)';
    diffCohMatrixF = arrayInRFCohMatrixF - arrayExRFCohMatrixF;
    longDistanceLogical = false(size(diffCohMatrixF));
    for i = 1:(size(longDistanceLogical, 1)-15)
        longDistanceLogical(i,(i+15):size(longDistanceLogical, 1)) = true;
    end
    meanDiff = mean(diffCohMatrixF(longDistanceLogical))
    seDiff = std(diffCohMatrixF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:)))
    means = [mean(arrayInRFCohMatrixF(longDistanceLogical)) ...
            mean(arrayExRFCohMatrixF(longDistanceLogical))]
    ses = [std(arrayInRFCohMatrixF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:))) ...
            std(arrayExRFCohMatrixF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:)))]
    
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.04, 'MB', 0.15);
    set(gcf, 'Color', 'white');
    hold on;
    bh(1) = bar(1, means(1));
    bh(2) = bar(2, means(2));
    errorbar(means, ses, '.k');
        box off;
%     bh = barwitherr(ses, means);
    bh(1).FaceColor = 'r';
    bh(2).FaceColor = 'b';
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', {'Attend-RF', 'Attend-Away'});
    ylabel({'Mean Beta' 'Coherence'});
    set(gca, 'FontSize', 12);
    
    plotFileName = sprintf('%s/coh-mat-all-FP%03d-FP%03d-%s-mean3mmSep-%s-arrayOnset-%0.1f-%0.1f.png', ...
        processedDataDir, channelsToProcess([1 end]), lower(ref), blockName, min(fArr2(fArrInds)), max(fArr2(fArrInds)));
    export_fig(plotFileName, '-nocrop');
    
    %%
    cohLineArrayOnsetHoldInRF = squeeze(cohLineArrayOnsetHoldAll(:,:,inRFLoc));
    cohLineArrayOnsetHoldExRF = squeeze(cohLineArrayOnsetHoldAll(:,:,exRFLoc));
    longDistanceLogical = false(size(cohLineArrayOnsetHoldInRF));
    minContactsDist = 15;
    for i = 1:(size(longDistanceLogical, 1)-minContactsDist)
        longDistanceLogical(i,(i+minContactsDist):size(longDistanceLogical, 1)) = true;
    end
    cohInRFLongDistance = cell2mat(cohLineArrayOnsetHoldInRF(longDistanceLogical)')';
    meanCohInRFLongDistance = mean(cohInRFLongDistance);
    seCohInRFLongDistance = std(cohInRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    cohExRFLongDistance = cell2mat(cohLineArrayOnsetHoldExRF(longDistanceLogical)')';
    meanCohExRFLongDistance = mean(cohExRFLongDistance);
    seCohExRFLongDistance = std(cohExRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.19, 'MB', 0.23);
    set(gcf, 'Color', 'white');
    box off;
    hold on;
    jbfill(fArr2, meanCohInRFLongDistance + seCohInRFLongDistance, meanCohInRFLongDistance - seCohInRFLongDistance, [1 0 0], [1 0 0], 1, 0.3);
    hold on;
    jbfill(fArr2, meanCohExRFLongDistance + seCohExRFLongDistance, meanCohExRFLongDistance - seCohExRFLongDistance, [0 0 1], [0 0 1], 1, 0.3);
    hold on;
    plot(fArr2, meanCohInRFLongDistance, 'r');
    plot(fArr2, meanCohExRFLongDistance, 'b');
    if strcmp(ref, 'BIP')
        ylim([0 0.1]);
    else
        ylim([0 0.1]);%ylim([0.1 0.3]);
    end
    
    xlim([50 100]);%5 50]);
    ylabel('Coherence');
    xlabel('Frequency (Hz)');
    
    plotFileName = sprintf('%s/coh-mat-all-FP%03d-FP%03d-%s-mean3mmSep-%s-arrayOnset-%0.1f-%0.1f-line.png', ...
            processedDataDir, channelsToProcess([1 end]), lower(ref), blockName, min(fArr2), max(fArr2));
    export_fig(plotFileName, '-nocrop');
    
    %%
    cohLineArrayOnsetHoldInRF = squeeze(cohLineArrayOnsetHoldAll(:,:,inRFLoc));
    cohLineArrayOnsetHoldExRF = squeeze(cohLineArrayOnsetHoldAll(:,:,exRFLoc));
    longDistanceLogical = false(size(cohLineArrayOnsetHoldInRF));
%     for i = 
    if strcmp(ref, 'BIP')
        if channelsToProcess(1) == 1
            longDistanceLogical(1:20,28:31) = true; % vprobe1
        else
            longDistanceLogical(1:21,29:31) = true; % vprobe2
        end
    else
        if channelsToProcess(1) == 1
            longDistanceLogical(1:20,28:32) = true; % vprobe1
        else
            longDistanceLogical(1:21,29:32) = true; % vprobe2
        end
    end
%      end
    cohInRFLongDistance = cell2mat(cohLineArrayOnsetHoldInRF(longDistanceLogical)')';
    meanCohInRFLongDistance = mean(cohInRFLongDistance);
    seCohInRFLongDistance = std(cohInRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    cohExRFLongDistance = cell2mat(cohLineArrayOnsetHoldExRF(longDistanceLogical)')';
    meanCohExRFLongDistance = mean(cohExRFLongDistance);
    seCohExRFLongDistance = std(cohExRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.19, 'MB', 0.23);
    set(gcf, 'Color', 'white');
    box off;
    hold on;
    jbfill(fArr2, meanCohInRFLongDistance + seCohInRFLongDistance, meanCohInRFLongDistance - seCohInRFLongDistance, [1 0 0], [1 0 0], 1, 0.3);
    hold on;
    jbfill(fArr2, meanCohExRFLongDistance + seCohExRFLongDistance, meanCohExRFLongDistance - seCohExRFLongDistance, [0 0 1], [0 0 1], 1, 0.3);
    hold on;
    plot(fArr2, meanCohInRFLongDistance, 'r');
    plot(fArr2, meanCohExRFLongDistance, 'b');
    if strcmp(ref, 'BIP')
        ylim([0 0.1]);
    else
        ylim([0 0.1]);%ylim([0.1 0.3]);
    end
    xlim([50 100]);%xlim([5 50]);
    ylabel('Coherence');
    xlabel('Frequency (Hz)');
    
    plotFileName = sprintf('%s/coh-mat-all-FP%03d-FP%03d-%s-PLdVSPI-%s-arrayOnset-%0.1f-%0.1f-line.png', ...
            processedDataDir, channelsToProcess([1 end]), lower(ref), blockName, min(fArr2), max(fArr2));
    export_fig(plotFileName, '-nocrop');
    
    tempSaveFileName = sprintf('%s/for-grant-temp-coh3.mat', processedDataDir);
    save(tempSaveFileName, 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
    
    %%
    cohLineCueOnsetInRF = squeeze(cohLineCueOnsetAll(:,:,inRFLoc));
    cohLineCueOnsetExRF = squeeze(cohLineCueOnsetAll(:,:,exRFLoc));
    cohLineArrayOnsetHoldInRF = squeeze(cohLineArrayOnsetHoldAll(:,:,inRFLoc));
    cohLineArrayOnsetHoldExRF = squeeze(cohLineArrayOnsetHoldAll(:,:,exRFLoc));
    longDistanceLogical = false(size(cohLineArrayOnsetHoldInRF));
%     for i = 
    if strcmp(ref, 'BIP')
        if channelsToProcess(1) == 1
            longDistanceLogical(1:20,28:31) = true; % vprobe1
        else
            longDistanceLogical(1:21,29:31) = true; % vprobe2
        end
    else
        if channelsToProcess(1) == 1
            longDistanceLogical(1:20,28:32) = true; % vprobe1
        else
            longDistanceLogical(1:21,29:32) = true; % vprobe2
        end
    end
%     end
    cueCohInRFLongDistance = cell2mat(cohLineCueOnsetInRF(longDistanceLogical)')';
    meanCueCohInRFLongDistance = mean(cueCohInRFLongDistance);
    seCueCohInRFLongDistance = std(cueCohInRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    cohInRFLongDistance = cell2mat(cohLineArrayOnsetHoldInRF(longDistanceLogical)')';
    meanCohInRFLongDistance = mean(cohInRFLongDistance);
    seCohInRFLongDistance = std(cohInRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    cohExRFLongDistance = cell2mat(cohLineArrayOnsetHoldExRF(longDistanceLogical)')';
    meanCohExRFLongDistance = mean(cohExRFLongDistance);
    seCohExRFLongDistance = std(cohExRFLongDistance)./sqrt(sum(longDistanceLogical(:)));
    
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.19, 'MB', 0.23);
    set(gcf, 'Color', 'white');
    box off;
    hold on;
    jbfill(fArr2, meanCohInRFLongDistance + seCohInRFLongDistance, meanCohInRFLongDistance - seCohInRFLongDistance, [1 0 0], [1 0 0], 1, 0.3);
    hold on;
    jbfill(fArr2, meanCueCohInRFLongDistance + seCueCohInRFLongDistance, meanCueCohInRFLongDistance - seCueCohInRFLongDistance, [0 1 0], [0 1 0], 1, 0.3);
    hold on;
    plot(fArr2, meanCohInRFLongDistance, 'r');
    plot(fArr2, meanCueCohInRFLongDistance, 'g');
    if strcmp(ref, 'BIP')
        ylim([0 0.1]);
    else
        ylim([0.1 0.3]);
    end
    xlim([5 50]);
    ylabel('Coherence');
    xlabel('Frequency (Hz)');
    
    plotFileName = sprintf('%s/coh-mat-all-FP%03d-FP%03d-%s-PLdVSPI-%s-arrayOnset-%0.1f-%0.1f-line-baseline.png', ...
            processedDataDir, channelsToProcess([1 end]), lower(ref), blockName, min(fArr2), max(fArr2));
    export_fig(plotFileName, '-nocrop');
    
    %%
    meanDiff = mean(diffCohMatrixF(longDistanceLogical))
    seDiff = std(diffCohMatrixF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:)))
    means = [mean(cohLineArrayOnsetHoldInRF(longDistanceLogical)) ...
            mean(cohLineArrayOnsetHoldExRF(longDistanceLogical))]
    ses = [std(cohLineArrayOnsetHoldInRF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:))) ...
            std(cohLineArrayOnsetHoldExRF(longDistanceLogical))/sqrt(sum(longDistanceLogical(:)))]
    
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.04, 'MB', 0.15);
    set(gcf, 'Color', 'white');
    hold on;
    bh(1) = bar(1, means(1));
    bh(2) = bar(2, means(2));
    errorbar(means, ses, '.k');
        box off;
%     bh = barwitherr(ses, means);
    bh(1).FaceColor = 'r';
    bh(2).FaceColor = 'b';
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', {'Attend-RF', 'Attend-Away'});
    ylabel({'Mean Beta' 'Coherence'});
    set(gca, 'FontSize', 12);
    
    
%%
a1 = load('M20170127\for-grant-temp1.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
a2 = load('M20170327\for-grant-temp2.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
a3 = load('M20170327\for-grant-temp3.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
% a4 = load('M20170529\for-grant-temp4.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
% a5 = load('M20170529\for-grant-temp5.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
a1 = load('M20170127\for-grant-temp-coh1.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
a2 = load('M20170327\for-grant-temp-coh2.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');
a3 = load('M20170327\for-grant-temp-coh3.mat', 'meanCohInRFLongDistance', 'meanCohExRFLongDistance');

% cohInRF = [a1.meanCohInRFLongDistance; a2.meanCohInRFLongDistance; a3.meanCohInRFLongDistance; a4.meanCohInRFLongDistance; a5.meanCohInRFLongDistance];
% cohExRF = [a1.meanCohExRFLongDistance; a2.meanCohExRFLongDistance; a3.meanCohExRFLongDistance; a4.meanCohExRFLongDistance; a5.meanCohExRFLongDistance];

cohInRF = [a1.meanCohInRFLongDistance; a2.meanCohInRFLongDistance; a3.meanCohInRFLongDistance];
cohExRF = [a1.meanCohExRFLongDistance; a2.meanCohExRFLongDistance; a3.meanCohExRFLongDistance];

meanCohInRFLongDistance = mean(cohInRF);
seCohInRFLongDistance = std(cohInRF)./sqrt(size(cohInRF, 1));
meanCohExRFLongDistance = mean(cohExRF);
seCohExRFLongDistance = std(cohExRF)./sqrt(size(cohInRF, 1));

figure_tr_inch(3, 2);
subaxis(1, 1, 1, 'MT', 0.05, 'ML', 0.22, 'MB', 0.26, 'MR', 0.04);
set(gcf, 'Color', 'white');
box off;
hold on;
jbfill(fArr2, meanCohInRFLongDistance + seCohInRFLongDistance, meanCohInRFLongDistance - seCohInRFLongDistance, [1 0 0], [1 0 0], 1, 0.3);
hold on;
jbfill(fArr2, meanCohExRFLongDistance + seCohExRFLongDistance, meanCohExRFLongDistance - seCohExRFLongDistance, [0 0 1], [0 0 1], 1, 0.3);
hold on;
plot(fArr2, meanCohInRFLongDistance, 'r', 'LineWidth', 2);
plot(fArr2, meanCohExRFLongDistance, 'b', 'LineWidth', 2);
xlim([5 50]);
if strcmp(ref, 'BIP')
    ylim([0 0.1]);
else
    ylim([0.12 0.27]);
end
if params.fpass(1) == 50
    xlim([50 100]);
    ylim([0.1 0.17]);
end
ylabel('Coherence');
xlabel('Frequency (Hz)');
textParams = {'Units', 'normalized'};
text(0.64, 0.95, 'Attend-RF', 'Color', [1 0 0], textParams{:});
text(0.64, 0.85, 'Attend-Away', 'Color', [0 0 1], textParams{:});
set(gca, 'FontSize', 12);

plotFileName = sprintf('%dprobes-coh-mat-all-FP%03d-FP%03d-%s-PLdVSPI-%s-arrayOnset-%0.1f-%0.1f-line.png', ...
        size(cohInRF, 1), channelsToProcess([1 end]), lower(ref), blockName, min(fArr2), max(fArr2));
export_fig(plotFileName, '-nocrop');


%% diff
for i = 1%1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
    for j = 1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        cohDiff = (cohLineArrayOnsetHoldAll{i,j,inRFLoc} - cohLineArrayOnsetHoldAll{i,j,exRFLoc});
        plot(fArr2, cohDiff, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('InRF - ExRF Coherence between Channel %d and All - Array Onset Hold', ch1));
    plotFileName = sprintf('%s/coh-FP%03d-all-car-%s-arrayOnset-diff.png', processedDataDir, ch1, blockName);
    export_fig(plotFileName, '-nocrop');
end

%% for grant
figure_tr_inch(3, 2.5);
subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.02, 'MB', 0.2);
box off;
set(gcf, 'Color', 'white');
hold on;
ch1 = 31;
legendEntries = cell(1, 1);
plotCount = 0;
channelsToProcess = 28:-3:1;%1:3:28;
cols = parula(numel(channelsToProcess)+2); % remove extremes
cols([1 end],:) = [];
cols2 = cols;% cols(end:-1:1,:);
ch2 = 1;
cohDiff = (cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc} - cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc});
plot(fArr2, cohDiff, 'LineWidth', 2);
plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));

xlim(fArr2([1 end]));
xlabel('Frequency (Hz)');
ylabel({'Difference in Coherence', 'with Attention'});
ylim([-0.11 0.1]);
zlabel('Distance between Electrodes');
plotFileName = sprintf('%s/for-grant-coh-FP%03d-all-car-%s-arrayOnset-diff.png', processedDataDir, ch1, blockName);
export_fig(plotFileName, '-nocrop');

%%
figure_tr_inch(3.5, 2.5);
subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.26, 'MB', 0.2);
box off;
set(gcf, 'Color', 'white');
hold on;
ch1 = 31;
legendEntries = cell(1, 1);
plotCount = 0;
channelsToProcess = 28:-3:1;%1:3:28;
cols = parula(numel(channelsToProcess)+2); % remove extremes
cols([1 end],:) = [];
cols2 = cols;% cols(end:-1:1,:);
plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
for j = 1:numel(channelsToProcess)
    ch2 = channelsToProcess(j);
    if ch1 == ch2
        continue;
    end
    plotCount = plotCount + 1;
    cohDiff = (cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc} - cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc});
    plot(fArr2, cohDiff, 'Color', cols2(j,:), 'LineWidth', 2);
end

xlim(fArr2([1 end]));
xlabel('Frequency (Hz)');
ylabel({'Difference in Coherence', 'with Attention'});
ylim([-0.2 0.2]);
zlabel('Distance between Electrodes');
colorH = colorbar;
colormap(cols);
set(colorH, 'Ticks', 0.05:0.3:0.95);
set(colorH, 'TickLabels', {'0.45mm', '1.80mm', '3.15mm', '4.5mm'});
set(colorH, 'FontSize', 10);
colorH.Label.String = 'Distance between Electrodes';
colorH.Label.Rotation = 270;
colorH.Label.Position = [7 0.5 0];
colorH.Label.FontSize = 10;
colorPos = get(colorH, 'Position');
colorPos(3) = colorPos(3)*3/4;
colorPos(1) = 0.76;
set(colorH, 'Position', colorPos);
plotFileName = sprintf('%s/for-grant-coh-FP%03d-all-car-%s-arrayOnset-diff.png', processedDataDir, ch1, blockName);
export_fig(plotFileName, '-nocrop');


%% for grant 2
% fArr2(20) == 24.4 Hz
% vidObj = VideoWriter(sprintf('%s/for-grant-test-movie', processedDataDir), 'Uncompressed AVI');
% vidObj.FrameRate = 1;
% open(vidObj);
for fInd = 1:numel(fArr2)
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.05, 'MB', 0.24, 'ML', 0.2);
    box off;
    set(gcf, 'Color', 'white');
    hold on;
    ch1 = 32;
    legendEntries = cell(1, 1);
    plotCount = 0;
    channelsToProcess = 31:-1:1;%1:3:28;
    cols = parula(numel(channelsToProcess)+2); % remove extremes
    cols([1 end],:) = [];
    cols2 = cols;% cols(end:-1:1,:);
    distance = 0.15:0.15:4.65; %0.15:0.45:4.2;
    cohBaselineByDistance = nan(size(distance));
    cohInRFByDistance = nan(size(distance));
    cohExRFByDistance = nan(size(distance));
    for j = 1:numel(channelsToProcess)
        ch2 = channelsToProcess(j);
        if ch1 == ch2
            continue;
        end
        plotCount = plotCount + 1;
%         cohBaselineByDistance(plotCount) = (cohLineCueOnsetAll{ch1,ch2,inRFLoc}(fInd) + cohLineCueOnsetAll{ch1,ch2,exRFLoc}(fInd)) / 2;
        cohInRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc}(fInd);
        cohExRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc}(fInd);
    end
%     plot(distance, cohBaselineByDistance, 'Color', 'g', 'LineWidth', 2);
    plot(distance, cohInRFByDistance, 'Color', 'r', 'LineWidth', 2);
    plot(distance, cohExRFByDistance, 'Color', 'b', 'LineWidth', 2);
    xlabel('Distance between Electrodes (mm)');
    ylabel(sprintf('Coherence at %0.1f Hz', fArr2(fInd)));
    xlim([0 5]);
    ylim([0 0.4]);%ylim([0.32 0.9]);
    textParams = {'Units', 'normalized'};
    text(0.04, 0.33, 'Baseline', 'Color', 'g', textParams{:});
    text(0.04, 0.22, 'Delay InRF', 'Color', 'r', textParams{:});
    text(0.04, 0.11, 'Delay ExRF', 'Color', 'b', textParams{:});
    plotFileName = sprintf('%s/for-grant-coh-FP%03d-all-car-%s-arrayOnset-f%04.1f.png', ...
            processedDataDir, ch1, blockName, fArr2(fInd));
    export_fig(plotFileName, '-nocrop');
    
%     writeVideo(vidObj, getframe(gcf));
end
% close(vidObj);

%% for grant: bipolar
figure_tr_inch(3.5, 2.5);
subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.26, 'MB', 0.2);
box off;
set(gcf, 'Color', 'white');
hold on;
ch1 = 31;
legendEntries = cell(1, 1);
plotCount = 0;
channelsToProcess = ch1-3:-3:1;%1:3:28;
cols = parula(numel(channelsToProcess)+2); % remove extremes
cols([1 end],:) = [];
cols2 = cols;% cols(end:-1:1,:);
plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
for j = 1:numel(channelsToProcess)
    ch2 = channelsToProcess(j);
    if ch1 == ch2
        continue;
    end
    plotCount = plotCount + 1;
    cohDiff = (cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc} - cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc});
    plot(fArr2, cohDiff, 'Color', cols2(j,:), 'LineWidth', 2);
end

xlim(fArr2([1 end]));
xlabel('Frequency (Hz)');
ylabel({'Difference in Coherence', 'with Attention'});
ylim([-0.11 0.1]);
zlabel('Distance between Electrodes');
colorH = colorbar;
colormap(cols);
set(colorH, 'Ticks', 0.05:0.3:0.95);
set(colorH, 'TickLabels', {'0.15mm', '1.50mm', '2.85mm', '4.2mm'});
set(colorH, 'FontSize', 10);
colorH.Label.String = 'Distance between Electrodes';
colorH.Label.Rotation = 270;
colorH.Label.Position = [7 0.5 0];
colorH.Label.FontSize = 10;
colorPos = get(colorH, 'Position');
colorPos(3) = colorPos(3)*3/4;
colorPos(1) = 0.76;
set(colorH, 'Position', colorPos);
plotFileName = sprintf('%s/for-grant-coh-FP%03d-bip-all-%s-arrayOnset-diff.png', processedDataDir, ch1, blockName);
export_fig(plotFileName, '-nocrop');

%% for grant 2 bipolar
for fInd = find(fArr2 > 24, 1, 'first')% 1:numel(fArr2)
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.05, 'MB', 0.24, 'ML', 0.2);
    box off;
    set(gcf, 'Color', 'white');
    hold on;
    ch1 = 31;
    legendEntries = cell(1, 1);
    plotCount = 0;
    channelsToProcess = ch1:-1:1;%1:3:28;
    cols = parula(numel(channelsToProcess)+2); % remove extremes
    cols([1 end],:) = [];
    cols2 = cols;% cols(end:-1:1,:);
    distance = 0.15:0.15:4.5; %0.15:0.45:4.2;
    cohBaselineByDistance = nan(size(distance));
    cohInRFByDistance = nan(size(distance));
    cohExRFByDistance = nan(size(distance));
    for j = 1:numel(channelsToProcess)
        ch2 = channelsToProcess(j);
        if ch1 == ch2
            continue;
        end
        plotCount = plotCount + 1;
        cohBaselineByDistance(plotCount) = (cohLineCueOnsetAll{ch1,ch2,inRFLoc}(fInd) + cohLineCueOnsetAll{ch1,ch2,exRFLoc}(fInd)) / 2;
        cohInRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc}(fInd);
        cohExRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc}(fInd);
    end
    plot(distance, cohBaselineByDistance, 'Color', 'g', 'LineWidth', 2);
    plot(distance, cohInRFByDistance, 'Color', 'r', 'LineWidth', 2);
    plot(distance, cohExRFByDistance, 'Color', 'b', 'LineWidth', 2);
    xlabel('Distance between Electrodes (mm)');
    ylabel(sprintf('Coherence at %0.1f Hz', fArr2(fInd)));
    xlim([0 5]);
    ylim([0 0.15]);%ylim([0.32 0.9]);
    textParams = {'Units', 'normalized'};
    text(0.04, 0.33, 'Baseline', 'Color', 'g', textParams{:});
    text(0.04, 0.22, 'Delay InRF', 'Color', 'r', textParams{:});
    text(0.04, 0.11, 'Delay ExRF', 'Color', 'b', textParams{:});
    plotFileName = sprintf('%s/for-grant-coh-FP%03d-bip-all-%s-arrayOnset-f%04.1f.png', ...
            processedDataDir, ch1, blockName, fArr2(fInd));
    export_fig(plotFileName, '-nocrop');
end

%% for grant 3
% fArr2(20) == 24.4 Hz
for fInd = [6 20]%1:numel(fArr2)
    figure_tr_inch(3, 2);
    subaxis(1, 1, 1, 'MT', 0.05, 'MR', 0.06, 'MB', 0.24, 'ML', 0.18);
    box off;
    set(gcf, 'Color', 'white');
    hold on;
    ch1 = 31;
    legendEntries = cell(1, 1);
    plotCount = 0;
    channelsToProcess = 28:-3:1;%1:3:28;
    cols = parula(numel(channelsToProcess)+2); % remove extremes
    cols([1 end],:) = [];
    cols2 = cols;% cols(end:-1:1,:);
    distance = 0.45:0.45:4.5;
    cohBaselineByDistance = nan(size(distance));
    cohInRFByDistance = nan(size(distance));
    cohExRFByDistance = nan(size(distance));
    for j = 1:numel(channelsToProcess)
        ch2 = channelsToProcess(j);
        if ch1 == ch2
            continue;
        end
        plotCount = plotCount + 1;
        cohBaselineByDistance(plotCount) = (cohLineCueOnsetAll{ch1,ch2,inRFLoc}(fInd) + cohLineCueOnsetAll{ch1,ch2,exRFLoc}(fInd)) / 2;
        cohInRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc}(fInd);
        cohExRFByDistance(plotCount) = cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc}(fInd);
    end
    plot(distance, cohBaselineByDistance, 'Color', 'g', 'LineWidth', 3);
    plot(distance, cohInRFByDistance, 'Color', 'r', 'LineWidth', 3);
    plot(distance, cohExRFByDistance, 'Color', 'b', 'LineWidth', 3);
    xlabel('Distance between Electrodes (mm)');
    ylabel(sprintf('Coherence at %0.1f Hz', fArr2(fInd)));
    xlim([0 4.5]);
    ylim([0.5 0.8]);
    textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
    text(0.04, 0.33, 'Baseline', 'Color', 'g', textParams{:});
    text(0.04, 0.22, 'Delay Attend-RF', 'Color', 'r', textParams{:});
    text(0.04, 0.11, 'Delay Attend-Away', 'Color', 'b', textParams{:});
    set(gca, 'LineWidth', 2);
    set(gca, 'TickLength', [0.02 0.02]);
    set(gca, 'FontWeight', 'bold');
    plotFileName = sprintf('%s/for-grant-coh-FP%03d-all-%s-arrayOnset-f%04.1f-v2.png', ...
            processedDataDir, ch1, blockName, fArr2(fInd));
    export_fig(plotFileName, '-nocrop');
end

%% for grant 4
ch1 = 30;%32;
for ch2 = 1%1:31
figure_tr_inch(3, 2.2);
subaxis(1, 1, 1, 'ML', 0.24, 'MT', 0.05, 'MR', 0.05, 'MB', 0.23);
box off;
set(gcf, 'Color', 'white');
hold on;
plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
baselineF = (cohLineCueOnsetAll{ch1,ch2,inRFLoc} + cohLineCueOnsetAll{ch1,ch2,exRFLoc})/2;
inRFMinusBaselinePct = (cohLineArrayOnsetHoldAll{ch1,ch2,inRFLoc} - baselineF)./baselineF * 100;
exRFMinusBaselinePct = (cohLineArrayOnsetHoldAll{ch1,ch2,exRFLoc} - baselineF)./baselineF * 100;
plot(fArr2, inRFMinusBaselinePct, 'LineWidth', 3, 'Color', 'r');
plot(fArr2, exRFMinusBaselinePct, 'LineWidth', 3, 'Color', 'b');
xlim(fArr2([1 end]));
% xlim([fArr2(1) 30]);
ylim([-25 11]);
xlabel('Frequency (Hz)');
ylabel({'% Change in Coherence', 'from Baseline'});
textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
text(0.65, 1.0, 'Attend-RF', 'Color', 'r', textParams{:});
text(0.65, 0.9, 'Attend-Away', 'Color', 'b', textParams{:});
text(0.2, 1.0, sprintf('Ch%d vs Ch%d', channelsToProcess(ch1), channelsToProcess(ch2)), 'Color', 'k', textParams{:});

set(gca, 'LineWidth', 2);
set(gca, 'TickLength', [0.02 0.02]);
set(gca, 'FontWeight', 'bold');
plotFileName = sprintf('%s/for-grant-coh-FP%03d-FP%03d-%s-arrayOnset-diff.png', ...
        processedDataDir, channelsToProcess(ch1), channelsToProcess(ch2), blockName);
export_fig(plotFileName, '-nocrop');
end







%% diff MI
for i = 1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        cohDiffMI = (cohLineArrayOnsetHoldAll{i,j,inRFLoc} - cohLineArrayOnsetHoldAll{i,j,exRFLoc}) ./ ...
                (cohLineArrayOnsetHoldAll{i,j,inRFLoc} + cohLineArrayOnsetHoldAll{i,j,exRFLoc});
        plot(fArr2, cohDiffMI, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('InRF - ExRF MI Coherence between Channel %d and All - Array Onset Hold', ch1));
    plotFileName = sprintf('%s/coh-FP%03d-all-%s-arrayOnset-diffMI.png', processedDataDir, ch1, blockName);
    export_fig(plotFileName, '-nocrop');
end

%% inrf, exrf
for i = 31%1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        plot(fArr2, cohLineArrayOnsetHoldAll{i,j,inRFLoc}, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('InRF Coherence between Channel %d and All - Array Onset Hold', ch1));
    plotFileName = sprintf('%s/coh-FP%03d-all-%s-arrayOnset-inrf.png', processedDataDir, ch1, blockName);
    export_fig(plotFileName, '-nocrop');
    
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fArr2, zeros(size(fArr2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        plot(fArr2, cohLineArrayOnsetHoldAll{i,j,exRFLoc}, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('ExRF Coherence between Channel %d and All - Array Onset Hold', ch1));
    plotFileName = sprintf('%s/coh-FP%03d-all-%s-arrayOnset-exrf.png', processedDataDir, ch1, blockName);
    export_fig(plotFileName, '-nocrop');
end

%% baseline
for i = 31%1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fCue2, zeros(size(fCue2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        hold on;
        baselineCoh = (cohLineCueOnsetAll{i,j,inRFLoc} + cohLineCueOnsetAll{i,j,exRFLoc}) / 2;
        plot(fCue2, baselineCoh, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('Coherence between Channel %d and All - Pre-Cue Baseline', ch1));
    plotFileName = sprintf('%s/coh-FP%03d-all-%s-preCueOnset.png', processedDataDir, ch1, blockName);
    export_fig(plotFileName, '-nocrop');
end

%% exrf - baseline
for i = 31%1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fCue2, zeros(size(fCue2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        hold on;
        exRFMinusBaselineCoh = cohLineArrayOnsetHoldAll{i,j,exRFLoc} - (cohLineCueOnsetAll{i,j,inRFLoc} + cohLineCueOnsetAll{i,j,exRFLoc}) / 2;
        plot(fCue2, exRFMinusBaselineCoh, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('Coherence between Channel %d and All - ExRF Delay - Baseline', ch1));
%     plotFileName = sprintf('%s/coh-FP%03d-all-%s-preCueOnset.png', processedDataDir, ch1, blockName);
%     export_fig(plotFileName, '-nocrop');
end

%% inrf - baseline
for i = 31%1:3:32
    figure;
    hold on;
    ch1 = channelsToProcess(i);
    legendEntries = cell(1, 1);
    plotCount = 0;
    cols = parula(numel(channelsToProcess));
    plot(fCue2, zeros(size(fCue2)), 'Color', 0.3*ones(3, 1));
    for j = 1:3:32%1:numel(channelsToProcess)
        if i == j
            continue;
        end
        plotCount = plotCount + 1;
        legendEntries{plotCount} = sprintf('Ch %d', channelsToProcess(j));
        hold on;
        inRFMinusBaselineCoh = cohLineArrayOnsetHoldAll{i,j,inRFLoc} - (cohLineCueOnsetAll{i,j,inRFLoc} + cohLineCueOnsetAll{i,j,exRFLoc}) / 2;
        plot(fCue2, inRFMinusBaselineCoh, 'Color', cols(j,:));
    end
    % legend(legendEntries);

    xlim(fArr2([1 end]));
    title(sprintf('Coherence between Channel %d and All - InRF Delay - Baseline', ch1));
%     plotFileName = sprintf('%s/coh-FP%03d-all-%s-preCueOnset.png', processedDataDir, ch1, blockName);
%     export_fig(plotFileName, '-nocrop');
end


%% compute power
load('params.mat');
params.fpass = [5 50];
params.tapers = [1 1];
params.pad = 1;

cohLinePreCueWindowOffset = [-0.3 0];
preCueLfpSelection = cueOnsetLfpT >= cohLinePreCueWindowOffset(1)*Fs & ...
        cueOnsetLfpT <= cohLinePreCueWindowOffset(2)*Fs;
powLineCueOnsetAll = cell(numel(channelsToProcess), nLoc);

cohLinePreArrayWindowOffset = [-0.3 0];
preArrayLfpSelection = arrayOnsetLfpT >= cohLinePreArrayWindowOffset(1)*Fs & ...
        arrayOnsetLfpT <= cohLinePreArrayWindowOffset(2)*Fs;
powLineArrayOnsetHoldAll = cell(numel(channelsToProcess), nLoc);

for i = 1:numel(channelsToProcess)
    ch1 = channelsToProcess(i);

    fprintf('Processing power of channel %d\n', ch1);
    for k = [inRFLoc exRFLoc]
        [powLineCueOnsetAll{i,k},fCue3,~] = mtspectrumc(...
                squeeze(lfpAroundCueOnsetByLocAll{k}(i,preCueLfpSelection,:)),...
                params);

        [powLineArrayOnsetHoldAll{i,k},fArr3,~] = mtspectrumc(...
                squeeze(lfpAroundArrayOnsetHoldByLocAll{k}(i,preArrayLfpSelection,:)),...
                params);
    end
end

saveFileName = sprintf('%s/all-chs-pow-preCueOnset.mat', processedDataDir);
save(saveFileName, 'powLineCueOnsetAll', 'fCue3');

saveFileName = sprintf('%s/all-chs-pow-preArrayOnset.mat', processedDataDir);
save(saveFileName, 'powLineArrayOnsetHoldAll', 'fArr3');


%% diff
figure;
hold on;
plotCount = 0;
cols = parula(numel(channelsToProcess));
plot(fArr3, zeros(size(fArr3)), 'Color', 0.3*ones(3, 1));
for i = 1:numel(channelsToProcess)
    ch1 = channelsToProcess(i);
    plotCount = plotCount + 1;
    powDiff = (10*log10(powLineArrayOnsetHoldAll{i,inRFLoc}) - 10*log10(powLineArrayOnsetHoldAll{i,exRFLoc}));
    plot(fArr3, powDiff, 'Color', cols(i,:), 'LineWidth', 2);
end
xlim(fArr3([1 end]));
xlabel('Frequency (Hz)');
ylabel('Difference in Power (dB)');
title(sprintf('InRF - ExRF Power - Array Onset Hold'));
plotFileName = sprintf('%s/pow-all-%s-arrayOnset-diff.png', processedDataDir, blockName);
export_fig(plotFileName, '-nocrop');

%% inrf, exrf
figure;
hold on;
i = find(channelsToProcess == 28);
ch1 = channelsToProcess(i);
plot(fCue3, 10*log10((powLineCueOnsetAll{i,inRFLoc} + powLineCueOnsetAll{i,exRFLoc})/2), 'Color', 'g', 'LineWidth', 2);
plot(fArr3, 10*log10(powLineArrayOnsetHoldAll{i,inRFLoc}), 'Color', 'r', 'LineWidth', 2);
plot(fArr3, 10*log10(powLineArrayOnsetHoldAll{i,exRFLoc}), 'Color', 'b', 'LineWidth', 2);
xlim(fArr3([1 end]));
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
textParams = {'Units', 'normalized', 'FontWeight', 'bold'};
text(0.04, 0.3, 'Baseline', 'Color', 'g', textParams{:});
text(0.04, 0.2, 'Delay InRF', 'Color', 'r', textParams{:});
text(0.04, 0.1, 'Delay ExRF', 'Color', 'b', textParams{:});
plotFileName = sprintf('%s/pow-FP%03d-all-%s-arrayOnset-all.png', processedDataDir, ch1, blockName);
export_fig(plotFileName, '-nocrop');


%%
% cueOnsetLfpXLim = [-0.4 0.7] * Fs;
% figure;
% hold on;
% for i = 1:size(lfpAroundCueOnset, 2)
%     plot(cueOnsetLfpT, lfpAroundCueOnset(:,i));
% end
% 
% plot(cueOnsetLfpT, mean(lfpAroundCueOnset, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(cueOnsetLfpXLim);
% xlabel('Time from Cue Onset (ms)');
% title(sprintf('Cue Evoked Potential, Raw and Average (N=%d)', size(lfpAroundCueOnset, 2)));

%%
% figure;
% hold on;
% 
% cueOnsetLfpXLim = [-0.4 0.7] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundCueOnsetByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(cueOnsetLfpT, mean(lfpAroundCueOnsetByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(cueOnsetByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(cueOnsetLfpXLim);
% xlabel('Time from Cue Onset (ms)');
% title(sprintf('Cue Evoked Potential by Location (N=%d)', size(lfpAroundCueOnset, 2)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%

% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% figure;
% hold on;
% for i = 1:size(lfpAroundArrayOnset, 2)
%     plot(arrayOnsetLfpT, lfpAroundArrayOnset(:,i));
% end
% 
% plot(arrayOnsetLfpT, mean(lfpAroundArrayOnset, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential (N=%d)', size(lfpAroundArrayOnset, 2)));

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundArrayOnsetHoldByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetHoldByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target (N=%d)', numel(arrayOnsetHold)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     assert(~any(any(isnan(lfpAroundArrayOnsetShortHoldByLoc{i}))));
%     if isempty(lfpAroundArrayOnsetShortHoldByLocFilt{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetShortHoldByLocFilt{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetShortHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target Short Dur (N=%d)', nTrialShortHoldDur));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 1] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     assert(~any(any(isnan(lfpAroundArrayOnsetLongHoldByLoc{i}))));
%     if isempty(lfpAroundArrayOnsetLongHoldByLocFilt{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetLongHoldByLocFilt{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetLongHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Hold Target Long Dur (N=%d)', nTrialLongHoldDur));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
% figure;
% hold on;
% 
% arrayOnsetLfpXLim = [-0.7 0.6] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundArrayOnsetRelByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(arrayOnsetLfpT, mean(lfpAroundArrayOnsetRelByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetRelByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(arrayOnsetLfpXLim);
% xlabel('Time from Array Onset (ms)');
% title(sprintf('Array Evoked Potential by Location -- Release Target (N=%d)', numel(arrayOnsetRel)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');
    
%%
% figure;
% hold on;
% for i = 1:size(lfpAroundTargetDim, 2)
%     plot(arrayOnsetLfpT, lfpAroundTargetDim(:,i));
% end
% 
% plot(targetDimLfpT, mean(lfpAroundTargetDim, 2), 'LineWidth', 5, 'Color', 'k');
% xlim(targetDimLfpXLim);
% xlabel('Time from Target Dimming (ms)');
% title(sprintf('Target Dimming Potential (N=%d)', size(lfpAroundTargetDim, 2)));

%%
% figure;
% hold on;
% 
% targetDimLfpXLim = [-0.7 0.4] * Fs;
% legendEntry = cell(nLoc, 1);
% fHandles = nan(nLoc, 1);
% for i = 1:nLoc
%     if isempty(lfpAroundTargetDimByLoc{i})
%         continue;
%     end
%     fHandles(i) = plot(targetDimLfpT, mean(lfpAroundTargetDimByLoc{i}, 2), 'LineWidth', 2);
%     legendEntry{i} = sprintf('P%d (N=%d)', i, numel(arrayOnsetHoldByLoc{i}));
% end
% legendEntry(isnan(fHandles)) = [];
% fHandles(isnan(fHandles)) = [];
% 
% origYLim = ylim();
% plot([0 0], [-100 100], '-', 'Color', [0.5 0.5 0.5]);
% ylim(origYLim);
% xlim(targetDimLfpXLim);
% xlabel('Time from Target Dimming (ms)');
% title(sprintf('Target Dimming Evoked Potential by Location (N=%d)', numel(targetDim)));
% legend(fHandles, legendEntry, 'Location', 'NorthWest');

%%
load('params.mat');
params.fpass = [5 50];
params.pad = 2;
specgramSlidingWindow = [0.3 0.05];

%% plot coherence between pul and v4 channels
lfpTs = eval(sprintf('FP%03d_ts', startChannel));
lfpInd = eval(sprintf('FP%03d_ind', startChannel));
lfpTsStep = eval(sprintf('FP%03d_ts_step', startChannel)); % sampling period
v4ChannelData = cell(nChannels, 1);
v4LfpAroundCueOnsetByLoc = cell(nChannels, nLoc);
v4LfpAroundArrayOnsetByLoc = cell(nChannels, nLoc);
v4LfpAroundArrayOnsetHoldByLoc = cell(nChannels, nLoc);
v4LfpAroundTargetDimByLoc = cell(nChannels, nLoc);

cueOnsetLfpWindowNew = [0.15 0.5];
arrayOnsetLfpWindowNew = [0.5 0.15];
targetDimLfpWindowNew = [0.5 0.15];

cueOnsetLfpTNew = -cueOnsetLfpWindowNew(1) * Fs : cueOnsetLfpWindowNew(2) * Fs - 1;
arrayOnsetLfpTNew = -arrayOnsetLfpWindowNew(1) * Fs : arrayOnsetLfpWindowNew(2) * Fs - 1;
targetDimLfpTNew = -targetDimLfpWindowNew(1) * Fs : targetDimLfpWindowNew(2) * Fs - 1;

supChannel = supChannelOrig + xChannelAlignmentOffset;
deepChannel = deepChannelOrig + xChannelAlignmentOffset;

for j = [supChannel deepChannel]%1:nChannels
    channelName1 = sprintf('FP%03d', startChannel + j - 1);
    channelName2 = sprintf('FP%03d', startChannel + j - 1 - 1); % one above
    v4ChannelData{j} = eval(channelName1) - eval(channelName2); % bipolar LFP <----------------
    v4ChannelData{j} = padNaNsToAdjustLfpOffset(v4ChannelData{j}, lfpTs, lfpInd, 1/lfpTsStep);
    for i = [inRFLoc exRFLoc]
        v4LfpAroundCueOnsetByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                cueOnsetByLoc{i}, Fs, cueOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundCueOnsetByLoc{j,i}))));
        v4LfpAroundArrayOnsetByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                arrayOnsetByLoc{i}, Fs, arrayOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundArrayOnsetByLoc{j,i}))));
        v4LfpAroundArrayOnsetHoldByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundArrayOnsetHoldByLoc{j,i}))));
        v4LfpAroundTargetDimByLoc{j,i} = createdatamatc(v4ChannelData{j}, ...
                targetDimByLoc{i}, Fs, targetDimLfpWindowNew);
        assert(~any(any(isnan(v4LfpAroundTargetDimByLoc{j,i}))));
    end
end
% pulvinar
lfpAroundCueOnsetByLocNew = cell(nLoc, 1);
lfpAroundArrayOnsetByLocNew = cell(nLoc, 1);
lfpAroundArrayOnsetHoldByLocNew = cell(nLoc, 1);
lfpAroundTargetDimByLocNew = cell(nLoc, 1);
for i = [inRFLoc exRFLoc]
    lfpAroundCueOnsetByLocNew{i} = createdatamatc(adjLfp, ...
            cueOnsetByLoc{i}, Fs, cueOnsetLfpWindowNew);
    lfpAroundArrayOnsetByLocNew{i} = createdatamatc(adjLfp, ...
            arrayOnsetByLoc{i}, Fs, arrayOnsetLfpWindowNew);
    lfpAroundArrayOnsetHoldByLocNew{i} = createdatamatc(adjLfp, ...
            arrayOnsetHoldByLoc{i}, Fs, arrayOnsetLfpWindowNew);
    lfpAroundTargetDimByLocNew{i} = createdatamatc(adjLfp, ...
            targetDimByLoc{i}, Fs, targetDimLfpWindowNew);
end

%%
load('params.mat');
params.fpass = [5 50];
params.pad = 2;

cohgramSlidingWindow = [0.3 0.025];
cohLinePreArrayWindowOffset = [-0.2 0];
cohLinePreDimWindowOffset = [-0.2 0];
preArrayLfpSelection = arrayOnsetLfpTNew >= cohLinePreArrayWindowOffset(1)*Fs & ...
        arrayOnsetLfpTNew <= cohLinePreArrayWindowOffset(2)*Fs;
preDimLfpSelection = targetDimLfpTNew >= cohLinePreDimWindowOffset(1)*Fs & ...
        targetDimLfpTNew <= cohLinePreDimWindowOffset(2)*Fs;

for k = [inRFLoc exRFLoc]
    for j = [supChannel deepChannel]%1:nChannels
        fprintf('Processing channel %d...\n', j);
%         [cohCueOnsetAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundCueOnsetByLocNew{k}, ...
%                 v4LfpAroundCueOnsetByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohArrayOnsetAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundArrayOnsetByLocNew{k}, ...
%                 v4LfpAroundArrayOnsetByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohArrayOnsetHoldAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundArrayOnsetHoldByLocNew{k}, ...
%                 v4LfpAroundArrayOnsetHoldByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
%         [cohTargetDimAll{l,k,j},phi,S12,~,~,t,f,confC,phistd,Cerr] = cohgramc(...
%                 lfpAroundTargetDimByLocNew{k}, ...
%                 v4LfpAroundTargetDimByLoc{j,k}, ...
%                 cohgramSlidingWindow, params);
        % spectrums only
        [cohLineArrayOnsetAll{l,k,j},phi,S12,S1,S2,fSpectrumCohLine] = coherencyc(...
                lfpAroundArrayOnsetByLocNew{k}(preArrayLfpSelection,:), ...
                v4LfpAroundArrayOnsetByLoc{j,k}(preArrayLfpSelection,:), ...
                params);
        [cohLineTargetDimAll{l,k,j},phi,S12,S1,S2,fSpectrumCohLine] = coherencyc(...
                lfpAroundTargetDimByLocNew{k}(preDimLfpSelection,:), ...
                v4LfpAroundTargetDimByLoc{j,k}(preDimLfpSelection,:), ...
                params);
        
%         figure;
%         hold on;
%         plot_matrix(cohAll{k,j}, t - arrayOnsetLfpWindowNew(1), f);
%         plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%         ylim(params.fpass);
%         title(sprintf('Coherence between Pulvinar and V4 Channel %d at P%d', j, k));
    end
end

%% inrf - exrf
% for j = [10 30]
%     figure;
%     hold on;
%     imagesc(t - cueOnsetLfpWindowNew(1), f, (cohCueOnsetAll{l,inRFLoc,j} - cohCueOnsetAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-postCue.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     imagesc(t - arrayOnsetLfpWindowNew(1), f, (cohArrayOnsetAll{l,inRFLoc,j} - cohArrayOnsetAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-preArray.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     imagesc(t - targetDimLfpWindowNew(1), f, (cohTargetDimAll{l,inRFLoc,j} - cohTargetDimAll{l,exRFLoc,j})');
%     axis xy;
%     plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
%     colorbar;
% %     caxis([-max(abs(caxis)) max(abs(caxis))]);
%     caxis([-0.15 0.15]);
%     colormap(getCoolWarmMap());
%     ylim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCoh-preDim.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     plot(fSpectrumCohLine, (cohLineArrayOnsetAll{l,inRFLoc,j} - cohLineArrayOnsetAll{l,exRFLoc,j}));
%     xlim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCohLine-preArray.png', sessionName, j);
%     export_fig(plotFileName);
%     
%     figure;
%     hold on;
%     plot(fSpectrumCohLine, (cohLineTargetDimAll{l,inRFLoc,j} - cohLineTargetDimAll{l,exRFLoc,j}));
%     xlim(params.fpass);
%     title(sprintf('%s Diff Coherence between Pulvinar and V4 Channel %d', sessionName, j));
%     
%     plotFileName = sprintf('%s-PUL-V4-ch%d-diffCohLine-preDim.png', sessionName, j);
%     export_fig(plotFileName);
% end


% end % end numSessions

%%
% % params.fpass = [30 50];
% % yAxis = [-35 -27];
% figure_tr_inch(8, 6);
% subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
% hold on;
% % plot_vector(mean(SPreCueInRFAll, 1), fCue, 'l', [], 'k');
% h1 = plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1)), ...
%         'Color', [1 0 0], ...
%         'LineWidth', 4);
% % plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% % plot(fSpectrum, 10*log10(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% h2 = plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1)), ...
%         'Color', [0 0 1], ...
%         'LineWidth', 4);
% % plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) + std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% % plot(fSpectrum, 10*log10(mean(SPreArrayExRFAll, 1) - std(SPreArrayExRFAll, 0, 1)/sqrt(numel(sessionRange) - 1)));
% h3 = plot(fSpectrum, 10*log10(mean(SPreDimInRFAll, 1)), ...
%         'Color', [1 0 1], ...
%         'LineWidth', 4);
% h4 = plot(fSpectrum, 10*log10(mean(SPreDimExRFAll, 1)), ...
%         'Color', [0 1 1], ...
%         'LineWidth', 4);
% % plot_vector(mean(SPreArrayInRFAll, 1), fSpectrum, 'l', [], 'r');
% % plot_vector(mean(SPreArrayInRFAll, 1) + std(SPreArrayInRFAll, 0, 1) / sqrt(size(SPreArrayInRFAll, 1)), fSpectrum, 'l', [], 'k');
% % plot_vector(mean(SPreArrayInRFAll, 1) - std(SPreArrayInRFAll, 0, 1) / sqrt(size(SPreArrayInRFAll, 1)), fSpectrum, 'l', [], 'k');
% % plot_vector(mean(SPreArrayExRFAll, 1), fSpectrum, 'l', [], 'b');
% % plot_vector(mean(SPreDimInRFAll, 1), fSpectrum, 'l', [], 'm');
% % plot_vector(mean(SPreDimExRFAll, 1), fSpectrum, 'l', [], 'c');
% xlim(params.fpass);
% ylim(yAxis);
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');
% legend([h1 h2 h3 h4], ...
%         {sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend In RF', num2str(h1.Color))...
%         sprintf('\\color[rgb]{%s}Pre-Stimulus, Attend Out RF', num2str(h2.Color)), ...
%         sprintf('\\color[rgb]{%s}Post-Stimulus, Attend In RF', num2str(h3.Color)), ...
%         sprintf('\\color[rgb]{%s}Post-Stimulus, Attend Out RF', num2str(h4.Color))}, ...
%         'Position', [0.5 0.78 0.4 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
% set(gca, 'FontName', 'Calibri')
% set(gca, 'FontSize', 26);
% set(gca, 'FontWeight', 'bold');
% set(gca, 'Box', 'off');
% set(gca, 'LineWidth', 3);
% set(gcf, 'Color', 'w');
% 
% plotFileName = sprintf('meanSessions-%s-powerComparison-%0.1f-%0.1fHz-n%d.png', ...
%         lfpVarName, params.fpass, numel(sessionRange));
% export_fig(plotFileName);
% 
% stop

%%

% TODO add post-cue 
l = sessionRange(1);
xChannelAlignmentOffset = sessions{l,4};
supChannel = supChannelOrig + xChannelAlignmentOffset;
deepChannel = deepChannelOrig + xChannelAlignmentOffset;

cohPreArrayDiffForAvgSup = nan(numSessions, size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohPreArrayDiffForAvgDeep = nan(numSessions, size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohPreDimDiffForAvgSup = nan(numSessions, size(cohTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohTargetDimAll{l,inRFLoc,supChannel}, 2));
cohPreDimDiffForAvgDeep = nan(numSessions, size(cohTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayDiffForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayDiffForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimDiffForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimDiffForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

cohLinePreArrayVsPreDimInRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimExRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimInRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayVsPreDimExRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

cohLinePreArrayInRFForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayExRFForAvgSup = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayInRFForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreArrayExRFForAvgDeep = nan(numSessions, size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 1), size(cohLineArrayOnsetAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimInRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimExRFForAvgSup = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimInRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));
cohLinePreDimExRFForAvgDeep = nan(numSessions, size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 1), size(cohLineTargetDimAll{l,inRFLoc,supChannel}, 2));

for l = sessionRange
    xChannelAlignmentOffset = sessions{l,4};
    supChannel = supChannelOrig + xChannelAlignmentOffset;
    deepChannel = deepChannelOrig + xChannelAlignmentOffset;

    cohPreArrayDiffForAvgSup(l,:,:) = (cohArrayOnsetAll{l,inRFLoc,supChannel} - cohArrayOnsetAll{l,exRFLoc,supChannel});
    cohPreArrayDiffForAvgDeep(l,:,:) = (cohArrayOnsetAll{l,inRFLoc,deepChannel} - cohArrayOnsetAll{l,exRFLoc,deepChannel});
    cohPreDimDiffForAvgSup(l,:,:) = (cohTargetDimAll{l,inRFLoc,supChannel} - cohTargetDimAll{l,exRFLoc,supChannel});
    cohPreDimDiffForAvgDeep(l,:,:) = (cohTargetDimAll{l,inRFLoc,deepChannel} - cohTargetDimAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayDiffForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,supChannel} - cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayDiffForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,deepChannel} - cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    cohLinePreDimDiffForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel} - cohLineTargetDimAll{l,exRFLoc,supChannel});
    cohLinePreDimDiffForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel} - cohLineTargetDimAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayVsPreDimInRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel} - cohLineArrayOnsetAll{l,inRFLoc,supChannel});
    cohLinePreArrayVsPreDimExRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,supChannel} - cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayVsPreDimInRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel} - cohLineArrayOnsetAll{l,inRFLoc,deepChannel});
    cohLinePreArrayVsPreDimExRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,deepChannel} - cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    
    cohLinePreArrayInRFForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,supChannel});
    cohLinePreArrayExRFForAvgSup(l,:,:) = (cohLineArrayOnsetAll{l,exRFLoc,supChannel});
    cohLinePreArrayInRFForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,inRFLoc,deepChannel});
    cohLinePreArrayExRFForAvgDeep(l,:,:) = (cohLineArrayOnsetAll{l,exRFLoc,deepChannel});
    cohLinePreDimInRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,supChannel});
    cohLinePreDimExRFForAvgSup(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,supChannel});
    cohLinePreDimInRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,inRFLoc,deepChannel});
    cohLinePreDimExRFForAvgDeep(l,:,:) = (cohLineTargetDimAll{l,exRFLoc,deepChannel});
    
end

meanCohPreArrayDiffSup = squeeze(nanmean(cohPreArrayDiffForAvgSup, 1));
meanCohPreArrayDiffDeep = squeeze(nanmean(cohPreArrayDiffForAvgDeep, 1));
meanCohPreDimDiffSup = squeeze(nanmean(cohPreDimDiffForAvgSup, 1));
meanCohPreDimDiffDeep = squeeze(nanmean(cohPreDimDiffForAvgDeep, 1));

meanCohLinePreArrayDiffSup = squeeze(nanmean(cohLinePreArrayDiffForAvgSup, 1));
meanCohLinePreArrayDiffDeep = squeeze(nanmean(cohLinePreArrayDiffForAvgDeep, 1));
meanCohLinePreDimDiffSup = squeeze(nanmean(cohLinePreDimDiffForAvgSup, 1));
meanCohLinePreDimDiffDeep = squeeze(nanmean(cohLinePreDimDiffForAvgDeep, 1));

meanCohLinePreArrayVsPreDimInRFSup = squeeze(nanmean(cohLinePreArrayVsPreDimInRFForAvgSup, 1));
meanCohLinePreArrayVsPreDimExRFSup = squeeze(nanmean(cohLinePreArrayVsPreDimExRFForAvgSup, 1));
meanCohLinePreArrayVsPreDimInRFDeep = squeeze(nanmean(cohLinePreArrayVsPreDimInRFForAvgDeep, 1));
meanCohLinePreArrayVsPreDimExRFDeep = squeeze(nanmean(cohLinePreArrayVsPreDimExRFForAvgDeep, 1));

meanCohLinePreArrayInRFSup = squeeze(nanmean(cohLinePreArrayInRFForAvgSup));
meanCohLinePreArrayExRFSup = squeeze(nanmean(cohLinePreArrayExRFForAvgSup));
meanCohLinePreArrayInRFDeep = squeeze(nanmean(cohLinePreArrayInRFForAvgDeep));
meanCohLinePreArrayExRFDeep = squeeze(nanmean(cohLinePreArrayExRFForAvgDeep));
meanCohLinePreDimInRFFSup = squeeze(nanmean(cohLinePreDimInRFForAvgSup));
meanCohLinePreDimExRFSup = squeeze(nanmean(cohLinePreDimExRFForAvgSup));
meanCohLinePreDimInRFDeep = squeeze(nanmean(cohLinePreDimInRFForAvgDeep));
meanCohLinePreDimExRFDeep = squeeze(nanmean(cohLinePreDimExRFForAvgDeep));

%%
cohLineYLim = [0 0.4];
cohLineXLim = [5 50];

for l = sessionRange
    sessionName = sessions{l,1};
    xChannelAlignmentOffset = sessions{l,4};
    supChannel = supChannelOrig + xChannelAlignmentOffset;
    deepChannel = deepChannelOrig + xChannelAlignmentOffset;
    
figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreArrayInRFForAvgSup(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreArrayExRFForAvgSup(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preArray-n%d.png', ...
        sessionName, lfpVarName, supChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreDimInRFForAvgSup(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreDimExRFForAvgSup(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preDim-n%d.png', ...
        sessionName, lfpVarName, supChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreArrayInRFForAvgDeep(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreArrayExRFForAvgDeep(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preArray-n%d.png', ...
        sessionName, lfpVarName, deepChannel, numel(sessionRange));
export_fig(plotFileName);


figure_tr_inch(8, 6);
subaxis(1,1,1,'MR',0.05,'ML',0.15,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, squeeze(cohLinePreDimInRFForAvgDeep(l,:,:)), ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, squeeze(cohLinePreDimExRFForAvgDeep(l,:,:)), ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
legend({'Attend In RF', ...
        'Attend Out RF'}, ...
        'Position', [0.66 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('%s-%s-V4-ch%d-cohLine-preDim-n%d.png', ...
        sessionName, lfpVarName, deepChannel, numel(sessionRange));
export_fig(plotFileName);
    
    
end % for each session

%%
cohLineYLim = [0 0.17];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreArrayInRFSup, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayExRFSup, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
% xlabel('Frequency (Hz)');
ylabel('Coherence');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArray-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
h1 = plot(fSpectrumCohLine, meanCohLinePreDimInRFFSup, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreDimExRFSup, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Attend In RF', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Attend Out RF', num2str(h2.Color))}, ...
        'Position', [0.65 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preDim-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreArrayInRFDeep, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayExRFDeep, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArray-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(fSpectrumCohLine, meanCohLinePreDimInRFDeep, ...
        'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreDimExRFDeep, ...
        'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preDim-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');

%%
cohLineYLim = [-0.04 0.12];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
h1 = plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimInRFSup, ...
        ':', 'Color', [1 0 0], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimExRFSup, ...
        ':', 'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Attend In RF (Post-Pre Stimulus)', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Attend Out RF (Post-Pre Stimulus)', num2str(h2.Color))}, ...
        'Position', [0.45 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArrayVsPreDim-n%d.png', ...
        lfpVarName, supChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimInRFDeep, ...
        ':', 'Color', [1 0 0], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayVsPreDimExRFDeep, ...
        ':', 'Color', [0 0 1], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence DIfference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-cohLine-preArrayVsPreDim-n%d.png', ...
        lfpVarName, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


%%
cohLineYLim = [-0.05 0.07];

figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
plot(fSpectrumCohLine, meanCohLinePreArrayDiffSup, ...
        'Color', [0.6 0.1 0.9], ...
        'LineWidth', 4);
plot(fSpectrumCohLine, meanCohLinePreArrayDiffDeep, ...
        'Color', [0 0.7 0.3], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
ylabel('Coherence Difference');
% legend({'Attend In RF', ...
%         'Attend Out RF'}, ...
%         'Position', [0.66 0.82 0.3 0.1], ...
%         'FontSize', 24, ...
%         'Box', 'off');
set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-chOff%d-diff-cohLine-preArray-n%d.png', ...
        lfpVarName, supChannelOrig, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


figure_tr_inch(7, 5);
subaxis(1,1,1,'MR',0.05,'ML',0.18,'MB',0.18,'MT',0.1); 
hold on;
plot(cohLineXLim, [0 0], 'Color', ones(3,1)*0.1, 'LineWidth', 1);
h1 = plot(fSpectrumCohLine, meanCohLinePreDimDiffSup, ...
        'Color', [0.6 0.1 0.9], ...
        'LineWidth', 4);
h2 = plot(fSpectrumCohLine, meanCohLinePreDimDiffDeep, ...
        'Color', [0 0.7 0.3], ...
        'LineWidth', 4);
xlim(cohLineXLim);
ylim(cohLineYLim);
xlabel('Frequency (Hz)');
% ylabel('Coherence');

set(gca, 'FontName', 'Calibri')
set(gca, 'FontSize', 26);
set(gca, 'FontWeight', 'bold');
set(gca, 'Box', 'off');
set(gca, 'LineWidth', 3);
set(gcf, 'Color', 'w');
legend([h1 h2], ...
        {sprintf('\\color[rgb]{%s}Pul & Superficial V4', num2str(h1.Color)), ...
        sprintf('\\color[rgb]{%s}Pul & Deep V4', num2str(h2.Color))}, ...
        'Position', [0.6 0.82 0.3 0.1], ...
        'FontSize', 24, ...
        'Box', 'off');

plotFileName = sprintf('meanSessions-%s-V4-chOff%d-chOff%d-diff-cohLine-preDim-n%d.png', ...
        lfpVarName, supChannelOrig, deepChannelOrig, numel(sessionRange));
export_fig(plotFileName, '-nocrop');


stop

%%
maxAbsCAxis = 0.1;
figure_tr_inch(12,6);

subaxis(1, 2, 1);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreArrayDiffSup');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        supChannel, numel(sessionRange)));

subaxis(1, 2, 2);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreArrayDiffDeep');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        deepChannel, numel(sessionRange)));

plotFileName = sprintf('meanSessions-PUL-V4-ch%d-ch%d-meanDiffCoh-preArray-n%d.png', ...
        supChannel, deepChannel, numel(sessionRange));
export_fig(plotFileName);

figure_tr_inch(12,6);

subaxis(1, 2, 1);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreDimDiffSup');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        supChannel, numel(sessionRange)));

subaxis(1, 2, 2);
hold on;
imagesc(t - arrayOnsetLfpWindowNew(1), f, meanCohPreDimDiffDeep');
axis xy;
plot([0 0], [0 100], '-', 'Color', [0.5 0.5 0.5]);
colorbar;
%     caxis([-max(abs(caxis)) max(abs(caxis))]);
caxis([-maxAbsCAxis maxAbsCAxis]);
colormap(getCoolWarmMap());
ylim(params.fpass);
title(sprintf('Mean Diff Coherence between Pulvinar and V4 Channel %d (n=%d)', ...
        deepChannel, numel(sessionRange)));

plotFileName = sprintf('meanSessions-PUL-V4-ch%d-ch%d-meanDiffCoh-preDim-n%d.png', ...
        supChannel, deepChannel, numel(sessionRange));
export_fig(plotFileName);

