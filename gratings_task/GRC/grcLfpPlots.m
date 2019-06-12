%%
clear;
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\LFP_GRATINGS_SUMMARY';
ref = 'CAR';
v = 12;

saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-%s-v%d.mat', outputDir, ref, v);
load(saveFileName);
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\GRC2';

%%
cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

if strcmp(ref, 'BIP')
    relPowYBounds = [-0.55 0.1];
    sfcCTDelayYBoundsLF = [0.035 0.065];
    sfcCTDelayYBoundsHF = [0.0125 0.03];
elseif strcmp(ref, 'CAR')
    relPowYBounds = [-1.5 0.2];
    sfcCTDelayYBoundsLF = [0.035 0.07];
    sfcCTDelayYBoundsHF = [0.015 0.03];
end


isInPulvinar = isInDPulvinar | isInVPulvinar;

%% plot power in cue-target delay allPul P3 vs P1
channelCond = isInPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelPowerLF = [cueTargetDelayPowerP3LF(channelCond,:); cueTargetDelayPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
cueTargetDelayRelPowerHF = [cueTargetDelayPowerP3HF(channelCond,:); cueTargetDelayPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(cueTargetDelayRelPowerLF, cueTargetDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', '', ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-pul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

% %% plot power in cue-target delay dPul P3 vs P1
% channelCond = isInDPulvinar;
% nChannel = sum(channelCond);
% cueTargetDelayRelPowerLF = [cueTargetDelayPowerP3LF(channelCond,:); cueTargetDelayPowerP1LF(channelCond,:)] ./ ...
%         [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
% cueTargetDelayRelPowerHF = [cueTargetDelayPowerP3HF(channelCond,:); cueTargetDelayPowerP1HF(channelCond,:)] ./ ...
%         [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
% p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];
% 
% grcPlotLfpPower2(cueTargetDelayRelPowerLF, cueTargetDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
%         'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
%         'yBounds', [relPowYBounds; relPowYBounds], ...
%         'cols', [p3Col; p1Col], ...
%         'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
%         'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
%         'titleText', '', ...
%         'doDB', 1);
% 
% plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
% fprintf('Saving to %s...\n', plotFileName);
% export_fig(plotFileName, '-nocrop');
% 
% %% plot power in cue-target delay vPul P3 vs P1
% channelCond = isInVPulvinar;
% nChannel = sum(channelCond);
% cueTargetDelayRelPowerLF = [cueTargetDelayPowerP3LF(channelCond,:); cueTargetDelayPowerP1LF(channelCond,:)] ./ ...
%         [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
% cueTargetDelayRelPowerHF = [cueTargetDelayPowerP3HF(channelCond,:); cueTargetDelayPowerP1HF(channelCond,:)] ./ ...
%         [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
% p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];
% 
% grcPlotLfpPower2(cueTargetDelayRelPowerLF, cueTargetDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
%         'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
%         'yBounds', [relPowYBounds; relPowYBounds], ...
%         'cols', [p3Col; p1Col], ...
%         'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
%         'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
%         'titleText', '', ...
%         'doDB', 1);
% 
% plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
% fprintf('Saving to %s...\n', plotFileName);
% export_fig(plotFileName, '-nocrop');

%% plot power in array response allPul P3 vs P1
channelCond = isInPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelPowerLF = [arrayResponseHoldPowerP3LF(channelCond,:); arrayResponseHoldPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
arrayResponseHoldRelPowerHF = [arrayResponseHoldPowerP3HF(channelCond,:); arrayResponseHoldPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(arrayResponseHoldRelPowerLF, arrayResponseHoldRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', '', ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldRelPowerDB-pul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

% %% plot power in array response dPul P3 vs P1
% channelCond = isInDPulvinar;
% nChannel = sum(channelCond);
% arrayResponseHoldRelPowerLF = [arrayResponseHoldPowerP3LF(channelCond,:); arrayResponseHoldPowerP1LF(channelCond,:)] ./ ...
%         [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
% arrayResponseHoldRelPowerHF = [arrayResponseHoldPowerP3HF(channelCond,:); arrayResponseHoldPowerP1HF(channelCond,:)] ./ ...
%         [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
% p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];
% 
% grcPlotLfpPower2(arrayResponseHoldRelPowerLF, arrayResponseHoldRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
%         'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
%         'yBounds', [relPowYBounds; relPowYBounds], ...
%         'cols', [p3Col; p1Col], ...
%         'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
%         'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
%         'titleText', '', ...
%         'doDB', 1);
% 
% plotFileName = sprintf('%s/allSessions-arrayResponseHoldRelPowerDB-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
% fprintf('Saving to %s...\n', plotFileName);
% export_fig(plotFileName, '-nocrop');
% 
% %% plot power in array response vPul P3 vs P1
% channelCond = isInVPulvinar;
% nChannel = sum(channelCond);
% arrayResponseHoldRelPowerLF = [arrayResponseHoldPowerP3LF(channelCond,:); arrayResponseHoldPowerP1LF(channelCond,:)] ./ ...
%         [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
% arrayResponseHoldRelPowerHF = [arrayResponseHoldPowerP3HF(channelCond,:); arrayResponseHoldPowerP1HF(channelCond,:)] ./ ...
%         [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
% p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];
% 
% grcPlotLfpPower2(arrayResponseHoldRelPowerLF, arrayResponseHoldRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
%         'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
%         'yBounds', [relPowYBounds; relPowYBounds], ...
%         'cols', [p3Col; p1Col], ...
%         'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
%         'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
%         'titleText', '', ...
%         'doDB', 1);
% 
% plotFileName = sprintf('%s/allSessions-arrayResponseHoldRelPowerDB-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
% fprintf('Saving to %s...\n', plotFileName);
% export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC dPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelSFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelayRelSFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC vPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelSFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelayRelSFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot array response SFC dPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelSFCLF = [arrayResponseHoldSFCP3LF(channelCond,:); arrayResponseHoldSFCP1LF(channelCond,:)];
arrayResponseHoldRelSFCHF = [arrayResponseHoldSFCP3HF(channelCond,:); arrayResponseHoldSFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldSFC-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot array response SFC vPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelSFCLF = [arrayResponseHoldSFCP3LF(channelCond,:); arrayResponseHoldSFCP1LF(channelCond,:)];
arrayResponseHoldRelSFCHF = [arrayResponseHoldSFCP3HF(channelCond,:); arrayResponseHoldSFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldSFC-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs dPul spike - vPul field coherence across conditions and periods
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF, 1);
sfcAllLF = [cueTargetDelaySFCDPulSpikeVPulFieldP3LF; cueTargetDelaySFCDPulSpikeVPulFieldP1LF];
sfcAllHF = [cueTargetDelaySFCDPulSpikeVPulFieldP3HF; cueTargetDelaySFCDPulSpikeVPulFieldP1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-ctDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF, 1);
sfcAllLF = [cueTargetDelaySFCVPulSpikeDPulFieldP3LF; cueTargetDelaySFCVPulSpikeDPulFieldP1LF];
sfcAllHF = [cueTargetDelaySFCVPulSpikeDPulFieldP3HF; cueTargetDelaySFCVPulSpikeDPulFieldP1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-ctDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs dPul spike - vPul field coherence across conditions and periods
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF, 1);
sfcAllLF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3LF; arrayResponseHoldSFCDPulSpikeVPulFieldP1LF];
sfcAllHF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3HF; arrayResponseHoldSFCDPulSpikeVPulFieldP1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-arrayResponseHoldSFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF, 1);
sfcAllLF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3LF; arrayResponseHoldSFCVPulSpikeDPulFieldP1LF];
sfcAllHF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3HF; arrayResponseHoldSFCVPulSpikeDPulFieldP1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-arrayResponseHoldSFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%%
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF, 1);
sfcAllLF = [cueTargetDelaySFCDPulSpikeVPulFieldP3LF; cueTargetDelaySFCDPulSpikeVPulFieldP1LF];
sfcAllHF = [cueTargetDelaySFCDPulSpikeVPulFieldP3HF; cueTargetDelaySFCDPulSpikeVPulFieldP1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

%%
sfcMeanBaselineDPulSpikeVPulField = mean(baselineSFCDPulSpikeVPulFieldHF, 2);
% sfcMeanHFP3 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF, 2);
% sfcMeanHFP1 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF, 2);

sIndDPulSpikeVPulField = sfcMeanBaselineDPulSpikeVPulField > 0.01;
nSession = sum(sIndDPulSpikeVPulField)
% dPulSpikeVPulFieldNames(sIndDPulSpikeVPulField)

fInd = fAxisHF >= 35 & fAxisHF <= 45;
sfcMeanLowGammaP3 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,fInd), 2);
sfcMeanLowGammaP1 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,fInd), 2);

figure_tr_inch(5, 5);
hold on;
plot([0 1], [0 1], '-', 'Color', 0.3*ones(3, 1));
plot(sfcMeanLowGammaP3, sfcMeanLowGammaP1, 'k.', 'MarkerSize', 20);
xlim([0 0.1]);
ylim([0 0.1]);

[h,p] = ttest(sfcMeanLowGammaP3, sfcMeanLowGammaP1)

figure;
hold on;
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,fInd), 1))
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,fInd), 1))
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,fInd), 1) - std(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,fInd), 0, 1) / sqrt(nSession));
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,fInd), 1) + std(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,fInd), 0, 1) / sqrt(nSession));

%%
sfcMeanBaselineVPulSpikeDPulField = mean(baselineSFCVPulSpikeDPulFieldHF, 2);
% sfcMeanHFP3 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF, 2);
% sfcMeanHFP1 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF, 2);

sIndVPulSpikeDPulField = sfcMeanBaselineVPulSpikeDPulField > 0.01;
nSession = sum(sIndVPulSpikeDPulField)
% vPulSpikeDPulFieldNames(sIndVPulSpikeDPulField)

fInd = fAxisHF >= 60 & fAxisHF <= 70;
sfcMeanHighGammaP3 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,fInd), 2);
sfcMeanHighGammaP1 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,fInd), 2);

figure_tr_inch(5, 5);
hold on;
plot([0 1], [0 1], '-', 'Color', 0.3*ones(3, 1));
plot(sfcMeanHighGammaP3, sfcMeanHighGammaP1, 'k.', 'MarkerSize', 20);
xlim([0 0.1]);
ylim([0 0.1]);

[h,p] = ttest(sfcMeanHighGammaP3, sfcMeanHighGammaP1)

figure;
hold on;
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,fInd), 1))
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,fInd), 1))
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,fInd), 1) - std(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,fInd), 0, 1) / sqrt(nSession));
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,fInd), 1) + std(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,fInd), 0, 1) / sqrt(nSession));


%%
sfcCTDelayYBoundsLF = [0.035 0.07];
sfcCTDelayYBoundsHF = [0.015 0.035];

%%
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(sIndDPulSpikeVPulField,:), 1);
sfcAllLF = [cueTargetDelaySFCDPulSpikeVPulFieldP3LF(sIndDPulSpikeVPulField,:); cueTargetDelaySFCDPulSpikeVPulFieldP1LF(sIndDPulSpikeVPulField,:)];
sfcAllHF = [cueTargetDelaySFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,:); cueTargetDelaySFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-ctDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(sIndVPulSpikeDPulField,:), 1);
sfcAllLF = [cueTargetDelaySFCVPulSpikeDPulFieldP3LF(sIndVPulSpikeDPulField,:); cueTargetDelaySFCVPulSpikeDPulFieldP1LF(sIndVPulSpikeDPulField,:)];
sfcAllHF = [cueTargetDelaySFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,:); cueTargetDelaySFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-ctDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs dPul spike - vPul field coherence across conditions and periods
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(sIndDPulSpikeVPulField,:), 1);
sfcAllLF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3LF(sIndDPulSpikeVPulField,:); arrayResponseHoldSFCDPulSpikeVPulFieldP1LF(sIndDPulSpikeVPulField,:)];
sfcAllHF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3HF(sIndDPulSpikeVPulField,:); arrayResponseHoldSFCDPulSpikeVPulFieldP1HF(sIndDPulSpikeVPulField,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-arrayResponseHoldSFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(sIndVPulSpikeDPulField,:), 1);
sfcAllLF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3LF(sIndVPulSpikeDPulField,:); arrayResponseHoldSFCVPulSpikeDPulFieldP1LF(sIndVPulSpikeDPulField,:)];
sfcAllHF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(sIndVPulSpikeDPulField,:); arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(sIndVPulSpikeDPulField,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end



grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-arrayResponseHoldSFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');










%%
p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

relPowYBounds = [-1.5 0.2];
sfcCTDelayYBoundsLF = [0.04 0.1];
sfcCTDelayYBoundsHF = [0.015 0.052];
sfcArrayRespYBoundsLF = [0.04 0.1];
sfcArrayRespYBoundsHF = [0.015 0.052];
sfcTDDelayYBoundsLF = [0.07 0.18];
sfcTDDelayYBoundsHF = [0.02 0.08];

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInDPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInDPulvinar & goodPairs';
nChannel = sum(channelCond);
cueTargetDelayRelSFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelayRelSFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-dPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInVPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInVPulvinar & goodPairs';
nChannel = sum(channelCond);
cueTargetDelayRelSFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelayRelSFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-vPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInDPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInDPulvinar & goodPairs';
nChannel = sum(channelCond);
arrayResponseHoldRelSFCLF = [arrayResponseHoldSFCP3LF(channelCond,:); arrayResponseHoldSFCP1LF(channelCond,:)];
arrayResponseHoldRelSFCHF = [arrayResponseHoldSFCP3HF(channelCond,:); arrayResponseHoldSFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcArrayRespYBoundsLF; sfcArrayRespYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldSFC-dPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInVPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInVPulvinar & goodPairs';
nChannel = sum(channelCond);
arrayResponseHoldRelSFCLF = [arrayResponseHoldSFCP3LF(channelCond,:); arrayResponseHoldSFCP1LF(channelCond,:)];
arrayResponseHoldRelSFCHF = [arrayResponseHoldSFCP3HF(channelCond,:); arrayResponseHoldSFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcArrayRespYBoundsLF; sfcArrayRespYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldSFC-vPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInDPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInDPulvinar & goodPairs';
nChannel = sum(channelCond);
targetDimDelaySFCLF = [targetDimDelaySFCP3LF(channelCond,:); targetDimDelaySFCP1LF(channelCond,:)];
targetDimDelaySFCHF = [targetDimDelaySFCP3HF(channelCond,:); targetDimDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(targetDimDelaySFCLF, targetDimDelaySFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcTDDelayYBoundsLF; sfcTDDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-tdDelaySFC-dPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot SFC for pairs with nonflat baseline SFC
% goodPairs = true(1, size(isInVPulvinar, 1));
goodPairs = ~any(isnan(baselineSFCLF'),1) & ~(all(baselineSFCLF' < 0.05) | all(baselineSFCHF' < 0.03));

channelCond = isInVPulvinar & goodPairs';
nChannel = sum(channelCond);
targetDimDelaySFCLF = [targetDimDelaySFCP3LF(channelCond,:); targetDimDelaySFCP1LF(channelCond,:)];
targetDimDelaySFCHF = [targetDimDelaySFCP3HF(channelCond,:); targetDimDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

grcPlotLfpPower2(targetDimDelaySFCLF, targetDimDelaySFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcTDDelayYBoundsLF; sfcTDDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-tdDelaySFC-vPul-P3vsP1-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');



%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.02));

nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(goodPairs,:), 1);
sfcAllLF = [cueTargetDelaySFCDPulSpikeVPulFieldP3LF(goodPairs,:); cueTargetDelaySFCDPulSpikeVPulFieldP1LF(goodPairs,:)];
sfcAllHF = [cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,:); cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-ctDelaySFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

fIndHF = fAxisHF >= 30 & fAxisHF <= 50;
sfcP3Sub = mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fIndHF), 2);
sfcP1Sub = mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fIndHF), 2);
[h,p] = ttest(sfcP3Sub, sfcP1Sub)


%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(goodPairs,:), 1);
sfcAllLF = [cueTargetDelaySFCVPulSpikeDPulFieldP3LF(goodPairs,:); cueTargetDelaySFCVPulSpikeDPulFieldP1LF(goodPairs,:)];
sfcAllHF = [cueTargetDelaySFCVPulSpikeDPulFieldP3HF(goodPairs,:); cueTargetDelaySFCVPulSpikeDPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcCTDelayYBoundsLF; sfcCTDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-ctDelaySFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(goodPairs,:), 1);
sfcAllLF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3LF(goodPairs,:); arrayResponseHoldSFCDPulSpikeVPulFieldP1LF(goodPairs,:)];
sfcAllHF = [arrayResponseHoldSFCDPulSpikeVPulFieldP3HF(goodPairs,:); arrayResponseHoldSFCDPulSpikeVPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcArrayRespYBoundsLF; sfcArrayRespYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-arrayResponseHoldSFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(goodPairs,:), 1);
sfcAllLF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3LF(goodPairs,:); arrayResponseHoldSFCVPulSpikeDPulFieldP1LF(goodPairs,:)];
sfcAllHF = [arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,:); arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcArrayRespYBoundsLF; sfcArrayRespYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-arrayResponseHoldSFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(goodPairs,:), 1);
sfcAllLF = [targetDimDelaySFCDPulSpikeVPulFieldP3LF(goodPairs,:); targetDimDelaySFCDPulSpikeVPulFieldP1LF(goodPairs,:)];
sfcAllHF = [targetDimDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,:); targetDimDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcTDDelayYBoundsLF; sfcTDDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-tdDelaySFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(goodPairs,:), 1);
sfcAllLF = [targetDimDelaySFCVPulSpikeDPulFieldP3LF(goodPairs,:); targetDimDelaySFCVPulSpikeDPulFieldP1LF(goodPairs,:)];
sfcAllHF = [targetDimDelaySFCVPulSpikeDPulFieldP3HF(goodPairs,:); targetDimDelaySFCVPulSpikeDPulFieldP1HF(goodPairs,:)];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

grcPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass(1) 80], ...
        'yBounds', [sfcTDDelayYBoundsLF; sfcTDDelayYBoundsHF], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'Attend-RF', 'Attend-Away'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', '', ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-tdDelaySFC-baselineRefined-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%%
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF(goodPairs,:), 1)
% dPulSpikeVPulFieldNames(goodPairs)

fInd = fAxisHF >= 35 & fAxisHF <= 45;
sfcCTMeanLowGammaP3 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fInd), 2);
sfcCTMeanLowGammaP1 = mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fInd), 2);

figure;
hold on;
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fInd), 1))
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fInd), 1))
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fInd), 1) - std(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fInd), 0, 1) / sqrt(nSubPairs));
plot(fAxisHF(fInd), mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fInd), 1) + std(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fInd), 0, 1) / sqrt(nSubPairs));

[~,p] = ttest(sfcCTMeanLowGammaP3, sfcCTMeanLowGammaP1)

%%
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF(goodPairs,:), 1);
% vPulSpikeDPulFieldNames(goodPairs)

fInd = fAxisHF >= 60 & fAxisHF <= 70;
sfcARMeanHighGammaP3 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,fInd), 2);
sfcARMeanHighGammaP1 = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,fInd), 2);

figure;
hold on;
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,fInd), 1))
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,fInd), 1))
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,fInd), 1) - std(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,fInd), 0, 1) / sqrt(nSubPairs));
plot(fAxisHF(fInd), mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,fInd), 1) + std(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,fInd), 0, 1) / sqrt(nSubPairs));

[~,p] = ttest(sfcARMeanHighGammaP3, sfcARMeanHighGammaP1)
