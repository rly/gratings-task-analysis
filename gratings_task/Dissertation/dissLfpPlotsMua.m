%%
clear;
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\LFP_GRATINGS_SUMMARY';
ref = 'CAR';
v = 12;

saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-%s-v%d.mat', outputDir, ref, v);
load(saveFileName);
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\Dissertation';

%%
p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

relPowYBounds = [-1.5 0.2];
sfcCTDelayYBoundsLF = [0.04 0.09];
sfcCTDelayYBoundsHF = [0.015 0.045];
sfcArrayRespYBoundsLF = [0.04 0.09];
sfcArrayRespYBoundsHF = [0.015 0.045];
sfcTDDelayYBoundsLF = [0.07 0.18];
sfcTDDelayYBoundsHF = [0.02 0.08];

baselineQuantileThresh = 0.8;

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

threshBaselineLF = quantile(mean(baselineSFCDPulSpikeVPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCDPulSpikeVPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ...
        any(baselineSFCDPulSpikeVPulFieldLF' > threshBaselineLF) & ...
        any(baselineSFCDPulSpikeVPulFieldHF' > threshBaselineHF);

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
    
plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-ctDelaySFC-baselineRefined-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

fIndHF = fAxisHF >= 30 & fAxisHF <= 50;
sfcP3Sub = mean(cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,fIndHF), 2);
sfcP1Sub = mean(cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,fIndHF), 2);
[h,p] = ttest(sfcP3Sub, sfcP1Sub)
p = signrank(sfcP3Sub, sfcP1Sub)

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

threshBaselineLF = quantile(mean(baselineSFCVPulSpikeDPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCVPulSpikeDPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ...
        any(baselineSFCVPulSpikeDPulFieldLF' > threshBaselineLF) & ...
        any(baselineSFCVPulSpikeDPulFieldHF' > threshBaselineHF);
    
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
    
plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-ctDelaySFC-baselineRefined-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

fIndHF = fAxisHF >= 50 & fAxisHF <= 70;
sfcP3Sub = mean(cueTargetDelaySFCVPulSpikeDPulFieldP3HF(goodPairs,fIndHF), 2);
sfcP1Sub = mean(cueTargetDelaySFCVPulSpikeDPulFieldP1HF(goodPairs,fIndHF), 2);
p = signrank(sfcP3Sub, sfcP1Sub)

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

threshBaselineLF = quantile(mean(baselineSFCDPulSpikeVPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCDPulSpikeVPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) & ...
        any(baselineSFCDPulSpikeVPulFieldLF' > threshBaselineLF) & ...
        any(baselineSFCDPulSpikeVPulFieldHF' > threshBaselineHF);

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
    
plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-arrayResponseHoldSFC-baselineRefined-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

fIndHF = fAxisHF >= 30 & fAxisHF <= 50;
sfcP3Sub = mean(arrayResponseHoldSFCDPulSpikeVPulFieldP3HF(goodPairs,fIndHF), 2);
sfcP1Sub = mean(arrayResponseHoldSFCDPulSpikeVPulFieldP1HF(goodPairs,fIndHF), 2);
p = signrank(sfcP3Sub, sfcP1Sub)

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

threshBaselineLF = quantile(mean(baselineSFCVPulSpikeDPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCVPulSpikeDPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) & ...
        any(baselineSFCVPulSpikeDPulFieldLF' > threshBaselineLF) & ...
        any(baselineSFCVPulSpikeDPulFieldHF' > threshBaselineHF);
    
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
    
plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-arrayResponseHoldSFC-baselineRefined-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

fIndHF = fAxisHF >= 50 & fAxisHF <= 70;
sfcP3Sub = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(goodPairs,fIndHF), 2);
sfcP1Sub = mean(arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(goodPairs,fIndHF), 2);
p = signrank(sfcP3Sub, sfcP1Sub)
