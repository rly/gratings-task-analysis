%%
clear;
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\LFP_GRATINGS_SUMMARY';
ref = 'CAR';
v = 13;

saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-%s-v%d.mat', outputDir, ref, v);
load(saveFileName);
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\Dissertation';

%%
cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];
baselineCol = [0 0.9 0];

if strcmp(ref, 'BIP')
    relPowYBounds = [-0.55 0.1];
    sfcCTDelayYBoundsLF = [0.035 0.065];
    sfcCTDelayYBoundsHF = [0.0125 0.03];
elseif strcmp(ref, 'CAR')
    relPowYBounds = [-1.5 0.2];
    sfcCTDelayYBoundsLF = [0.035 0.07];
    sfcCTDelayYBoundsHF = [0.015 0.03];
end

baselineQuantileThresh = 0.8;

%% plot cue-target delay SFC dPul P3 vs P1 (relative to baseline is NOISY)
    
ixRemove = any(isnan(cueTargetDelaySFCDPulSpikeDPulFieldP3LF), 2) | any(isnan(cueTargetDelaySFCDPulSpikeDPulFieldP3HF), 2) | ...
        any(isnan(cueTargetDelaySFCDPulSpikeDPulFieldP1LF), 2) | any(isnan(cueTargetDelaySFCDPulSpikeDPulFieldP1HF), 2) | ...
        all(cueTargetDelaySFCDPulSpikeDPulFieldP3LF < 0.025, 2) | all(cueTargetDelaySFCDPulSpikeDPulFieldP3HF < 0.025, 2) | ...
        all(cueTargetDelaySFCDPulSpikeDPulFieldP1LF < 0.025, 2) | all(cueTargetDelaySFCDPulSpikeDPulFieldP1HF < 0.025, 2) | ...
        any(cueTargetDelaySFCDPulSpikeDPulFieldP3LF > 0.5, 2) | any(cueTargetDelaySFCDPulSpikeDPulFieldP3HF > 0.3, 2) | ...
        any(cueTargetDelaySFCDPulSpikeDPulFieldP1LF > 0.5, 2) | any(cueTargetDelaySFCDPulSpikeDPulFieldP1HF > 0.3, 2);
p3LF = cueTargetDelaySFCDPulSpikeDPulFieldP3LF(~ixRemove,:);
p3HF = cueTargetDelaySFCDPulSpikeDPulFieldP3HF(~ixRemove,:);
p1LF = cueTargetDelaySFCDPulSpikeDPulFieldP1LF(~ixRemove,:);
p1HF = cueTargetDelaySFCDPulSpikeDPulFieldP1HF(~ixRemove,:);

nSubPairs = size(p3LF, 1);
sfcAllLF = [p3LF; p1LF];
sfcAllHF = [p3HF; p1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.12];
dissPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'CT Delay P3', 'CT Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('dPul Spike - dPul Field Coherence (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-dPulField-SFC-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC vPul P3 vs P1 (relative to baseline is NOISY)
% clean data

% ixRemove = any(isnan(cueTargetDelaySFCVPulSpikeVPulFieldP3LF), 2) | any(isnan(cueTargetDelaySFCVPulSpikeVPulFieldP3HF), 2) | ...
%         any(isnan(cueTargetDelaySFCVPulSpikeVPulFieldP1LF), 2) | any(isnan(cueTargetDelaySFCVPulSpikeVPulFieldP1HF), 2) | ...
%         all(cueTargetDelaySFCVPulSpikeVPulFieldP3LF < 0.03, 2) | all(cueTargetDelaySFCVPulSpikeVPulFieldP3HF < 0.03, 2) | ...
%         all(cueTargetDelaySFCVPulSpikeVPulFieldP1LF < 0.03, 2) | all(cueTargetDelaySFCVPulSpikeVPulFieldP1HF < 0.03, 2) | ...
%         any(cueTargetDelaySFCVPulSpikeVPulFieldP3LF > 0.5, 2) | any(cueTargetDelaySFCVPulSpikeVPulFieldP3HF > 0.3, 2) | ...
%         any(cueTargetDelaySFCVPulSpikeVPulFieldP1LF > 0.5, 2) | any(cueTargetDelaySFCVPulSpikeVPulFieldP1HF > 0.3, 2);
p3LF = cueTargetDelaySFCVPulSpikeVPulFieldP3LF(~ixRemove,:);
p3HF = cueTargetDelaySFCVPulSpikeVPulFieldP3HF(~ixRemove,:);
p1LF = cueTargetDelaySFCVPulSpikeVPulFieldP1LF(~ixRemove,:);
p1HF = cueTargetDelaySFCVPulSpikeVPulFieldP1HF(~ixRemove,:);

nSubPairs = size(p3LF, 1);
sfcAllLF = [p3LF; p1LF];
sfcAllHF = [p3HF; p1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.12];
baselineCol = [0 0.9 0];
dissPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'CT Delay P3', 'CT Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('vPul Spike - vPul Field Coherence (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-vPulField-SFC-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay dPul spike - vPul LFP P3 vs P1 (relative to baseline is NOISY)
threshBaselineLF = quantile(mean(baselineSFCDPulSpikeVPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCDPulSpikeVPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCDPulSpikeVPulFieldLF, 1), 1);
ixRemove = any(isnan(baselineSFCDPulSpikeVPulFieldLF'),1) | ...
        all(baselineSFCDPulSpikeVPulFieldLF' < threshBaselineLF) | ...
        all(baselineSFCDPulSpikeVPulFieldHF' < threshBaselineHF);

% ixRemove = any(isnan(baselineSFCDPulSpikeVPulFieldLF), 2) | any(isnan(baselineSFCDPulSpikeVPulFieldHF), 2) | ...
%         all(baselineSFCDPulSpikeVPulFieldLF < 0.02, 2) | all(baselineSFCDPulSpikeVPulFieldHF < 0.01, 2);
p3LF = cueTargetDelaySFCDPulSpikeVPulFieldP3LF(~ixRemove,:);
p3HF = cueTargetDelaySFCDPulSpikeVPulFieldP3HF(~ixRemove,:);
p1LF = cueTargetDelaySFCDPulSpikeVPulFieldP1LF(~ixRemove,:);
p1HF = cueTargetDelaySFCDPulSpikeVPulFieldP1HF(~ixRemove,:);

nSubPairs = size(p3LF, 1);
sfcAllLF = [p3LF; p1LF];
sfcAllHF = [p3HF; p1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.12];
dissPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; 25 80], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'CT Delay P3', 'CT Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('dPul Spike - vPul Field Coherence (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-SFC-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC vPul spike - dPul LFP P3 vs P1 (relative to baseline is NOISY)
% clean data
threshBaselineLF = quantile(mean(baselineSFCVPulSpikeDPulFieldLF, 2), baselineQuantileThresh);
threshBaselineHF = quantile(mean(baselineSFCVPulSpikeDPulFieldHF, 2), baselineQuantileThresh);

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
ixRemove = any(isnan(baselineSFCVPulSpikeDPulFieldLF'),1) | ...
        all(baselineSFCVPulSpikeDPulFieldLF' < threshBaselineLF) | ...
        all(baselineSFCVPulSpikeDPulFieldHF' < threshBaselineHF);

p3LF = cueTargetDelaySFCVPulSpikeDPulFieldP3LF(~ixRemove,:);
p3HF = cueTargetDelaySFCVPulSpikeDPulFieldP3HF(~ixRemove,:);
p1LF = cueTargetDelaySFCVPulSpikeDPulFieldP1LF(~ixRemove,:);
p1HF = cueTargetDelaySFCVPulSpikeDPulFieldP1HF(~ixRemove,:);

nSubPairs = size(p3LF, 1);
sfcAllLF = [p3LF; p1LF];
sfcAllHF = [p3HF; p1HF];
condLogical = false(nSubPairs * 2, 2);
for i = 1:2
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.12];
baselineCol = [0 0.9 0];
dissPlotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; 25 80], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'CT Delay P3', 'CT Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('vPul Spike - dPul Field Coherence (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-SFC-%s-diss-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
