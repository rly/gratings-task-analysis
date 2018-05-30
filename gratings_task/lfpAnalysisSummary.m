function lfpAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds, ref)

% clear;
% readDataLocally;
% sessionInds = 8;%1:23;

v = 12;

%%
recordingInfo = readRecordingInfo(recordingInfoFileName);

outputDir = sprintf('%s/%s/', processedDataRootDir, 'LFP_GRATINGS_SUMMARY');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

nSessions = numel(sessionInds);
nLfpsApprox = nSessions * 16; % should be equal or an underestimate
lfpCount = 0;

lfpNames = cell(nLfpsApprox, 1);
isInDPulvinar = false(nLfpsApprox, 1);
isInVPulvinar = false(nLfpsApprox, 1);

fAxisLF = NaN;
nFAxisLF = 20;
baselinePowerLF = nan(nLfpsApprox, nFAxisLF);
cueResponsePowerLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerAllLocsLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerP3LF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerP1LF = nan(nLfpsApprox, nFAxisLF);
baselineSFCLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP3LF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP1LF = nan(nLfpsApprox, nFAxisLF);

fAxisHF = NaN;
nFAxisHF = 77;
baselinePowerHF = nan(nLfpsApprox, nFAxisHF);
cueResponsePowerHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerAllLocsHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP1HF = nan(nLfpsApprox, nFAxisHF);
baselineSFCHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP1HF = nan(nLfpsApprox, nFAxisHF);

baselineWindowOffset = [-0.25 0];
cueResponseOffset = [0 0.25];
cueTargetDelayOffset = [-0.25 0];

% chronux parameters
paramsLF.tapers = [1 1];
paramsLF.fpass = [5 25];
paramsLF.pad = 2;
paramsLF.Fs = 1000;
paramsLF.trialave = 1;

paramsHF.tapers = [2 3];
paramsHF.fpass = [25 100];
paramsHF.pad = 2;
paramsHF.Fs = 1000;
paramsHF.trialave = 1;

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');

%% session loop
for i = 1:nSessions
    % TODO
    % session 15, 'M20170329', 1-32
    % session 17, 'M20170329', 1-32
    % have abnormally high power around 40 Hz
    if i == 15 || i == 17
        continue;
    end
    sessionInd = sessionInds(i);
    R = recordingInfo(sessionInd);
    sessionName = R.sessionName;
    lfpChannelsToLoad = R.lfpChannelsToLoad;
    blockName = strjoin(R.blockNames(R.gratingsTask3DIndices), '-'); % temp, only 3D
    processedDataDir = sprintf('%s/%s/LFP_GRATINGS/', processedDataRootDir, sessionName);
    fileNamePrefix = sprintf('%s-ch%d-ch%d-%s', sessionName, lfpChannelsToLoad([1 end]), blockName);
    saveFileName = sprintf('%s/%s-evokedLfps-v%d.mat', processedDataDir, fileNamePrefix, v);
    fprintf('Loading file %s...\n', saveFileName);
    EL = load(saveFileName);
    
    for j = 1:numel(EL.channelInds)
        if strcmp(ref, 'BIP') && j == numel(EL.channelInds)
            continue;
        end
        fprintf('Processing channel %d...\n', EL.channelInds(j));
        
        lfpCount = lfpCount + 1;
        lfpNames{lfpCount,:} = sprintf('%s_FP%03d', sessionName, EL.channelInds(j));
        isInVPulvinar(lfpCount) = 0;
        isInDPulvinar(lfpCount) = 0;
        % TODO decide how to split up end channels for bipolar reference
        if ismember(EL.channelInds(j), R.vPulChannels)
            isInVPulvinar(lfpCount) = 1;
        end
        if ismember(EL.channelInds(j), R.dPulChannels)
            isInDPulvinar(lfpCount) = 1;
        end
        
        % read MUA from a nearby channel to account for MUA possibly
        % contributing to LFP on the same channel
        % use +2/-2 to be compatible with bipolar reference
        if j < numel(EL.channelInds)-1
            muaNearbyChannelTs = EL.allMUAStructs{j+2}.ts;
        else
            muaNearbyChannelTs = EL.allMUAStructs{j-2}.ts;
        end
        
        baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);
        cueResponseInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, cueResponseOffset);
        cueTargetDelayInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, cueTargetDelayOffset);

        % compute power in baseline - should probably be bipolar reference
        % all locations together
        if strcmp(ref, 'RAW')
            cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:))';
        elseif strcmp(ref, 'BIP')
            cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j+1,:,:) - EL.cueOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j+1,:,:) - EL.arrayOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j+1,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j+1,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:))';
        end
        % ignore array shapes because we are looking at pre-array delay
        % period
        
        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        [baselinePowerLF(lfpCount,:),fAxisLF] = mtspectrumc(preCueBaselineLfps, paramsLF);
        [baselinePowerHF(lfpCount,:),fAxisHF] = mtspectrumc(preCueBaselineLfps, paramsHF);
        
%         cueResponseLfps = cueOnsetLfpCurrent(cueResponseInd,:); % currently unused
%         [cueResponsePowerLF(lfpCount,:),fAxisLF] = mtspectrumc(cueResponseLfps, paramsLF);
%         [cueResponsePowerHF(lfpCount,:),fAxisHF] = mtspectrumc(cueResponseLfps, paramsHF);
        
%         cueTargetDelayLfps = arrayOnsetLfpCurrent(cueTargetDelayInd,:); % currently unused
%         [cueTargetDelayPowerAllLocsLF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfps, paramsLF);
%         [cueTargetDelayPowerAllLocsHF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfps, paramsHF);
        
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        [cueTargetDelayPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP3, paramsLF);
        [cueTargetDelayPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP3, paramsHF);
        
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        [cueTargetDelayPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP1, paramsLF);
        [cueTargetDelayPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP1, paramsHF);
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
        [baselineSFCLF(lfpCount,:),phi,S12,S1,S2,fAxisLF] = coherencycpt(preCueBaselineLfps, alignedSpikeTs, paramsLF);
        [baselineSFCHF(lfpCount,:),phi,S12,S1,S2,fAxisHF] = coherencycpt(preCueBaselineLfps, alignedSpikeTs, paramsHF);
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
        [cueTargetDelaySFCP3LF(lfpCount,:),phi,S12,S1,S2,fAxisLF] = coherencycpt(cueTargetDelayLfpsP3, alignedSpikeTs, paramsLF);
        [cueTargetDelaySFCP3HF(lfpCount,:),phi,S12,S1,S2,fAxisHF] = coherencycpt(cueTargetDelayLfpsP3, alignedSpikeTs, paramsHF);
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
        [cueTargetDelaySFCP1LF(lfpCount,:),phi,S12,S1,S2,fAxisLF] = coherencycpt(cueTargetDelayLfpsP1, alignedSpikeTs, paramsLF);
        [cueTargetDelaySFCP1HF(lfpCount,:),phi,S12,S1,S2,fAxisHF] = coherencycpt(cueTargetDelayLfpsP1, alignedSpikeTs, paramsHF);
    end
    
    clear EL;
end

%%
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-v%d.mat', outputDir, v);
save(saveFileName);

%%
outputDirOrig = outputDir;
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-v%d.mat', outputDir, v);
load(saveFileName);
outputDir = outputDirOrig;

%%
cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

maxSDsPower = 5;

%% plot baseline power LF
fAxis = fAxisLF;
baselinePowerDPul = (baselinePowerLF(isInDPulvinar,:));
baselinePowerVPul = (baselinePowerLF(isInVPulvinar,:));
xBounds = paramsLF.fpass;
yBounds = [0 0.01];

plotLfpPower(baselinePowerDPul, baselinePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Power');
title('Pre-Cue Baseline Power');

plotFileName = sprintf('%s/allSessions-baselinePower-LF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot baseline power HF
fAxis = fAxisHF;
baselinePowerDPul = (baselinePowerHF(isInDPulvinar,:));
baselinePowerVPul = (baselinePowerHF(isInVPulvinar,:));
xBounds = paramsHF.fpass;
yBounds = [0 0.003];

% outlier removal
baselinePowerDPul = cleanRowsBySDsFromMean(baselinePowerDPul, maxSDsPower);
baselinePowerVPul = cleanRowsBySDsFromMean(baselinePowerVPul, maxSDsPower);

plotLfpPower(baselinePowerDPul, baselinePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

title('Pre-Cue Baseline Power');

plotFileName = sprintf('%s/allSessions-baselinePower-HF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline LF
fAxis = fAxisLF;
cueTargetDelayRelativePowerDPul = (((cueTargetDelayPowerP3LF(isInDPulvinar,:))) ./ ((baselinePowerLF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerVPul = (((cueTargetDelayPowerP3LF(isInVPulvinar,:))) ./ ((baselinePowerLF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsLF.fpass;
yBounds = [-25 0];

% outlier removal
cueTargetDelayRelativePowerDPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerDPul, maxSDsPower);
cueTargetDelayRelativePowerVPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerVPul, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerDPul, cueTargetDelayRelativePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

title('Cue-Target Delay Power P3 dPul vs vPul');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-LF-P3-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline HF
fAxis = fAxisHF;
cueTargetDelayRelativePowerDPul = (((cueTargetDelayPowerP3HF(isInDPulvinar,:))) ./ ((baselinePowerHF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerVPul = (((cueTargetDelayPowerP3HF(isInVPulvinar,:))) ./ ((baselinePowerHF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsHF.fpass;
yBounds = [-10 5];

% outlier removal
cueTargetDelayRelativePowerDPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerDPul, maxSDsPower);
cueTargetDelayRelativePowerVPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerVPul, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerDPul, cueTargetDelayRelativePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

title('Cue-Target Delay Power P3 dPul vs vPul');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-HF-P3-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline LF
fAxis = fAxisLF;
cueTargetDelayRelativePowerDPul = (((cueTargetDelayPowerP1LF(isInDPulvinar,:))) ./ ((baselinePowerLF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerVPul = (((cueTargetDelayPowerP1LF(isInVPulvinar,:))) ./ ((baselinePowerLF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsLF.fpass;
yBounds = [-25 0];

% outlier removal
cueTargetDelayRelativePowerDPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerDPul, maxSDsPower);
cueTargetDelayRelativePowerVPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerVPul, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerDPul, cueTargetDelayRelativePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

title('Cue-Target Delay Power P3 dPul vs vPul');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-LF-P1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline HF
fAxis = fAxisHF;
cueTargetDelayRelativePowerDPul = (((cueTargetDelayPowerP1HF(isInDPulvinar,:))) ./ ((baselinePowerHF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerVPul = (((cueTargetDelayPowerP1HF(isInVPulvinar,:))) ./ ((baselinePowerHF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsHF.fpass;
yBounds = [-10 5];

% outlier removal
cueTargetDelayRelativePowerDPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerDPul, maxSDsPower);
cueTargetDelayRelativePowerVPul = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerVPul, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerDPul, cueTargetDelayRelativePowerVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

title('Cue-Target Delay Power P3 dPul vs vPul');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-HF-P1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
fAxis = fAxisLF;
cueTargetDelayRelativePowerP3 = (((cueTargetDelayPowerP3LF(isInDPulvinar,:))) ./ ((baselinePowerLF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerP1 = (((cueTargetDelayPowerP1LF(isInDPulvinar,:))) ./ ((baselinePowerLF(isInDPulvinar,:))) - 1) * 100;
xBounds = paramsLF.fpass;
yBounds = [-25 0];

% outlier removal
cueTargetDelayRelativePowerP3 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP3, maxSDsPower);
cueTargetDelayRelativePowerP1 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP1, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerP3, cueTargetDelayRelativePowerP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')
title('Cue-Target Delay Power dPul P3 vs P1');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-dPul-P3vsP1-LF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
fAxis = fAxisHF;
cueTargetDelayRelativePowerP3 = (((cueTargetDelayPowerP3HF(isInDPulvinar,:))) ./ ((baselinePowerHF(isInDPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerP1 = (((cueTargetDelayPowerP1HF(isInDPulvinar,:))) ./ ((baselinePowerHF(isInDPulvinar,:))) - 1) * 100;
xBounds = paramsHF.fpass;
yBounds = [-10 5];

% outlier removal
cueTargetDelayRelativePowerP3 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP3, maxSDsPower);
cueTargetDelayRelativePowerP1 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP1, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerP3, cueTargetDelayRelativePowerP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')
title('Cue-Target Delay Power dPul P3 vs P1');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-dPul-P3vsP1-HF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
fAxis = fAxisLF;
cueTargetDelayRelativePowerP3 = (((cueTargetDelayPowerP3LF(isInVPulvinar,:))) ./ ((baselinePowerLF(isInVPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerP1 = (((cueTargetDelayPowerP1LF(isInVPulvinar,:))) ./ ((baselinePowerLF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsLF.fpass;
yBounds = [-25 0];

% outlier removal
cueTargetDelayRelativePowerP3 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP3, maxSDsPower);
cueTargetDelayRelativePowerP1 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP1, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerP3, cueTargetDelayRelativePowerP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')
title('Cue-Target Delay Power vPul P3 vs P1');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-vPul-P3vsP1-LF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline
fAxis = fAxisHF;
cueTargetDelayRelativePowerP3 = (((cueTargetDelayPowerP3HF(isInVPulvinar,:))) ./ ((baselinePowerHF(isInVPulvinar,:))) - 1) * 100;
cueTargetDelayRelativePowerP1 = (((cueTargetDelayPowerP1HF(isInVPulvinar,:))) ./ ((baselinePowerHF(isInVPulvinar,:))) - 1) * 100;
xBounds = paramsHF.fpass;
yBounds = [-10 5];

% outlier removal
cueTargetDelayRelativePowerP3 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP3, maxSDsPower);
cueTargetDelayRelativePowerP1 = cleanRowsBySDsFromMean(cueTargetDelayRelativePowerP1, maxSDsPower);

plotLfpPower(cueTargetDelayRelativePowerP3, cueTargetDelayRelativePowerP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')
title('Cue-Target Delay Power vPul P3 vs P1');

plotFileName = sprintf('%s/allSessions-cueTargetDelayPower-vPul-P3vsP1-HF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot baseline SFC LF
fAxis = fAxisLF;
baselineSFCDPul = (baselineSFCLF(isInDPulvinar,:));
baselineSFCVPul = (baselineSFCLF(isInVPulvinar,:));
xBounds = paramsLF.fpass;
yBounds = [0 0.1];

plotLfpPower(baselineSFCDPul, baselineSFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Local Spike-Field Coherence');
title('Pre-Cue Baseline SFC (No Ref)');

plotFileName = sprintf('%s/allSessions-baselineSFC-LF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot baseline SFC LF
fAxis = fAxisHF;
baselineSFCDPul = (baselineSFCHF(isInDPulvinar,:));
baselineSFCVPul = (baselineSFCHF(isInVPulvinar,:));
xBounds = paramsHF.fpass;
yBounds = [0 0.1];

plotLfpPower(baselineSFCDPul, baselineSFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Local Spike-Field Coherence');
title('Pre-Cue Baseline SFC (No Ref)');

plotFileName = sprintf('%s/allSessions-baselineSFC-HF-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline LF
fAxis = fAxisLF;
cueTargetDelaySFCDPul = cueTargetDelaySFCP3LF(isInDPulvinar,:);
cueTargetDelaySFCVPul = cueTargetDelaySFCP3LF(isInVPulvinar,:);
xBounds = paramsLF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCDPul, cueTargetDelaySFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC P3 (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-LF-P3-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline HF
fAxis = fAxisHF;
cueTargetDelaySFCDPul = cueTargetDelaySFCP3HF(isInDPulvinar,:);
cueTargetDelaySFCVPul = cueTargetDelaySFCP3HF(isInVPulvinar,:);
xBounds = paramsHF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCDPul, cueTargetDelaySFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC P3 (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-HF-P3-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay LF
fAxis = fAxisLF;
cueTargetDelaySFCDPul = cueTargetDelaySFCP1LF(isInDPulvinar,:);
cueTargetDelaySFCVPul = cueTargetDelaySFCP1LF(isInVPulvinar,:);
xBounds = paramsLF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCDPul, cueTargetDelaySFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC P1 (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-LF-P1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay HF
fAxis = fAxisHF;
cueTargetDelaySFCDPul = cueTargetDelaySFCP1HF(isInDPulvinar,:);
cueTargetDelaySFCVPul = cueTargetDelaySFCP1HF(isInVPulvinar,:);
xBounds = paramsHF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCDPul, cueTargetDelaySFCVPul, fAxis, xBounds, yBounds, dPulCol, vPulCol, 'dPul', 'vPul')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC P1 (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-HF-P1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline LF
fAxis = fAxisLF;
cueTargetDelaySFCP3 = cueTargetDelaySFCP3LF(isInDPulvinar,:);
cueTargetDelaySFCP1 = cueTargetDelaySFCP1LF(isInDPulvinar,:);
xBounds = paramsLF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCP3, cueTargetDelaySFCP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC dPul (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-LF-dPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline LF
fAxis = fAxisHF;
cueTargetDelaySFCP3 = cueTargetDelaySFCP3HF(isInDPulvinar,:);
cueTargetDelaySFCP1 = cueTargetDelaySFCP1HF(isInDPulvinar,:);
xBounds = paramsHF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCP3, cueTargetDelaySFCP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC dPul (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-HF-dPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline LF
fAxis = fAxisLF;
cueTargetDelaySFCP3 = cueTargetDelaySFCP3LF(isInVPulvinar,:);
cueTargetDelaySFCP1 = cueTargetDelaySFCP1LF(isInVPulvinar,:);
xBounds = paramsLF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCP3, cueTargetDelaySFCP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC vPul (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-LF-vPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot SFC in cue-target delay relative to baseline LF
fAxis = fAxisHF;
cueTargetDelaySFCP3 = cueTargetDelaySFCP3HF(isInVPulvinar,:);
cueTargetDelaySFCP1 = cueTargetDelaySFCP1HF(isInVPulvinar,:);
xBounds = paramsHF.fpass;
yBounds = [0 0.1];

plotLfpPower(cueTargetDelaySFCP3, cueTargetDelaySFCP1, fAxis, xBounds, yBounds, p3Col, p1Col, 'P3', 'P1')

ylabel('Spike-Field Coherence');
title('Cue-Target Delay SFC vPul (No Ref)');

plotFileName = sprintf('%s/allSessions-cueTargetDelaySFC-HF-vPul-P3vsP1-v%d.png', outputDir, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');