function lfpAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds, ref)

% clear;
% readDataLocally;
% sessionInds = 8;%1:23;
% ref = 'CAR';

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
dPulSpikeVPulFieldCount = 0;
vPulSpikeDPulFieldCount = 0;

lfpNames = cell(nLfpsApprox, 1);
isInDPulvinar = false(nLfpsApprox, 1);
isInVPulvinar = false(nLfpsApprox, 1);

subPair = nan(nSessions, 2);
subPairCount = 0;

fAxisLF = NaN;
nFAxisLF = 17;
baselinePowerLF = nan(nLfpsApprox, nFAxisLF);
cueResponsePowerLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerAllLocsLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerP3LF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelayPowerP1LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelayPowerP3LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelayPowerP1LF = nan(nLfpsApprox, nFAxisLF);

baselineSFCLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP3LF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP1LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelaySFCP3LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelaySFCP1LF = nan(nLfpsApprox, nFAxisLF);

fAxisHF = NaN;
nFAxisHF = 77;
baselinePowerHF = nan(nLfpsApprox, nFAxisHF);
cueResponsePowerHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerAllLocsHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP1HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelayPowerP3HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelayPowerP1HF = nan(nLfpsApprox, nFAxisHF);

baselineSFCHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP1HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelaySFCP3HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelaySFCP1HF = nan(nLfpsApprox, nFAxisHF);

nSubPairsApprox = round(nSessions / 3);
baselineSubPairCohLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySubPairCohP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySubPairCohP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySubPairCohP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySubPairCohP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSubPairCohHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySubPairCohP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySubPairCohP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySubPairCohP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySubPairCohP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCSingleDPulSpikeVPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCSingleDPulSpikeVPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCSingleVPulSpikeDPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCSingleVPulSpikeDPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCDPulSpikeVPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCDPulSpikeVPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCVPulSpikeDPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCVPulSpikeDPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineWindowOffset = [-0.25 0];
cueResponseOffset = [0 0.25];
cueTargetDelayOffset = [-0.25 0];
targetDimDelayOffset = [-0.25 0];

% chronux parameters
paramsLF.tapers = [1 1];
paramsLF.fpass = [8 25];
paramsLF.pad = 2;
paramsLF.Fs = 1000;
paramsLF.trialave = 1;

paramsHF.tapers = [2 3];
paramsHF.fpass = [25 100];
paramsHF.pad = 2;
paramsHF.Fs = 1000;
paramsHF.trialave = 1;

adjCohNormRate = 5; % rate in Hz to normalize spike rates to -- too high leads to C > 1 and atanh(C) = undef
% TODO normalize to the other rate in comparison, not to general rate

fprintf('\n-------------------------------------------------------\n');
fprintf('Across Session Analysis\n');
fprintf('Reference: %s\n', ref);

%% session loop
for i = 1:nSessions
    sessionInd = sessionInds(i);
    % TODO
    % session 15, 'M20170329', 1-32
    % session 17, 'M20170329', 1-32
    % have abnormally high power around 40 Hz
    if sessionInd == 15 || sessionInd == 17
        fprintf('Excluding session %d...\n', sessionInd);
        continue;
    end
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
        lfpNames{lfpCount} = sprintf('%s_FP%03d', sessionName, EL.channelInds(j));
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
        targetDimDelayInd = getTimeLogicalWithTolerance(EL.targetDimBalLfp.t, targetDimDelayOffset);

        % compute power in baseline - should probably be bipolar reference
        % all locations together
        if strcmp(ref, 'RAW')
            cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:))';
            targetDimLfpP3Current = squeeze(EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 3,:))';
            targetDimLfpP1Current = squeeze(EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 1,:))';
        elseif strcmp(ref, 'BIP')
            cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j+1,:,:) - EL.cueOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j+1,:,:) - EL.arrayOnsetLfp.lfp(j,:,:))';
            arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j+1,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j+1,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:))';
            targetDimLfpP3Current = squeeze(EL.targetDimBalLfp.lfp(j+1,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 3,:))';
            targetDimLfpP1Current = squeeze(EL.targetDimBalLfp.lfp(j+1,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 1,:))';
        elseif strcmp(ref, 'CAR')
            caCh = numel(EL.channelInds) + 1;
            cueOnsetLfpCurrent = squeeze(EL.cueOnsetLfp.lfp(j,:,:) - EL.cueOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpCurrent = squeeze(EL.arrayOnsetLfp.lfp(j,:,:) - EL.arrayOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpP3Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(EL.arrayOnsetLfp.lfp(j,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 1,:))';
            targetDimLfpP3Current = squeeze(EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 3,:))';
            targetDimLfpP1Current = squeeze(EL.targetDimBalLfp.lfp(j,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 1,:))';
        end
        % ignore array shapes because we are looking at pre-array delay
        % period
        
        numTrials = size(cueOnsetLfpCurrent, 1);
        numTrialsP3 = size(arrayOnsetLfpP3Current, 1);
        numTrialsP1 = size(arrayOnsetLfpP1Current, 1);
        numTrialsBalP3 = size(targetDimLfpP3Current, 1);
        numTrialsBalP1 = size(targetDimLfpP1Current, 1);
        
        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);

        [baselinePowerLF(lfpCount,:),fAxisLF] = mtspectrumc(preCueBaselineLfps, paramsLF);
        [baselinePowerHF(lfpCount,:),fAxisHF] = mtspectrumc(preCueBaselineLfps, paramsHF);
        [cueTargetDelayPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP3, paramsLF);
        [cueTargetDelayPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP3, paramsHF);
        [cueTargetDelayPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP1, paramsLF);
        [cueTargetDelayPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP1, paramsHF);
        [targetDimDelayPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(targetDimDelayLfpsP3, paramsLF);
        [targetDimDelayPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(targetDimDelayLfpsP3, paramsHF);
        [targetDimDelayPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(targetDimDelayLfpsP1, paramsLF);
        [targetDimDelayPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(targetDimDelayLfpsP1, paramsHF);
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        baselineSFCLF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrials)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        baselineSFCHF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrials)-2)); % adjust for num trials
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        cueTargetDelaySFCP3LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        cueTargetDelaySFCP3HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        cueTargetDelaySFCP1LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        cueTargetDelaySFCP1HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.targetDimBalByLoc{3}, targetDimDelayOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        targetDimDelaySFCP3LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        targetDimDelaySFCP3HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.targetDimBalByLoc{1}, targetDimDelayOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        targetDimDelaySFCP1LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        targetDimDelaySFCP1HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
        
        % compute coherence across subdivisions, one pair per session
        % TODO use all pairs per session unless this double dips too much
        % due to within session correlation. should be fine in the bipolar
        % reference
        if any(R.dPulChannels) && any(R.vPulChannels)
            if EL.channelInds(j) == floor(mean(R.dPulChannels))
                subPair(i,1) = j;
            elseif EL.channelInds(j) == floor(mean(R.vPulChannels))
                subPair(i,2) = j;
                if all(~isnan(subPair(i,:))) % ventral always comes after dorsal so put this code in ventral
                    subPairCount = subPairCount + 1;
                    j1 = subPair(i,1);
                    j2 = subPair(i,2);
                    if strcmp(ref, 'RAW')
                        cueOnsetLfpCurrent1 = squeeze(EL.cueOnsetLfp.lfp(j1,:,:))';
                        arrayOnsetLfpCurrent1 = squeeze(EL.arrayOnsetLfp.lfp(j1,:,:))';
                        arrayOnsetLfpP3Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current1 = squeeze(EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current1 = squeeze(EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 1,:))';
                        cueOnsetLfpCurrent2 = squeeze(EL.cueOnsetLfp.lfp(j2,:,:))';
                        arrayOnsetLfpCurrent2 = squeeze(EL.arrayOnsetLfp.lfp(j2,:,:))';
                        arrayOnsetLfpP3Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current2 = squeeze(EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current2 = squeeze(EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 1,:))';
                    elseif strcmp(ref, 'BIP')
                        cueOnsetLfpCurrent1 = squeeze(EL.cueOnsetLfp.lfp(j1+1,:,:) - EL.cueOnsetLfp.lfp(j1,:,:))';
                        arrayOnsetLfpCurrent1 = squeeze(EL.arrayOnsetLfp.lfp(j1+1,:,:) - EL.arrayOnsetLfp.lfp(j1,:,:))';
                        arrayOnsetLfpP3Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1+1,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1+1,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current1 = squeeze(EL.targetDimBalLfp.lfp(j1+1,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current1 = squeeze(EL.targetDimBalLfp.lfp(j1+1,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 1,:))';
                        cueOnsetLfpCurrent2 = squeeze(EL.cueOnsetLfp.lfp(j2+1,:,:) - EL.cueOnsetLfp.lfp(j2,:,:))';
                        arrayOnsetLfpCurrent2 = squeeze(EL.arrayOnsetLfp.lfp(j2+1,:,:) - EL.arrayOnsetLfp.lfp(j2,:,:))';
                        arrayOnsetLfpP3Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2+1,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2+1,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current2 = squeeze(EL.targetDimBalLfp.lfp(j2+1,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current2 = squeeze(EL.targetDimBalLfp.lfp(j2+1,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 1,:))';
                    elseif strcmp(ref, 'CAR')
                        cueOnsetLfpCurrent1 = squeeze(EL.cueOnsetLfp.lfp(j1,:,:) - EL.cueOnsetLfp.lfp(caCh,:,:))';
                        arrayOnsetLfpCurrent1 = squeeze(EL.arrayOnsetLfp.lfp(j1,:,:) - EL.arrayOnsetLfp.lfp(caCh,:,:))';
                        arrayOnsetLfpP3Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current1 = squeeze(EL.arrayOnsetLfp.lfp(j1,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current1 = squeeze(EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current1 = squeeze(EL.targetDimBalLfp.lfp(j1,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 1,:))';
                        cueOnsetLfpCurrent2 = squeeze(EL.cueOnsetLfp.lfp(j2,:,:) - EL.cueOnsetLfp.lfp(caCh,:,:))';
                        arrayOnsetLfpCurrent2 = squeeze(EL.arrayOnsetLfp.lfp(j2,:,:) - EL.arrayOnsetLfp.lfp(caCh,:,:))';
                        arrayOnsetLfpP3Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 3,:))';
                        arrayOnsetLfpP1Current2 = squeeze(EL.arrayOnsetLfp.lfp(j2,EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 1,:))';
                        targetDimLfpP3Current2 = squeeze(EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 3,:))';
                        targetDimLfpP1Current2 = squeeze(EL.targetDimBalLfp.lfp(j2,EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 1,:))';
                    end
                    
                    numTrials = size(cueOnsetLfpCurrent1, 1);
                    numTrialsP3 = size(arrayOnsetLfpP3Current1, 1);
                    numTrialsP1 = size(arrayOnsetLfpP1Current1, 1);
                    numTrialsBalP3 = size(targetDimLfpP3Current1, 1);
                    numTrialsBalP1 = size(targetDimLfpP1Current1, 1);

                    preCueBaselineLfps1 = cueOnsetLfpCurrent1(baselineInd,:);
                    preCueBaselineLfps2 = cueOnsetLfpCurrent2(baselineInd,:);
                    cueTargetDelayLfps1P3 = arrayOnsetLfpP3Current1(cueTargetDelayInd,:); 
                    cueTargetDelayLfps2P3 = arrayOnsetLfpP3Current2(cueTargetDelayInd,:); 
                    cueTargetDelayLfps1P1 = arrayOnsetLfpP1Current1(cueTargetDelayInd,:);
                    cueTargetDelayLfps2P1 = arrayOnsetLfpP1Current2(cueTargetDelayInd,:);
                    targetDimDelayLfps1P3 = targetDimLfpP3Current1(targetDimDelayInd,:); 
                    targetDimDelayLfps2P3 = targetDimLfpP3Current2(targetDimDelayInd,:); 
                    targetDimDelayLfps1P1 = targetDimLfpP1Current1(targetDimDelayInd,:);
                    targetDimDelayLfps2P1 = targetDimLfpP1Current2(targetDimDelayInd,:);

                    baselineSubPairCohLF(subPairCount,:) = coherencyc(preCueBaselineLfps1, preCueBaselineLfps2, paramsLF);
                    baselineSubPairCohHF(subPairCount,:) = coherencyc(preCueBaselineLfps1, preCueBaselineLfps2, paramsHF);
                    cueTargetDelaySubPairCohP3LF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P3, cueTargetDelayLfps2P3, paramsLF);
                    cueTargetDelaySubPairCohP3HF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P3, cueTargetDelayLfps2P3, paramsHF);
                    cueTargetDelaySubPairCohP1LF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P1, cueTargetDelayLfps2P1, paramsLF);
                    cueTargetDelaySubPairCohP1HF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P1, cueTargetDelayLfps2P1, paramsHF);
                    targetDimDelaySubPairCohP3LF(subPairCount,:) = coherencyc(targetDimDelayLfps1P3, targetDimDelayLfps2P3, paramsLF);
                    targetDimDelaySubPairCohP3HF(subPairCount,:) = coherencyc(targetDimDelayLfps1P3, targetDimDelayLfps2P3, paramsHF);
                    targetDimDelaySubPairCohP1LF(subPairCount,:) = coherencyc(targetDimDelayLfps1P1, targetDimDelayLfps2P1, paramsLF);
                    targetDimDelaySubPairCohP1HF(subPairCount,:) = coherencyc(targetDimDelayLfps1P1, targetDimDelayLfps2P1, paramsHF);
                    
                    muaChannel1Ts = EL.allMUAStructs{j1}.ts;
                    muaChannel2Ts = EL.allMUAStructs{j2}.ts;
                    
                    % dPul spike, vPul field - one each per session
                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps2, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleDPulSpikeVPulFieldLF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrials)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps2, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleDPulSpikeVPulFieldHF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrials)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfps2P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleDPulSpikeVPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfps2P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleDPulSpikeVPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfps2P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleDPulSpikeVPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfps2P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleDPulSpikeVPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.targetDimBalByLoc{3}, targetDimDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfps2P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleDPulSpikeVPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfps2P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleDPulSpikeVPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.targetDimBalByLoc{1}, targetDimDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfps2P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleDPulSpikeVPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfps2P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleDPulSpikeVPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
                    
                    % vPul spike, dPul field - one each per session
                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleVPulSpikeDPulFieldLF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrials)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleVPulSpikeDPulFieldHF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrials)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfps1P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleVPulSpikeDPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfps1P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleVPulSpikeDPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfps1P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleVPulSpikeDPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfps1P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    cueTargetDelaySFCSingleVPulSpikeDPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.targetDimBalByLoc{3}, targetDimDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfps1P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleVPulSpikeDPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfps1P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleVPulSpikeDPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.targetDimBalByLoc{1}, targetDimDelayOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfps1P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleVPulSpikeDPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfps1P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    targetDimDelaySFCSingleVPulSpikeDPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials                    
                end
            end                
        end
    end
    
    % across subdivision spike field coherence
    if any(R.dPulChannels) && any(R.vPulChannels)
        
        % dPul spikes separately to mean vPul LFP
        fprintf('Computing dorsal pulvinar spikes - ventral pulvinar LFPs coherence...\n');
        js = arrayfun(@(x) find(x == EL.channelInds), R.vPulChannels);
        if strcmp(ref, 'RAW')
            cueOnsetLfpCurrent = squeeze(mean(EL.cueOnsetLfp.lfp(js,:,:), 1))';
            arrayOnsetLfpCurrent = squeeze(mean(EL.arrayOnsetLfp.lfp(js,:,:), 1))';
            arrayOnsetLfpP3Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 3,:), 1))';
            arrayOnsetLfpP1Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 1,:), 1))';
            targetDimLfpP3Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 3,:), 1))';
            targetDimLfpP1Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 1,:), 1))';
        elseif strcmp(ref, 'BIP')
            if numel(js) == 1
                continue;
            end
            % mean of bipolar signals: ((a_2-a_1) + (a_3-a_2) + (a_n-a_{n-1})) + ...)/n 
            % is equal to (a_n-a_1)/n
            cueOnsetLfpCurrent = squeeze((EL.cueOnsetLfp.lfp(js(end),:,:) - EL.cueOnsetLfp.lfp(js(1),:,:))/numel(js))';
            arrayOnsetLfpCurrent = squeeze((EL.arrayOnsetLfp.lfp(js(end),:,:) - EL.arrayOnsetLfp.lfp(js(1),:,:))/numel(js))';
            arrayOnsetLfpP3Current = squeeze((EL.arrayOnsetLfp.lfp(js(end),EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(js(1),EL.UE.cueLoc == 3,:))/numel(js))';
            arrayOnsetLfpP1Current = squeeze((EL.arrayOnsetLfp.lfp(js(end),EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(js(1),EL.UE.cueLoc == 1,:))/numel(js))';
            targetDimLfpP3Current = squeeze((EL.targetDimBalLfp.lfp(js(end),EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(js(1),EL.UE.cueLocHoldBal == 3,:))/numel(js))';
            targetDimLfpP1Current = squeeze((EL.targetDimBalLfp.lfp(js(end),EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(js(1),EL.UE.cueLocHoldBal == 1,:))/numel(js))';
        elseif strcmp(ref, 'CAR')
            caCh = numel(EL.channelInds) + 1;
            cueOnsetLfpCurrent = squeeze(mean(EL.cueOnsetLfp.lfp(js,:,:), 1) - EL.cueOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpCurrent = squeeze(mean(EL.arrayOnsetLfp.lfp(js,:,:), 1) - EL.arrayOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpP3Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 3,:), 1) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 1,:), 1) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 1,:))';
            targetDimLfpP3Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 3,:), 1) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 3,:))';
            targetDimLfpP1Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 1,:), 1) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 1,:))';
        end
        
        numTrials = size(cueOnsetLfpCurrent, 1);
        numTrialsP3 = size(arrayOnsetLfpP3Current, 1);
        numTrialsP1 = size(arrayOnsetLfpP1Current, 1);
        numTrialsBalP3 = size(targetDimLfpP3Current, 1);
        numTrialsBalP1 = size(targetDimLfpP1Current, 1);

        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);
        
        for k = 1:numel(R.dPulChannels)
            j = find(EL.channelInds == R.dPulChannels(k));
            fprintf('Processing channel %d...\n', EL.channelInds(j));
            spikeTs = EL.allMUAStructs{j}.ts;
            
            dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount + 1;
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCDPulSpikeVPulFieldLF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrials)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCDPulSpikeVPulFieldHF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrials)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCDPulSpikeVPulFieldP3LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCDPulSpikeVPulFieldP3HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCDPulSpikeVPulFieldP1LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCDPulSpikeVPulFieldP1HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{3}, targetDimDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCDPulSpikeVPulFieldP3LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCDPulSpikeVPulFieldP3HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{1}, targetDimDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCDPulSpikeVPulFieldP1LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCDPulSpikeVPulFieldP1HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
        end
        
        % vPul spikes all together to dPul LFP separately
        fprintf('Computing ventral pulvinar spikes - dorsal pulvinar LFPs coherence...\n');
        js = arrayfun(@(x) find(x == EL.channelInds), R.dPulChannels);
        if strcmp(ref, 'RAW')
            cueOnsetLfpCurrent = squeeze(mean(EL.cueOnsetLfp.lfp(js,:,:), 1))';
            arrayOnsetLfpCurrent = squeeze(mean(EL.arrayOnsetLfp.lfp(js,:,:), 1))';
            arrayOnsetLfpP3Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 3,:), 1))';
            arrayOnsetLfpP1Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 1,:), 1))';
            targetDimLfpP3Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 3,:), 1))';
            targetDimLfpP1Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 1,:), 1))';
        elseif strcmp(ref, 'BIP')
            if numel(js) == 1
                continue;
            end
            % mean of bipolar signals: ((a_2-a_1) + (a_3-a_2) + (a_n-a_{n-1})) + ...)/n 
            % is equal to (a_n-a_1)/n
            cueOnsetLfpCurrent = squeeze((EL.cueOnsetLfp.lfp(js(end),:,:) - EL.cueOnsetLfp.lfp(js(1),:,:))/numel(js))';
            arrayOnsetLfpCurrent = squeeze((EL.arrayOnsetLfp.lfp(js(end),:,:) - EL.arrayOnsetLfp.lfp(js(1),:,:))/numel(js))';
            arrayOnsetLfpP3Current = squeeze((EL.arrayOnsetLfp.lfp(js(end),EL.UE.cueLoc == 3,:) - EL.arrayOnsetLfp.lfp(js(1),EL.UE.cueLoc == 3,:))/numel(js))';
            arrayOnsetLfpP1Current = squeeze((EL.arrayOnsetLfp.lfp(js(end),EL.UE.cueLoc == 1,:) - EL.arrayOnsetLfp.lfp(js(1),EL.UE.cueLoc == 1,:))/numel(js))';
            targetDimLfpP3Current = squeeze((EL.targetDimBalLfp.lfp(js(end),EL.UE.cueLocHoldBal == 3,:) - EL.targetDimBalLfp.lfp(js(1),EL.UE.cueLocHoldBal == 3,:))/numel(js))';
            targetDimLfpP1Current = squeeze((EL.targetDimBalLfp.lfp(js(end),EL.UE.cueLocHoldBal == 1,:) - EL.targetDimBalLfp.lfp(js(1),EL.UE.cueLocHoldBal == 1,:))/numel(js))';
        elseif strcmp(ref, 'CAR')
            caCh = numel(EL.channelInds) + 1;
            cueOnsetLfpCurrent = squeeze(mean(EL.cueOnsetLfp.lfp(js,:,:), 1) - EL.cueOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpCurrent = squeeze(mean(EL.arrayOnsetLfp.lfp(js,:,:), 1) - EL.arrayOnsetLfp.lfp(caCh,:,:))';
            arrayOnsetLfpP3Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 3,:), 1) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 3,:))';
            arrayOnsetLfpP1Current = squeeze(mean(EL.arrayOnsetLfp.lfp(js,EL.UE.cueLoc == 1,:), 1) - EL.arrayOnsetLfp.lfp(caCh,EL.UE.cueLoc == 1,:))';
            targetDimLfpP3Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 3,:), 1) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 3,:))';
            targetDimLfpP1Current = squeeze(mean(EL.targetDimBalLfp.lfp(js,EL.UE.cueLocHoldBal == 1,:), 1) - EL.targetDimBalLfp.lfp(caCh,EL.UE.cueLocHoldBal == 1,:))';
        end
        
        numTrials = size(cueOnsetLfpCurrent, 1);
        numTrialsP3 = size(arrayOnsetLfpP3Current, 1);
        numTrialsP1 = size(arrayOnsetLfpP1Current, 1);
        numTrialsBalP3 = size(targetDimLfpP3Current, 1);
        numTrialsBalP1 = size(targetDimLfpP1Current, 1);

        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);
        
        for k = 1:numel(R.vPulChannels)
            j = find(EL.channelInds == R.vPulChannels(k));
            fprintf('Processing channel %d...\n', EL.channelInds(j));
            spikeTs = EL.allMUAStructs{j}.ts;
            
            vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount + 1;
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCVPulSpikeDPulFieldLF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrials)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCVPulSpikeDPulFieldHF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrials)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{3}, cueTargetDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCVPulSpikeDPulFieldP3LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCVPulSpikeDPulFieldP3HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{1}, cueTargetDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, cueTargetDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCVPulSpikeDPulFieldP1LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(cueTargetDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            cueTargetDelaySFCVPulSpikeDPulFieldP1HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{3}, targetDimDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCVPulSpikeDPulFieldP3LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCVPulSpikeDPulFieldP3HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.targetDimBalByLoc{1}, targetDimDelayOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, targetDimDelayOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCVPulSpikeDPulFieldP1LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(targetDimDelayLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            targetDimDelaySFCVPulSpikeDPulFieldP1HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsBalP1)-2)); % adjust for num trials
        end
    end

    
    clear EL;
end


%%
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-%s-v%d.mat', outputDir, ref, v);
save(saveFileName);

%%
outputDirOrig = outputDir;
saveFileName = sprintf('%s/allSessions-lfpAnalysisSummary-%s-v%d.mat', outputDir, ref, v);
load(saveFileName);
outputDir = outputDirOrig;

%%
cols = lines(6);
dPulCol = cols(3,:);
vPulCol = cols(5,:);

p3Col = [0.9 0 0];
p1Col = [0 0 0.9];

if strcmp(ref, 'BIP')
    relPowYBounds = [-0.35 0.1];
elseif strcmp(ref, 'CAR')
    relPowYBounds = [-2 0.25];
end

sfcBaselineYBounds = [0 0.06];
sfcCTDelayYBounds = [0 0.08];
sfcTDDelayYBounds = [0 0.12];

%% plot power in pre-cue baseline dPul vs vPul
plotLfpPower2(baselinePowerLF, baselinePowerHF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [-35 -17; -35 -17], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Power (dB/Hz)', ...
        'titleText', sprintf('Pre-Cue Baseline Power (%s)', ref), ...
        'doDB', 1);
    
plotFileName = sprintf('%s/allSessions-baselinePowerDB-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay dPul vs vPul
plotLfpPower2(cueTargetDelayPowerP3LF, cueTargetDelayPowerP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [-35 -17; -35 -17], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Power (dB/Hz)', ...
        'titleText', sprintf('Cue-Target Delay Power (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-ctDelayPowerDB-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay relative to baseline dPul vs vPul
cueTargetDelayRelPowerP3LF = cueTargetDelayPowerP3LF ./ baselinePowerLF;
cueTargetDelayRelPowerP3HF = cueTargetDelayPowerP3HF ./ baselinePowerHF;

plotLfpPower2(cueTargetDelayRelPowerP3LF, cueTargetDelayRelPowerP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Cue-Target Delay Relative Power (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay dPul P3 vs P1
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelPowerLF = [cueTargetDelayPowerP3LF(channelCond,:); cueTargetDelayPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
cueTargetDelayRelPowerHF = [cueTargetDelayPowerP3HF(channelCond,:); cueTargetDelayPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(cueTargetDelayRelPowerLF, cueTargetDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Cue-Target Delay Relative Power - Dorsal Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in cue-target delay vPul P3 vs P1
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelPowerLF = [cueTargetDelayPowerP3LF(channelCond,:); cueTargetDelayPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
cueTargetDelayRelPowerHF = [cueTargetDelayPowerP3HF(channelCond,:); cueTargetDelayPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(cueTargetDelayRelPowerLF, cueTargetDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Cue-Target Delay Relative Power - Ventral Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-ctDelayRelPowerDB-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in target-dim delay dPul vs vPul
plotLfpPower2(targetDimDelayPowerP3LF, targetDimDelayPowerP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [-35 -17; -35 -17], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Power (dB/Hz)', ...
        'titleText', sprintf('Target-Dim Delay Power (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-tdDelayPowerDB-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in target-dim delay relative to baseline dPul vs vPul
targetDimDelayRelPowerP3LF = targetDimDelayPowerP3LF ./ baselinePowerLF;
targetDimDelayRelPowerP3HF = targetDimDelayPowerP3HF ./ baselinePowerHF;

plotLfpPower2(targetDimDelayRelPowerP3LF, targetDimDelayRelPowerP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Target-Dim Delay Relative Power (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-tdDelayRelPowerDB-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in target-dim delay dPul P3 vs P1
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
targetDimDelayRelPowerLF = [targetDimDelayPowerP3LF(channelCond,:); targetDimDelayPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
targetDimDelayRelPowerHF = [targetDimDelayPowerP3HF(channelCond,:); targetDimDelayPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(targetDimDelayRelPowerLF, targetDimDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Target-Dim Delay Relative Power - Dorsal Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-tdDelayRelPowerDB-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in target-dim delay vPul P3 vs P1
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
targetDimDelayRelPowerLF = [targetDimDelayPowerP3LF(channelCond,:); targetDimDelayPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
targetDimDelayRelPowerHF = [targetDimDelayPowerP3HF(channelCond,:); targetDimDelayPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(targetDimDelayRelPowerLF, targetDimDelayRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Target-Dim Delay Relative Power - Ventral Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-tdDelayRelPowerDB-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%% plot baseline SFC
plotLfpPower2(baselineSFCLF, baselineSFCHF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcBaselineYBounds; sfcBaselineYBounds], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Pre-Cue Baseline Local SFC (%s)', ref), ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-baselineSFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC
plotLfpPower2(cueTargetDelaySFCP3LF, cueTargetDelaySFCP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Cue-Target Delay Local SFC (%s)', ref), ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-ctDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot cue-target delay SFC dPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
cueTargetDelayRelSFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelayRelSFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Cue-Target Delay SFC - Dorsal Pulvinar (%s)', ref), ...
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

plotLfpPower2(cueTargetDelayRelSFCLF, cueTargetDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Cue-Target Delay SFC - Ventral Pulvinar (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot target-dim delay SFC
plotLfpPower2(targetDimDelaySFCP3LF, targetDimDelaySFCP3HF, fAxisLF, fAxisHF, [isInDPulvinar isInVPulvinar], ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcTDDelayYBounds; sfcTDDelayYBounds], ...
        'cols', [dPulCol; vPulCol], ...
        'lineLabels', {'dPul', 'vPul'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Target-Dim Delay Local SFC (%s)', ref), ...
        'doDB', 0);
    
plotFileName = sprintf('%s/allSessions-tdDelaySFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot target-dim delay SFC dPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
targetDimDelayRelSFCLF = [targetDimDelaySFCP3LF(channelCond,:); targetDimDelaySFCP1LF(channelCond,:)];
targetDimDelayRelSFCHF = [targetDimDelaySFCP3HF(channelCond,:); targetDimDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(targetDimDelayRelSFCLF, targetDimDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcTDDelayYBounds; sfcTDDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Target-Dim Delay SFC - Dorsal Pulvinar (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-tdDelaySFC-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot target-dim delay SFC vPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
targetDimDelayRelSFCLF = [targetDimDelaySFCP3LF(channelCond,:); targetDimDelaySFCP1LF(channelCond,:)];
targetDimDelayRelSFCHF = [targetDimDelaySFCP3HF(channelCond,:); targetDimDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(targetDimDelayRelSFCLF, targetDimDelayRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcTDDelayYBounds; sfcTDDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Target-Dim Delay SFC - Ventral Pulvinar (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-tdDelaySFC-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision pair coherence across conditions and periods
% using CAR is tricky here because CAR is a common signal to the dorsal and
% ventral channels and any fluctuations in CAR will be picked up here
nSubPairs = size(baselineSubPairCohLF, 1);
cohAllLF = [baselineSubPairCohLF; cueTargetDelaySubPairCohP3LF; cueTargetDelaySubPairCohP1LF; targetDimDelaySubPairCohP3LF; targetDimDelaySubPairCohP1LF];
cohAllHF = [baselineSubPairCohHF; cueTargetDelaySubPairCohP3HF; cueTargetDelaySubPairCohP1HF; targetDimDelaySubPairCohP3HF; targetDimDelaySubPairCohP1HF];
condLogical = false(nSubPairs * 5, 5);
for i = 1:5
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

if strcmp(ref, 'BIP')
    yBounds = [0 0.15];
elseif strcmp(ref, 'CAR')
    yBounds = [0 0.8]; % !!!
end
baselineCol = [0 0.9 0];
plotLfpPower2(cohAllLF, cohAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [baselineCol; p3Col; p1Col; p3Col * 0.5; p1Col * 0.5], ...
        'lineLabels', {'Baseline', 'CT Delay P3', 'CT Delay P1', 'TD Delay P3', 'TD Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Across Subdivision Coherence (N=%d; %s)', nSubPairs, ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-subPairCoh-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision pair dPul spike - vPul field coherence across conditions and periods
nSubPairs = size(baselineSFCSingleDPulSpikeVPulFieldLF, 1);
sfcAllLF = [baselineSFCSingleDPulSpikeVPulFieldLF; cueTargetDelaySFCSingleDPulSpikeVPulFieldP3LF; cueTargetDelaySFCSingleDPulSpikeVPulFieldP1LF; targetDimDelaySFCSingleDPulSpikeVPulFieldP3LF; targetDimDelaySFCSingleDPulSpikeVPulFieldP1LF];
sfcAllHF = [baselineSFCSingleDPulSpikeVPulFieldHF; cueTargetDelaySFCSingleDPulSpikeVPulFieldP3HF; cueTargetDelaySFCSingleDPulSpikeVPulFieldP1HF; targetDimDelaySFCSingleDPulSpikeVPulFieldP3HF; targetDimDelaySFCSingleDPulSpikeVPulFieldP1HF];
condLogical = false(nSubPairs * 5, 5);
for i = 1:5
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.13];
baselineCol = [0 0.9 0];
plotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [baselineCol; p3Col; p1Col; p3Col * 0.5; p1Col * 0.5], ...
        'lineLabels', {'Baseline', 'CT Delay P3', 'CT Delay P1', 'TD Delay P3', 'TD Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('dPul Spike - vPul Field Coherence (N=%d; %s)', nSubPairs, ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-SFCSingle-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision pair vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCSingleVPulSpikeDPulFieldLF, 1);
sfcAllLF = [baselineSFCSingleVPulSpikeDPulFieldLF; cueTargetDelaySFCSingleVPulSpikeDPulFieldP3LF; cueTargetDelaySFCSingleVPulSpikeDPulFieldP1LF; targetDimDelaySFCSingleVPulSpikeDPulFieldP3LF; targetDimDelaySFCSingleVPulSpikeDPulFieldP1LF];
sfcAllHF = [baselineSFCSingleVPulSpikeDPulFieldHF; cueTargetDelaySFCSingleVPulSpikeDPulFieldP3HF; cueTargetDelaySFCSingleVPulSpikeDPulFieldP1HF; targetDimDelaySFCSingleVPulSpikeDPulFieldP3HF; targetDimDelaySFCSingleVPulSpikeDPulFieldP1HF];
condLogical = false(nSubPairs * 5, 5);
for i = 1:5
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.13];
baselineCol = [0 0.9 0];
plotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [baselineCol; p3Col; p1Col; p3Col * 0.5; p1Col * 0.5], ...
        'lineLabels', {'Baseline', 'CT Delay P3', 'CT Delay P1', 'TD Delay P3', 'TD Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('vPul Spike - dPul Field Coherence (N=%d; %s)', nSubPairs, ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-SFCSingle-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs dPul spike - vPul field coherence across conditions and periods
nSubPairs = size(baselineSFCDPulSpikeVPulFieldLF, 1);
sfcAllLF = [baselineSFCDPulSpikeVPulFieldLF; cueTargetDelaySFCDPulSpikeVPulFieldP3LF; cueTargetDelaySFCDPulSpikeVPulFieldP1LF; targetDimDelaySFCDPulSpikeVPulFieldP3LF; targetDimDelaySFCDPulSpikeVPulFieldP1LF];
sfcAllHF = [baselineSFCDPulSpikeVPulFieldHF; cueTargetDelaySFCDPulSpikeVPulFieldP3HF; cueTargetDelaySFCDPulSpikeVPulFieldP1HF; targetDimDelaySFCDPulSpikeVPulFieldP3HF; targetDimDelaySFCDPulSpikeVPulFieldP1HF];
condLogical = false(nSubPairs * 5, 5);
for i = 1:5
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.06];
baselineCol = [0 0.9 0];
plotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [baselineCol; p3Col; p1Col; p3Col * 0.5; p1Col * 0.5], ...
        'lineLabels', {'Baseline', 'CT Delay P3', 'CT Delay P1', 'TD Delay P3', 'TD Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('dPul Spike - vPul Field Coherence (N=%d; %s)', nSubPairs, ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-dPulSpike-vPulField-SFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot subdivision all pairs vPul spike - dPul field coherence across conditions and periods
nSubPairs = size(baselineSFCVPulSpikeDPulFieldLF, 1);
sfcAllLF = [baselineSFCVPulSpikeDPulFieldLF; cueTargetDelaySFCVPulSpikeDPulFieldP3LF; cueTargetDelaySFCVPulSpikeDPulFieldP1LF; targetDimDelaySFCVPulSpikeDPulFieldP3LF; targetDimDelaySFCVPulSpikeDPulFieldP1LF];
sfcAllHF = [baselineSFCVPulSpikeDPulFieldHF; cueTargetDelaySFCVPulSpikeDPulFieldP3HF; cueTargetDelaySFCVPulSpikeDPulFieldP1HF; targetDimDelaySFCVPulSpikeDPulFieldP3HF; targetDimDelaySFCVPulSpikeDPulFieldP1HF];
condLogical = false(nSubPairs * 5, 5);
for i = 1:5
    condLogical(((i-1)*nSubPairs+1):(i*nSubPairs), i) = 1;
end

yBounds = [0 0.06];
baselineCol = [0 0.9 0];
plotLfpPower2(sfcAllLF, sfcAllHF, fAxisLF, fAxisHF, condLogical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [yBounds; yBounds], ...
        'cols', [baselineCol; p3Col; p1Col; p3Col * 0.5; p1Col * 0.5], ...
        'lineLabels', {'Baseline', 'CT Delay P3', 'CT Delay P1', 'TD Delay P3', 'TD Delay P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('vPul Spike - dPul Field Coherence (N=%d; %s)', nSubPairs, ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-vPulSpike-dPulField-SFC-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');


%%
return; % end plots, rest is live testing

%% inspect cue-target power in dPul and vPul at 8-12 Hz
fBounds = [8 12];
fInd = fAxisLF >= fBounds(1) & fAxisLF <= fBounds(2); 

powP3 = mean(cueTargetDelayPowerP3LF(isInDPulvinar | isInVPulvinar, fInd), 2); % mean over f
powP1 = mean(cueTargetDelayPowerP1LF(isInDPulvinar | isInVPulvinar, fInd), 2); % mean over f
powDiff = powP3 - powP1;

histBins = 0:0.002:0.05;

figure_tr_inch(5, 10);
subaxis(2, 1, 1);
histogram(powP3, histBins);
subaxis(2, 1, 2);
histogram(powP1, histBins);

histDiffBins = -0.02:0.001:0.02;
figure_tr_inch(6, 6);
histogram(powDiff, histDiffBins);

figure;
hold on;
plot([0 0.1], [0 0.1], 'Color', 0.3*ones(3, 1));
plot(powP3, powP1, '.', 'Color', lines(1), 'MarkerSize', 10);
xlim([0 0.03]);
ylim([0 0.03]);

[h,p,stats] = ttest(powDiff);
fprintf('Median diff power = %0.3f (N = %d), p = %0.3f\n', median(powDiff), size(powDiff, 1), p);

%% inspect cue-target power in dPul at 75-85 Hz
fBounds = [75 85];
fInd = fAxisHF >= fBounds(1) & fAxisHF <= fBounds(2); 

powP3 = mean(cueTargetDelayPowerP3HF(isInDPulvinar, fInd), 2); % mean over f
powP1 = mean(cueTargetDelayPowerP1HF(isInDPulvinar, fInd), 2); % mean over f
powDiff = powP3 - powP1;

histBins = 0:0.0002:0.005;

figure_tr_inch(5, 10);
subaxis(2, 1, 1);
histogram(powP3, histBins);
subaxis(2, 1, 2);
histogram(powP1, histBins);

histDiffBins = -0.002:0.0001:0.002;
figure_tr_inch(6, 6);
histogram(powDiff, histDiffBins);

figure;
hold on;
plot([0 0.1], [0 0.1], 'Color', 0.3*ones(3, 1));
plot(powP3, powP1, '.', 'Color', lines(1), 'MarkerSize', 10);
xlim([0 0.003]);
ylim([0 0.003]);

[h,p,stats] = ttest(powDiff);
fprintf('Median diff power = %0.3f (N = %d), p = %0.3f\n', median(powDiff), size(powDiff, 1), p);

%% inspect cue-target SFC in dPul at 75-85 Hz
fBounds = [75 85];
fInd = fAxisHF >= fBounds(1) & fAxisHF <= fBounds(2); 

sfcP3 = mean(cueTargetDelaySFCP3HF(isInDPulvinar, fInd), 2); % mean over f
sfcP1 = mean(cueTargetDelaySFCP1HF(isInDPulvinar, fInd), 2); % mean over f
sfcDiff = sfcP3 - sfcP1;

histBins = 0:0.01:0.2;

figure_tr_inch(5, 10);
subaxis(2, 1, 1);
histogram(sfcP3, histBins);
subaxis(2, 1, 2);
histogram(sfcP1, histBins);

histDiffBins = -0.15:0.01:0.15;
figure_tr_inch(6, 6);
histogram(sfcDiff, histDiffBins);

figure;
hold on;
plot([0 0.1], [0 0.1], 'Color', 0.3*ones(3, 1));
plot(sfcP3, sfcP1, '.', 'Color', lines(1), 'MarkerSize', 10);
xlim([0 0.1]);
ylim([0 0.1]);

[h,p,stats] = ttest(sfcDiff);
fprintf('Median diff SFC = %0.3f (N = %d), p = %0.3f\n', median(sfcDiff), size(sfcDiff, 1), p);

%% inspect cue-target SFC in vPul at 18-22 Hz
fBounds = [18 22];
fInd = fAxisLF >= fBounds(1) & fAxisLF <= fBounds(2); 

sfcP3 = mean(cueTargetDelaySFCP3LF(isInVPulvinar, fInd), 2); % mean over f
sfcP1 = mean(cueTargetDelaySFCP1LF(isInVPulvinar, fInd), 2); % mean over f
sfcDiff = sfcP3 - sfcP1;

histBins = 0:0.01:0.2;

figure_tr_inch(5, 10);
subaxis(2, 1, 1);
histogram(sfcP3, histBins);
subaxis(2, 1, 2);
histogram(sfcP1, histBins);

histDiffBins = -0.15:0.01:0.15;
figure_tr_inch(6, 6);
histogram(sfcDiff, histDiffBins);

figure;
hold on;
plot([0 0.1], [0 0.1], 'Color', 0.3*ones(3, 1));
plot(sfcP3, sfcP1, '.', 'Color', lines(1), 'MarkerSize', 10);
xlim([0 0.1]);
ylim([0 0.1]);

[h,p,stats] = ttest(sfcDiff);
fprintf('Median diff SFC = %0.3f (N = %d), p = %0.3f\n', median(sfcDiff), size(sfcDiff, 1), p);
