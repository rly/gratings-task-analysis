function lfpAnalysisSummary(processedDataRootDir, recordingInfoFileName, sessionInds, ref)

% clear;
% readDataLocally;
% sessionInds = 1:23;
% ref = 'CAR';

v = 13;

%%
recordingInfo = readRecordingInfo(recordingInfoFileName);

outputDir = sprintf('%s/%s/', processedDataRootDir, 'LFP_GRATINGS_SUMMARY');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
muaOutputDir = sprintf('%s/%s/', processedDataRootDir, 'SUA_MUA_GRATINGS_SUMMARY');
muaUnitNamesPulFileName = sprintf('%s/unitNamesPul-v%d.mat', muaOutputDir, v);
load(muaUnitNamesPulFileName);

if isempty(sessionInds)
    sessionInds = 1:numel(recordingInfo);
end

nSessions = numel(sessionInds);
nLfpsApprox = nSessions * 8; % should be equal or an underestimate
lfpCount = 0;
dPulSpikeVPulFieldCount = 0;
vPulSpikeDPulFieldCount = 0;

lfpNames = cell(nLfpsApprox, 1);
sfcNames = cell(nLfpsApprox, 1);
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
arrayResponseHoldPowerP3LF = nan(nLfpsApprox, nFAxisLF);
arrayResponseHoldPowerP1LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelayPowerP3LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelayPowerP1LF = nan(nLfpsApprox, nFAxisLF);

baselineSFCLF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP3LF = nan(nLfpsApprox, nFAxisLF);
cueTargetDelaySFCP1LF = nan(nLfpsApprox, nFAxisLF);
arrayResponseHoldSFCP3LF = nan(nLfpsApprox, nFAxisLF);
arrayResponseHoldSFCP1LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelaySFCP3LF = nan(nLfpsApprox, nFAxisLF);
targetDimDelaySFCP1LF = nan(nLfpsApprox, nFAxisLF);

fAxisHF = NaN;
nFAxisHF = 77;
baselinePowerHF = nan(nLfpsApprox, nFAxisHF);
cueResponsePowerHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerAllLocsHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelayPowerP1HF = nan(nLfpsApprox, nFAxisHF);
arrayResponseHoldPowerP3HF = nan(nLfpsApprox, nFAxisHF);
arrayResponseHoldPowerP1HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelayPowerP3HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelayPowerP1HF = nan(nLfpsApprox, nFAxisHF);

baselineSFCHF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP3HF = nan(nLfpsApprox, nFAxisHF);
cueTargetDelaySFCP1HF = nan(nLfpsApprox, nFAxisHF);
arrayResponseHoldSFCP3HF = nan(nLfpsApprox, nFAxisHF);
arrayResponseHoldSFCP1HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelaySFCP3HF = nan(nLfpsApprox, nFAxisHF);
targetDimDelaySFCP1HF = nan(nLfpsApprox, nFAxisHF);

nSubPairsApprox = round(nSessions / 3);
baselineSubPairCohLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySubPairCohP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySubPairCohP1LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSubPairCohP3LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSubPairCohP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySubPairCohP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySubPairCohP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSubPairCohHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySubPairCohP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySubPairCohP1HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSubPairCohP3HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSubPairCohP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySubPairCohP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySubPairCohP1HF = nan(nSubPairsApprox, nFAxisHF);

dPulSpikeVPulFieldNames = cell(nSubPairsApprox, 1);
vPulSpikeDPulFieldNames = cell(nSubPairsApprox, 1);

baselineSFCSingleDPulSpikeVPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCSingleDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCSingleDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCSingleDPulSpikeVPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCSingleDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCSingleDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCSingleVPulSpikeDPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCSingleVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCSingleVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCSingleVPulSpikeDPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCSingleVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCSingleVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCSingleVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCSingleVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCDPulSpikeVPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCDPulSpikeVPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCDPulSpikeVPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCDPulSpikeVPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCDPulSpikeVPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCDPulSpikeVPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineSFCVPulSpikeDPulFieldLF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
cueTargetDelaySFCVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
arrayResponseHoldSFCVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCVPulSpikeDPulFieldP3LF = nan(nSubPairsApprox, nFAxisLF);
targetDimDelaySFCVPulSpikeDPulFieldP1LF = nan(nSubPairsApprox, nFAxisLF);

baselineSFCVPulSpikeDPulFieldHF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
cueTargetDelaySFCVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
arrayResponseHoldSFCVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCVPulSpikeDPulFieldP3HF = nan(nSubPairsApprox, nFAxisHF);
targetDimDelaySFCVPulSpikeDPulFieldP1HF = nan(nSubPairsApprox, nFAxisHF);

baselineWindowOffset = [-0.25 0];
cueResponseOffset = [0 0.25];
cueTargetDelayOffset = [-0.4 0];
arrayResponseOffset = [0.2 0.6];
targetDimDelayOffset = [-0.4 0];

% chronux parameters
paramsLF.tapers = [1 1];
paramsLF.fpass = [8 25];
paramsLF.pad = 1;
paramsLF.Fs = 1000;
paramsLF.trialave = 1;
paramsBaselineLF = paramsLF;
paramsBaselineLF.pad = 2; % less padding for smaller window so that faxis the same

paramsHF.tapers = [4 7];
paramsHF.fpass = [25 100];
paramsHF.pad = 1;
paramsHF.Fs = 1000;
paramsHF.trialave = 1;
paramsBaselineHF = paramsHF;
paramsBaselineHF.pad = 2; % less padding for smaller window so that faxis the same

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
        if ~ismember(EL.channelInds(j), R.dPulChannels) && ~ismember(EL.channelInds(j), R.vPulChannels)
            fprintf('Skipping channel %d - not in pulvinar...\n', EL.channelInds(j));
            continue; % only process channels in pulvinar for now
        end
        
        % read MUA from another channel in the same subdivision
        if ismember(EL.channelInds(j), R.dPulChannels)
            isFoundMUA = 0;
            offset = -1; % -1, 1, -2, 2, -3, 3, ...
            while ~isFoundMUA
                if j + offset > 0 && j + offset <= numel(EL.channelInds)
                    tempUnitName = sprintf('%s_PUL_%dM', sessionName, EL.channelInds(j + offset));
                    if any(strcmp(unitNamesDPul, tempUnitName)) % units that have significant cue response based on saved list
                        isFoundMUA = 1;
                        muaNearbyChannelTs = EL.allMUAStructs{j + offset}.ts;
                    end
                end
                if offset < 0
                    offset = -1 * offset;
                else
                    offset = -1 * offset - 1;
                end
                if abs(offset) == 32 % no possible units found
                    fprintf('No possible MUA found to match channel %d...\n', EL.channelInds(j));
                    break;
                end
            end
            if ~isFoundMUA
                continue;
            end
        elseif ismember(EL.channelInds(j), R.vPulChannels)
            isFoundMUA = 0;
            offset = -1; % -1, 1, -2, 2, -3, 3, ...
            while ~isFoundMUA
                if j + offset > 0 && j + offset <= numel(EL.channelInds)
                    tempUnitName = sprintf('%s_PUL_%dM', sessionName, EL.channelInds(j + offset));
                    if any(strcmp(unitNamesVPul, tempUnitName)) % units that have significant cue response based on saved list
                        isFoundMUA = 1;
                        muaNearbyChannelTs = EL.allMUAStructs{j + offset}.ts;
                    end
                end
                if offset < 0
                    offset = -1 * offset;
                else
                    offset = -1 * offset - 1;
                end
                if abs(offset) == 32 % no possible units found
                    fprintf('No possible MUA found to match channel %d...\n', EL.channelInds(j));
                    break;
                end
            end
            if ~isFoundMUA
                continue;
            end
        end
        
        fprintf('Processing channel %d...\n', EL.channelInds(j));
        
        lfpCount = lfpCount + 1;
        lfpNames{lfpCount} = sprintf('%s_%d_FP%03d', sessionName, sessionInd, EL.channelInds(j));
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
%         if j < numel(EL.channelInds)-1
%             muaNearbyChannelTs = EL.allMUAStructs{j+2}.ts;
%             sfcNames{lfpCount} = sprintf('%s_%d_%dM_FP%03d', sessionName, sessionInd, j+2, EL.channelInds(j));
%         else
%             muaNearbyChannelTs = EL.allMUAStructs{j-2}.ts;
%             sfcNames{lfpCount} = sprintf('%s_%d_%dM_FP%03d', sessionName, sessionInd, j-2, EL.channelInds(j));
%         end

        
                    
                
            
            
        
        baselineInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, baselineWindowOffset);
        cueResponseInd = getTimeLogicalWithTolerance(EL.cueOnsetLfp.t, cueResponseOffset);
        cueTargetDelayInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, cueTargetDelayOffset);
        arrayResponseInd = getTimeLogicalWithTolerance(EL.arrayOnsetHoldBalLfp.t, arrayResponseOffset);
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
        
        numTrials = size(cueOnsetLfpCurrent, 1);
        numTrialsP3 = size(arrayOnsetLfpP3Current, 1);
        numTrialsP1 = size(arrayOnsetLfpP1Current, 1);
        numTrialsBalP3 = size(targetDimLfpP3Current, 1);
        numTrialsBalP1 = size(targetDimLfpP1Current, 1);
        
        preCueBaselineLfps = cueOnsetLfpCurrent(baselineInd,:);
        cueTargetDelayLfpsP3 = arrayOnsetLfpP3Current(cueTargetDelayInd,:);
        cueTargetDelayLfpsP1 = arrayOnsetLfpP1Current(cueTargetDelayInd,:);
        arrayResponseHoldLfpsP3 = arrayOnsetLfpP3Current(arrayResponseInd,:);
        arrayResponseHoldLfpsP1 = arrayOnsetLfpP1Current(arrayResponseInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);

        [baselinePowerLF(lfpCount,:),fAxisLF] = mtspectrumc(preCueBaselineLfps, paramsBaselineLF);
        [baselinePowerHF(lfpCount,:),fAxisHF] = mtspectrumc(preCueBaselineLfps, paramsBaselineHF);
        [cueTargetDelayPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP3, paramsLF);
        [cueTargetDelayPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP3, paramsHF);
        [cueTargetDelayPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(cueTargetDelayLfpsP1, paramsLF);
        [cueTargetDelayPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(cueTargetDelayLfpsP1, paramsHF);
        [arrayResponseHoldPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(arrayResponseHoldLfpsP3, paramsLF);
        [arrayResponseHoldPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(arrayResponseHoldLfpsP3, paramsHF);
        [arrayResponseHoldPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(arrayResponseHoldLfpsP1, paramsLF);
        [arrayResponseHoldPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(arrayResponseHoldLfpsP1, paramsHF);
        [targetDimDelayPowerP3LF(lfpCount,:),fAxisLF] = mtspectrumc(targetDimDelayLfpsP3, paramsLF);
        [targetDimDelayPowerP3HF(lfpCount,:),fAxisHF] = mtspectrumc(targetDimDelayLfpsP3, paramsHF);
        [targetDimDelayPowerP1LF(lfpCount,:),fAxisLF] = mtspectrumc(targetDimDelayLfpsP1, paramsLF);
        [targetDimDelayPowerP1HF(lfpCount,:),fAxisHF] = mtspectrumc(targetDimDelayLfpsP1, paramsHF);
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        baselineSFCLF(lfpCount,:) = atanh(C)-(1/((2*paramsBaselineLF.tapers(2)*numTrials)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        baselineSFCHF(lfpCount,:) = atanh(C)-(1/((2*paramsBaselineHF.tapers(2)*numTrials)-2)); % adjust for num trials
        
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
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{3}, arrayResponseOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        arrayResponseHoldSFCP3LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        arrayResponseHoldSFCP3HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
        
        alignedSpikeTs = createnonemptydatamatpt(muaNearbyChannelTs, EL.UE.arrayOnsetByLoc{1}, arrayResponseOffset .* [-1 1]);
        meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
        [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        arrayResponseHoldSFCP1LF(lfpCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
        [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
        if any(C(:) > 0.8), lfpCount = lfpCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
        arrayResponseHoldSFCP1HF(lfpCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
        
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
        
        % TODO across subdivision array response SFC
        
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
                    arrayResponseHoldLfps1P3 = arrayOnsetLfpP3Current1(arrayResponseInd,:); 
                    arrayResponseHoldLfps2P3 = arrayOnsetLfpP3Current2(arrayResponseInd,:); 
                    arrayResponseHoldLfps1P1 = arrayOnsetLfpP1Current1(arrayResponseInd,:);
                    arrayResponseHoldLfps2P1 = arrayOnsetLfpP1Current2(arrayResponseInd,:);
                    targetDimDelayLfps1P3 = targetDimLfpP3Current1(targetDimDelayInd,:); 
                    targetDimDelayLfps2P3 = targetDimLfpP3Current2(targetDimDelayInd,:); 
                    targetDimDelayLfps1P1 = targetDimLfpP1Current1(targetDimDelayInd,:);
                    targetDimDelayLfps2P1 = targetDimLfpP1Current2(targetDimDelayInd,:);

                    baselineSubPairCohLF(subPairCount,:) = coherencyc(preCueBaselineLfps1, preCueBaselineLfps2, paramsBaselineLF);
                    baselineSubPairCohHF(subPairCount,:) = coherencyc(preCueBaselineLfps1, preCueBaselineLfps2, paramsBaselineHF);
                    cueTargetDelaySubPairCohP3LF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P3, cueTargetDelayLfps2P3, paramsLF);
                    cueTargetDelaySubPairCohP3HF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P3, cueTargetDelayLfps2P3, paramsHF);
                    cueTargetDelaySubPairCohP1LF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P1, cueTargetDelayLfps2P1, paramsLF);
                    cueTargetDelaySubPairCohP1HF(subPairCount,:) = coherencyc(cueTargetDelayLfps1P1, cueTargetDelayLfps2P1, paramsHF);
                    arrayResponseHoldSubPairCohP3LF(subPairCount,:) = coherencyc(arrayResponseHoldLfps1P3, arrayResponseHoldLfps2P3, paramsLF);
                    arrayResponseHoldSubPairCohP3HF(subPairCount,:) = coherencyc(arrayResponseHoldLfps1P3, arrayResponseHoldLfps2P3, paramsHF);
                    arrayResponseHoldSubPairCohP1LF(subPairCount,:) = coherencyc(arrayResponseHoldLfps1P1, arrayResponseHoldLfps2P1, paramsLF);
                    arrayResponseHoldSubPairCohP1HF(subPairCount,:) = coherencyc(arrayResponseHoldLfps1P1, arrayResponseHoldLfps2P1, paramsHF);
                    targetDimDelaySubPairCohP3LF(subPairCount,:) = coherencyc(targetDimDelayLfps1P3, targetDimDelayLfps2P3, paramsLF);
                    targetDimDelaySubPairCohP3HF(subPairCount,:) = coherencyc(targetDimDelayLfps1P3, targetDimDelayLfps2P3, paramsHF);
                    targetDimDelaySubPairCohP1LF(subPairCount,:) = coherencyc(targetDimDelayLfps1P1, targetDimDelayLfps2P1, paramsLF);
                    targetDimDelaySubPairCohP1HF(subPairCount,:) = coherencyc(targetDimDelayLfps1P1, targetDimDelayLfps2P1, paramsHF);
                    
                    muaChannel1Ts = EL.allMUAStructs{j1}.ts;
                    muaChannel2Ts = EL.allMUAStructs{j2}.ts;
                    
                    % dPul spike, vPul field - one each per session
                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps2, alignedSpikeTs, paramsBaselineLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleDPulSpikeVPulFieldLF(subPairCount,:) = atanh(C)-(1/((2*paramsBaselineLF.tapers(2)*numTrials)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps2, alignedSpikeTs, paramsBaselineHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleDPulSpikeVPulFieldHF(subPairCount,:) = atanh(C)-(1/((2*paramsBaselineHF.tapers(2)*numTrials)-2)); % adjust for num trials

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
                    
                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.arrayOnsetByLoc{3}, arrayResponseOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfps2P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleDPulSpikeVPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfps2P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleDPulSpikeVPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel1Ts, EL.UE.arrayOnsetByLoc{1}, arrayResponseOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfps2P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleDPulSpikeVPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfps2P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleDPulSpikeVPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

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
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps1, alignedSpikeTs, paramsBaselineLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleVPulSpikeDPulFieldLF(subPairCount,:) = atanh(C)-(1/((2*paramsBaselineLF.tapers(2)*numTrials)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps1, alignedSpikeTs, paramsBaselineHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    baselineSFCSingleVPulSpikeDPulFieldHF(subPairCount,:) = atanh(C)-(1/((2*paramsBaselineHF.tapers(2)*numTrials)-2)); % adjust for num trials

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
                    
                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.arrayOnsetByLoc{3}, arrayResponseOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfps1P3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleVPulSpikeDPulFieldP3LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfps1P3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleVPulSpikeDPulFieldP3HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

                    alignedSpikeTs = createnonemptydatamatpt(muaChannel2Ts, EL.UE.arrayOnsetByLoc{1}, arrayResponseOffset .* [-1 1]);
                    meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
                    [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfps1P1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleVPulSpikeDPulFieldP1LF(subPairCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
                    [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfps1P1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
                    if any(C(:) > 0.8), subPairCount = subPairCount - 1; fprintf('Abnormally high coherence; skipping channel %d - %d\n', EL.channelInds([j1 j2])); continue; end; % skip this channel
                    arrayResponseHoldSFCSingleVPulSpikeDPulFieldP1HF(subPairCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

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
        arrayResponseHoldLfpsP3 = arrayOnsetLfpP3Current(arrayResponseInd,:);
        arrayResponseHoldLfpsP1 = arrayOnsetLfpP1Current(arrayResponseInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);
        
        for k = 1:numel(R.dPulChannels)
            tempUnitName = sprintf('%s_PUL_%dM', sessionName, R.dPulChannels(k));
            if ~any(strcmp(unitNamesDPul, tempUnitName)) % skip all units that do not have significant cue response based on saved list
                continue;
            end
            
            j = find(EL.channelInds == R.dPulChannels(k));
            fprintf('Processing channel %d...\n', EL.channelInds(j));
            spikeTs = EL.allMUAStructs{j}.ts;
            
            dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount + 1;
            dPulSpikeVPulFieldNames{dPulSpikeVPulFieldCount} = sprintf('%s_%d_%dM-vPulLfp', sessionName, sessionInd, EL.channelInds(j));
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCDPulSpikeVPulFieldLF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsBaselineLF.tapers(2)*numTrials)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCDPulSpikeVPulFieldHF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsBaselineHF.tapers(2)*numTrials)-2)); % adjust for num trials

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
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{3}, arrayResponseOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCDPulSpikeVPulFieldP3LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCDPulSpikeVPulFieldP3HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{1}, arrayResponseOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCDPulSpikeVPulFieldP1LF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), dPulSpikeVPulFieldCount = dPulSpikeVPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCDPulSpikeVPulFieldP1HF(dPulSpikeVPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

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
        arrayResponseHoldLfpsP3 = arrayOnsetLfpP3Current(arrayResponseInd,:);
        arrayResponseHoldLfpsP1 = arrayOnsetLfpP1Current(arrayResponseInd,:);
        targetDimDelayLfpsP3 = targetDimLfpP3Current(targetDimDelayInd,:);
        targetDimDelayLfpsP1 = targetDimLfpP1Current(targetDimDelayInd,:);
        
        for k = 1:numel(R.vPulChannels)
            tempUnitName = sprintf('%s_PUL_%dM', sessionName, R.vPulChannels(k));
            if ~any(strcmp(unitNamesVPul, tempUnitName)) % skip all units that do not have significant cue response based on saved list
                continue;
            end
            
            j = find(EL.channelInds == R.vPulChannels(k));
            fprintf('Processing channel %d...\n', EL.channelInds(j));
            spikeTs = EL.allMUAStructs{j}.ts;
            
            vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount + 1;
            vPulSpikeDPulFieldNames{vPulSpikeDPulFieldCount} = sprintf('%s_%d_%dM-dPulLfp', sessionName, sessionInd, EL.channelInds(j));
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.cueOnset, baselineWindowOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, baselineWindowOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCVPulSpikeDPulFieldLF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsBaselineLF.tapers(2)*numTrials)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(preCueBaselineLfps, alignedSpikeTs, paramsBaselineHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            baselineSFCVPulSpikeDPulFieldHF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsBaselineHF.tapers(2)*numTrials)-2)); % adjust for num trials

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
            
            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{3}, arrayResponseOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCVPulSpikeDPulFieldP3LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP3)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP3, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCVPulSpikeDPulFieldP3HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP3)-2)); % adjust for num trials

            alignedSpikeTs = createnonemptydatamatpt(spikeTs, EL.UE.arrayOnsetByLoc{1}, arrayResponseOffset .* [-1 1]);
            meanFR = computeMeanFiringRateFromSpikeTimesMat(alignedSpikeTs, arrayResponseOffset);
            [C,~,~,~,~,fAxisLF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsLF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCVPulSpikeDPulFieldP1LF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsLF.tapers(2)*numTrialsP1)-2)); % adjust for num trials
            [C,~,~,~,~,fAxisHF] = Adjcoherencycpt_faster(arrayResponseHoldLfpsP1, alignedSpikeTs, paramsHF, 0, [], adjCohNormRate, meanFR);
            if any(C(:) > 0.8), vPulSpikeDPulFieldCount = vPulSpikeDPulFieldCount - 1; fprintf('Abnormally high coherence; skipping channel %d\n', EL.channelInds(j)); continue; end; % skip this channel
            arrayResponseHoldSFCVPulSpikeDPulFieldP1HF(vPulSpikeDPulFieldCount,:) = atanh(C)-(1/((2*paramsHF.tapers(2)*numTrialsP1)-2)); % adjust for num trials

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

%% plot power in array response dPul P3 vs P1
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelPowerLF = [arrayResponseHoldPowerP3LF(channelCond,:); arrayResponseHoldPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
arrayResponseHoldRelPowerHF = [arrayResponseHoldPowerP3HF(channelCond,:); arrayResponseHoldPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(arrayResponseHoldRelPowerLF, arrayResponseHoldRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Array Response Relative Power - Dorsal Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldRelPowerDB-dPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% plot power in array response vPul P3 vs P1
channelCond = isInVPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelPowerLF = [arrayResponseHoldPowerP3LF(channelCond,:); arrayResponseHoldPowerP1LF(channelCond,:)] ./ ...
        [baselinePowerLF(channelCond,:); baselinePowerLF(channelCond,:)];
arrayResponseHoldRelPowerHF = [arrayResponseHoldPowerP3HF(channelCond,:); arrayResponseHoldPowerP1HF(channelCond,:)] ./ ...
        [baselinePowerHF(channelCond,:); baselinePowerHF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(arrayResponseHoldRelPowerLF, arrayResponseHoldRelPowerHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [relPowYBounds; relPowYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Power Rel. to Baseline (dB/Hz)', ...
        'titleText', sprintf('Array Response Relative Power - Ventral Pulvinar (%s)', ref), ...
        'doDB', 1);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldRelPowerDB-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
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

%% plot array response SFC dPul P3 vs P1 (relative to baseline is NOISY)
channelCond = isInDPulvinar;
nChannel = sum(channelCond);
arrayResponseHoldRelSFCLF = [arrayResponseHoldSFCP3LF(channelCond,:); arrayResponseHoldSFCP1LF(channelCond,:)];
arrayResponseHoldRelSFCHF = [arrayResponseHoldSFCP3HF(channelCond,:); arrayResponseHoldSFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Array Response SFC - Dorsal Pulvinar (%s)', ref), ...
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

plotLfpPower2(arrayResponseHoldRelSFCLF, arrayResponseHoldRelSFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Array Response SFC - Ventral Pulvinar (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-arrayResponseHoldSFC-vPul-P3vsP1-%s-v%d.png', outputDir, ref, v);
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

%% plot subdivision pair lfp coherence across conditions and periods
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

yBounds = [0 0.12];
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

yBounds = [0 0.12];
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

%% make tiny population plots of LFP power by channel
% TODO convert to db
powYBounds = [0 0.02];

powP3 = cueTargetDelayPowerP3LF(isInDPulvinar,:);
powP1 = cueTargetDelayPowerP1LF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelayPower-dPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        powYBounds, sprintf('Cue-Target Delay Power - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = cueTargetDelayPowerP3LF(isInVPulvinar,:);
powP1 = cueTargetDelayPowerP1LF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelayPower-vPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        powYBounds, sprintf('Cue-Target Delay Power - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelayPowerP3LF(isInDPulvinar,:);
powP1 = targetDimDelayPowerP1LF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelayPower-dPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        powYBounds, sprintf('Target-Dim Delay Power - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelayPowerP3LF(isInVPulvinar,:);
powP1 = targetDimDelayPowerP1LF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelayPower-vPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        powYBounds, sprintf('Target-Dim Delay Power - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

powYBounds = [0 0.01];

powP3 = cueTargetDelayPowerP3HF(isInDPulvinar,:);
powP1 = cueTargetDelayPowerP1HF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelayPower-dPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        powYBounds, sprintf('Cue-Target Delay Power - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = cueTargetDelayPowerP3HF(isInVPulvinar,:);
powP1 = cueTargetDelayPowerP1HF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelayPower-vPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        powYBounds, sprintf('Cue-Target Delay Power - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelayPowerP3HF(isInDPulvinar,:);
powP1 = targetDimDelayPowerP1HF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelayPower-dPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        powYBounds, sprintf('Target-Dim Delay Power - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelayPowerP3HF(isInVPulvinar,:);
powP1 = targetDimDelayPowerP1HF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelayPower-vPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        powYBounds, sprintf('Target-Dim Delay Power - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
%% make tiny population plots of SFC by channel
sfcYBounds = [0 0.3];

powP3 = baselineSFCLF(isInDPulvinar,:);
powP1 = nan(size(powP3));
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-dPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = baselineSFCLF(isInVPulvinar,:);
powP1 = nan(size(powP3));
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-vPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

powP3 = cueTargetDelaySFCP3LF(isInDPulvinar,:);
powP1 = cueTargetDelaySFCP1LF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-dPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = cueTargetDelaySFCP3LF(isInVPulvinar,:);
powP1 = cueTargetDelaySFCP1LF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-vPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelaySFCP3LF(isInDPulvinar,:);
powP1 = targetDimDelaySFCP1LF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-dPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelaySFCP3LF(isInVPulvinar,:);
powP1 = targetDimDelaySFCP1LF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-vPul-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

powP3 = baselineSFCHF(isInDPulvinar,:);
powP1 = nan(size(powP3));
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-dPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = baselineSFCHF(isInVPulvinar,:);
powP1 = nan(size(powP3));
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-vPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = cueTargetDelaySFCP3HF(isInDPulvinar,:);
powP1 = cueTargetDelaySFCP1HF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-dPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = cueTargetDelaySFCP3HF(isInVPulvinar,:);
powP1 = cueTargetDelaySFCP1HF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-vPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelaySFCP3HF(isInDPulvinar,:);
powP1 = targetDimDelaySFCP1HF(isInDPulvinar,:);
chanNames = lfpNames(isInDPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-dPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - Dorsal Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
powP3 = targetDimDelaySFCP3HF(isInVPulvinar,:);
powP1 = targetDimDelaySFCP1HF(isInVPulvinar,:);
chanNames = lfpNames(isInVPulvinar);

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-vPul-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(powP3, zeros(size(powP3)), powP1, zeros(size(powP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - Ventral Pulvinar (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);

%% make tiny population plots of across subdivision SFC
sfcYBounds = [0 0.3];

sfcP3 = baselineSFCDPulSpikeVPulFieldLF;
sfcP1 = nan(size(sfcP3));
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-dPulSpikeVPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

sfcP3 = baselineSFCVPulSpikeDPulFieldLF;
sfcP1 = nan(size(sfcP3));
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-vPulSpikeDPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

sfcP3 = cueTargetDelaySFCDPulSpikeVPulFieldP3LF;
sfcP1 = cueTargetDelaySFCDPulSpikeVPulFieldP1LF;
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-dPulSpikeVPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

sfcP3 = cueTargetDelaySFCVPulSpikeDPulFieldP3LF;
sfcP1 = cueTargetDelaySFCVPulSpikeDPulFieldP1LF;
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-vPulSpikeDPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);
    
sfcP3 = targetDimDelaySFCDPulSpikeVPulFieldP3LF;
sfcP1 = targetDimDelaySFCDPulSpikeVPulFieldP1LF;
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-dPulSpikeVPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

sfcP3 = targetDimDelaySFCVPulSpikeDPulFieldP3LF;
sfcP1 = targetDimDelaySFCVPulSpikeDPulFieldP1LF;
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-vPulSpikeDPulField-P3VsP1-LF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisLF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisLF([1 end])), ref), plotFileBaseName);

sfcP3 = baselineSFCDPulSpikeVPulFieldHF;
sfcP1 = nan(size(sfcP3));
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-dPulSpikeVPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);

sfcP3 = baselineSFCVPulSpikeDPulFieldHF;
sfcP1 = nan(size(sfcP3));
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-baselineSFC-vPulSpikeDPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Pre-Cue Baseline SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
sfcP3 = cueTargetDelaySFCDPulSpikeVPulFieldP3HF;
sfcP1 = cueTargetDelaySFCDPulSpikeVPulFieldP1HF;
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-dPulSpikeVPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);

sfcP3 = cueTargetDelaySFCVPulSpikeDPulFieldP3HF;
sfcP1 = cueTargetDelaySFCVPulSpikeDPulFieldP1HF;
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-vPulSpikeDPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);
    
sfcP3 = targetDimDelaySFCDPulSpikeVPulFieldP3HF;
sfcP1 = targetDimDelaySFCDPulSpikeVPulFieldP1HF;
chanNames = dPulSpikeVPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-dPulSpikeVPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);

sfcP3 = targetDimDelaySFCVPulSpikeDPulFieldP3HF;
sfcP1 = targetDimDelaySFCVPulSpikeDPulFieldP1HF;
chanNames = vPulSpikeDPulFieldNames;

plotFileBaseName = sprintf('%s/allSessions-tdDelaySFC-vPulSpikeDPulField-P3VsP1-HF-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Target-Dim Delay SFC - vPul Spikes - dPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);

%%
return; % end plots, rest is live testing

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCLF' < 0.05) | all(baselineSFCLF' < 0.03));

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
goodPairs = ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

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

%% only plot across-subdiv SFC for pairs with nonflat baseline SFC

% goodPairs = true(size(baselineSFCVPulSpikeDPulFieldLF, 1), 1);
goodPairs = ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

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
goodPairs = ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

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
goodPairs = ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

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
goodPairs = ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

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
goodPairs = ~(all(baselineSFCVPulSpikeDPulFieldLF' < 0.05) | all(baselineSFCVPulSpikeDPulFieldHF' < 0.03));

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
goodPairs = ~(all(baselineSFCDPulSpikeVPulFieldLF' < 0.05) | all(baselineSFCDPulSpikeVPulFieldHF' < 0.03));

sfcP3 = cueTargetDelaySFCDPulSpikeVPulFieldP3HF(goodPairs,:);
sfcP1 = cueTargetDelaySFCDPulSpikeVPulFieldP1HF(goodPairs,:);
chanNames = dPulSpikeVPulFieldNames(goodPairs,:);

plotFileBaseName = sprintf('%s/allSessions-ctDelaySFC-dPulSpikeVPulField-P3VsP1-HF-baselineRefined-%s-v%d', outputDir, ref, v);
makeTinyPlotsOfPopulationPower(sfcP3, zeros(size(sfcP3)), sfcP1, zeros(size(sfcP1)), fAxisHF, chanNames, ...
        sfcYBounds, sprintf('Cue-Target Delay SFC - dPul Spikes - vPul Field (%d-%d Hz, %s)', round(fAxisHF([1 end])), ref), plotFileBaseName);




%% inspect cue-target power in vPul at 8-12 Hz
fBounds = [8 12];
fInd = fAxisLF >= fBounds(1) & fAxisLF <= fBounds(2); 

powP3 = mean(cueTargetDelayPowerP3LF(isInVPulvinar,fInd), 2); % mean over f
powP1 = mean(cueTargetDelayPowerP1LF(isInVPulvinar,fInd), 2); % mean over f
powDiff = powP3 - powP1;

histBins = 0:0.002:0.02;

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
xlim([0 0.02]);
ylim([0 0.02]);

[h,p,stats] = ttest(powDiff);
fprintf('Median diff power = %0.3f (N = %d), p = %0.3f\n', median(powDiff), size(powDiff, 1), p);

%% inspect cue-target power in vPul at 75-80 Hz
fBounds = [75 80];
fInd = fAxisHF >= fBounds(1) & fAxisHF <= fBounds(2); 

powP3 = mean(cueTargetDelayPowerP3HF(isInVPulvinar,fInd), 2); % mean over f
powP1 = mean(cueTargetDelayPowerP1HF(isInVPulvinar,fInd), 2); % mean over f
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

sfcP3 = mean(cueTargetDelaySFCP3HF(isInDPulvinar,fInd), 2); % mean over f
sfcP1 = mean(cueTargetDelaySFCP1HF(isInDPulvinar,fInd), 2); % mean over f
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

sfcP3 = mean(cueTargetDelaySFCP3LF(isInVPulvinar,fInd), 2); % mean over f
sfcP1 = mean(cueTargetDelaySFCP1LF(isInVPulvinar,fInd), 2); % mean over f
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


%% inspect target-dim all pairs vPul spike - dPul field coherence at 27-35 Hz
fBounds = [27 35];
fInd = fAxisHF >= fBounds(1) & fAxisHF <= fBounds(2); 

sfcP3 = mean(targetDimDelaySFCVPulSpikeDPulFieldP3HF(:,fInd), 2); % mean over f
sfcP1 = mean(targetDimDelaySFCVPulSpikeDPulFieldP1HF(:,fInd), 2); % mean over f
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

%% mean SFC of units with strong attentional modulation in ct delay
lfpNamesStrongMUACTDelay = {
        'M20170311_7_FP004' %'M20170311_PUL_6M'
        'M20170327_13_FP009' %'M20170327_PUL_11M'
        'M20170311_8_FP046' %'M20170311_PUL_48M'
        'M20170311_8_FP043' %'M20170311_PUL_45M'
        'M20170327_14_FP046' %'M20170327_PUL_48M'
        'M20170331_22_FP039' %'M20170331_PUL_41M'
        'M20170311_7_FP005' %'M20170311_PUL_7M'
        'M20170327_14_FP047' %'M20170327_PUL_49M'
        'M20170327_14_FP052' %'M20170327_PUL_54M'
        'M20170211_4_FP018' %'M20170211_PUL_20M'
        'M20170211_4_FP014' %'M20170211_PUL_16M'
        'M20170327_14_FP045' %'M20170327_PUL_47M'
        'M20170331_22_FP034' %'M20170331_PUL_36M'
%         'M20170331_22_FP032'
%         'M20170329_PUL_FP011'
% vPul below
        'M20170311_8_FP052' %'M20170311_PUL_52M'
        'M20170308_5_FP022' %'M20170308_PUL_22M'
        'M20170320_10_FP043' %'M20170320_PUL_43M'
        'M20170324_12_FP048' %'M20170324_PUL_48M'
        'M20170311_8_FP058' %'M20170311_PUL_58M'
        'M20170311_8_FP054' %'M20170311_PUL_54M'
        'M20170127_1_FP005' %'M20170127_PUL_5M'
        'M20170308_6_FP050' %'M20170308_PUL_50M'
        'M20170308_6_FP053' %'M20170308_PUL_53M'
        'M20170407_23_FP010' %'M20170407_PUL_10M'
        'M20170311_8_FP060' %'M20170311_PUL_60M'
        'M20170311_8_FP055' %'M20170311_PUL_55M'
        'M20170311_8_FP053' %'M20170311_PUL_53M'
        'M20170331_20_FP044' %'M20170331_PUL_44M'
        'M20170201_3_FP032' %'M20170201_PUL_32M'
        };

matchInds = nan(numel(lfpNamesStrongMUACTDelay), 1);
for i = 1:numel(lfpNamesStrongMUACTDelay)
    ind = find(strcmp(lfpNames, lfpNamesStrongMUACTDelay{i}));
    if ind > 0
        matchInds(i) = ind;
    end
end
matchInds(isnan(matchInds)) = [];

matchLog = false(size(lfpNames));
matchLog(matchInds) = 1;

channelCond = matchLog;
nChannel = sum(channelCond);
cueTargetDelaySFCLF = [cueTargetDelaySFCP3LF(channelCond,:); cueTargetDelaySFCP1LF(channelCond,:)];
cueTargetDelaySFCHF = [cueTargetDelaySFCP3HF(channelCond,:); cueTargetDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(cueTargetDelaySFCLF, cueTargetDelaySFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcCTDelayYBounds; sfcCTDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Cue-Target Delay SFC - Units with Most Attn Modulation in FR (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-ctDelaySFC-top15MUAMod-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');

%% mean SFC of units with strong attentional modulation in td delay
lfpNamesStrongMUATDDelay = {
        'M20170211_4_FP020' %'M20170211_PUL_20M'
        'M20170308_5_FP011' %'M20170308_PUL_11M'
        'M20170407_23_FP001' %'M20170407_PUL_1M'
%         'M20170329_P_17' %'M20170329_PUL_17M'
        'M20170211_4_FP021' %'M20170211_PUL_21M'
        'M20170324_11_FP003' %'M20170324_PUL_3M'
        'M20170311_8_FP045' %'M20170311_PUL_45M'
        'M20170331_20_FP039' %'M20170331_PUL_39M' % or 22
        'M20170308_5_FP008' %'M20170308_PUL_8M'
        'M20170327_14_FP054' %'M20170327_PUL_54M'
        'M20170331_20_FP033' %'M20170331_PUL_33M'
        'M20170320_9_FP015' %'M20170320_PUL_15M'
        'M20170308_8_FP044' %'M20170308_PUL_44M'
%         'M20170329_PUL_13' %'M20170329_PUL_13M'
        'M20170327_13_FP011' %'M20170327_PUL_11M'
        % vPul below
        'M20170320_10_FP049' %'M20170320_PUL_49M'
        'M20170311_8_FP060' %'M20170311_PUL_60M'
        'M20170308_6_FP053' %'M20170308_PUL_53M'
        'M20170320_10_FP048' %'M20170320_PUL_48M'
        'M20170311_8_FP054' %'M20170311_PUL_54M'
%         'M20170329_PUL_26' %'M20170329_PUL_26M'
        'M20170407_23_FP010' %'M20170407_PUL_10M'
        'M20170311_8_FP061' %'M20170311_PUL_61M'
        'M20170331_20_FP046' %'M20170331_PUL_46M' % or 22
        'M20170127_1_FP005' %'M20170127_PUL_5M'
%         'M20170329_PUL_26' %'M20170329_PUL_26M'
        'M20170331_20_FP047' %'M20170331_PUL_47M'
        'M20170311_8_FP053' %'M20170311_PUL_53M'
        'M20170327_13_FP013' %'M20170327_PUL_13M'
%         'M20170329_PUL_23' %'M20170329_PUL_23M'
        };

matchInds = nan(numel(lfpNamesStrongMUATDDelay), 1);
for i = 1:numel(lfpNamesStrongMUATDDelay)
    ind = find(strcmp(lfpNames, lfpNamesStrongMUATDDelay{i}));
    if ind > 0
        matchInds(i) = ind;
    end
end
matchInds(isnan(matchInds)) = [];

matchLog = false(size(lfpNames));
matchLog(matchInds) = 1;

channelCond = matchLog;
nChannel = sum(channelCond);
targetDimDelaySFCLF = [targetDimDelaySFCP3LF(channelCond,:); targetDimDelaySFCP1LF(channelCond,:)];
targetDimDelaySFCHF = [targetDimDelaySFCP3HF(channelCond,:); targetDimDelaySFCP1HF(channelCond,:)];
p3p1Logical = [true(nChannel, 1) false(nChannel, 1); false(nChannel, 1) true(nChannel, 1)];

plotLfpPower2(targetDimDelaySFCLF, targetDimDelaySFCHF, fAxisLF, fAxisHF, p3p1Logical, ...
        'xBounds', [paramsLF.fpass; paramsHF.fpass], ...
        'yBounds', [sfcTDDelayYBounds; sfcTDDelayYBounds], ...
        'cols', [p3Col; p1Col], ...
        'lineLabels', {'P3', 'P1'}, ...
        'ylabelText', 'Coherence', ...
        'titleText', sprintf('Target-Dim Delay SFC - Units with Most Attn Modulation in FR (%s)', ref), ...
        'doDB', 0);

plotFileName = sprintf('%s/allSessions-tdDelaySFC-top15MUAMod-P3vsP1-%s-v%d.png', outputDir, ref, v);
fprintf('Saving to %s...\n', plotFileName);
export_fig(plotFileName, '-nocrop');
