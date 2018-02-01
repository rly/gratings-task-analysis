function EL = computeEvokedLfps(saveFileName, lfp, Fs, nLoc, UE)
% lfp = 1-D vector
% UE = useful events struct
% assert(iscolumn(lfp));

nTrials = numel(UE.cueOnset);

%% align spikes to cue onset
cueOnsetWindowOffset = [-0.7 0.7];
[cueOnsetLfp,cueOnsetT] = alignLfpToEvents(lfp, UE.cueOnset, Fs, cueOnsetWindowOffset);

% don't split by location -- just index into cueOnsetLfp with logical check
% for location -- this saves disk space by 2x
% cueOnsetLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     cueOnsetLfpByLoc{i} = alignLfpToEvents(lfp, UE.cueOnsetByLoc{i}, Fs, cueOnsetWindowOffset);
% end

preCueBaselineWindowOffset = [-0.3 0];
cueResponseWindowOffset = [0.025 0.2];

%% align spikes to array onset
arrayOnsetWindowOffset = [-0.7 0.7];
[arrayOnsetLfp,arrayOnsetT] = alignLfpToEvents(lfp, UE.arrayOnset, Fs, arrayOnsetWindowOffset);
[arrayOnsetRelLfp,~] = alignLfpToEvents(lfp, UE.arrayOnsetRel, Fs, arrayOnsetWindowOffset);
[arrayOnsetHoldLfp,~] = alignLfpToEvents(lfp, UE.arrayOnsetHold, Fs, arrayOnsetWindowOffset);

% arrayOnsetLfpByLoc = cell(nLoc, 1);
% arrayOnsetRelLfpByLoc = cell(nLoc, 1);
% arrayOnsetHoldLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     arrayOnsetLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetByLoc{i}, Fs, arrayOnsetWindowOffset);
%     arrayOnsetRelLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetRelByLoc{i}, Fs, arrayOnsetWindowOffset);
%     arrayOnsetHoldLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetHoldByLoc{i}, Fs, arrayOnsetWindowOffset);
% end

cueTargetDelayWindowOffset = [-0.175 0];
cueTargetDelayLongWindowOffset = [-0.4 0];
arrayResponseWindowOffset = [0.025 0.2];

arrayOnsetRelToJuiceEventTime = UE.firstJuiceEvent(~UE.isHoldTrial) - UE.arrayOnsetRel; % almost same as UE.rt(~isHoldTrial)?


%% align spikes to target dimming
targetDimWindowOffset = [-0.7 0.7];
[targetDimLfp,targetDimT] = alignLfpToEvents(lfp, UE.targetDim, Fs, targetDimWindowOffset);

% targetDimLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     targetDimLfpByLoc{i} = alignLfpToEvents(lfp, UE.targetDimByLoc{i}, Fs, targetDimWindowOffset);
% end

targetDimDelayWindowOffset = [-0.175 0];
targetDimDelayLongWindowOffset = [-0.4 0];
targetDimResponseWindowOffset = [0.025 0.2]; 
postTargetDimMotorResponseWindowOffset = [0.425 0.6];

targetDimToJuiceEventTime = UE.firstJuiceEvent(UE.isHoldTrial) - UE.targetDim; % almost same as UE.rt(isHoldTrial)?

%% look at lever-locked responses and saccade-locked responses
% note that this includes both hold and release trials
assert(~isempty(UE.fixationAndLeverTimes))

%% align spikes to enter/exit fixation
enterFixationWindowOffset = [-0.7 0.7];
[enterFixationLfp,enterFixationT] = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue, Fs, enterFixationWindowOffset);

% enterFixationLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     enterFixationLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCueByLoc{i}, Fs, enterFixationWindowOffset);
% end

exitFixationWindowOffset = [-0.7 0.7];
[exitFixationLfp,exitFixationT] = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice, Fs, exitFixationWindowOffset);

% exitFixationLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     exitFixationLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuiceByLoc{i}, Fs, exitFixationWindowOffset);
% end

preEnterFixationWindowOffset = [-0.175 0];
postEnterFixationWindowOffset = [0.025 0.2];
postEnterFixationLateWindowOffset = [0.15 0.325];
preExitFixationWindowOffset = [-0.175 0];
postExitFixationWindowOffset = [0.025 0.2];

%% align spikes to lever press and release
leverPressWindowOffset = [-0.7 0.7];
[leverPressLfp,leverPressT] = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverPressTimesPreCue, Fs, leverPressWindowOffset);

% leverPressLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     leverPressLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverPressTimesPreCueByLoc{i}, Fs, leverPressWindowOffset);
% end

leverReleaseWindowOffset = [-0.7 0.7];
[leverReleaseLfp,leverReleaseT] = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice, Fs, leverReleaseWindowOffset);

% leverReleaseLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     leverReleaseLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuiceByLoc{i}, Fs, leverReleaseWindowOffset);
% end

preLeverReleaseWindowOffset = [-0.175 0];
postLeverReleaseWindowOffset = [0.025 0.2];
preLeverReleaseWindowOffset = [-0.175 0];
postLeverReleaseWindowOffset = [0.025 0.2];

%% calc time between motor events
assert(numel(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue) == numel(UE.fixationAndLeverTimes.firstLeverPressTimesPreCue));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice));

enterFixationToLeverPressTime = UE.fixationAndLeverTimes.firstLeverPressTimesPreCue - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
enterFixationToCueOnsetTime = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
exitFixationToJuiceEventTime = UE.firstJuiceEvent - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
exitFixationToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
arrayOnsetRelToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;
arrayOnsetRelToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;

%%
clear lfp;
tic;
fprintf('Saving to file: %s...\n', saveFileName);
save(saveFileName);
fileInfo = dir(saveFileName);
fprintf('%d MB file created. Took %0.1f minutes.\n', round(fileInfo.bytes/1024^2), toc/60);

end
% 
% %% average number of spikes in each window
% averageFiringRatesByCount = struct();
% averageFiringRatesByCount.preEnterFixation = computeAverageFiringRateByCount(...
%         preEnterFixationWindowOffset, enterFixationWindow(1), enterFixationSpikeTimes, enterFixationSpikeTimesByLoc);
% averageFiringRatesByCount.postEnterFixation = computeAverageFiringRateByCount(...
%         postEnterFixationWindowOffset, enterFixationWindow(1), enterFixationSpikeTimes, enterFixationSpikeTimesByLoc);
% averageFiringRatesByCount.postEnterFixationLate = computeAverageFiringRateByCount(...
%         postEnterFixationLateWindowOffset, enterFixationWindow(1), enterFixationSpikeTimes, enterFixationSpikeTimesByLoc);
% averageFiringRatesByCount.preCueBaseline = computeAverageFiringRateByCount(...
%         preCueBaselineWindowOffset, cueOnsetWindow(1), cueOnsetSpikeTimes, cueOnsetSpikeTimesByLoc);
% averageFiringRatesByCount.cueResponse = computeAverageFiringRateByCount(...
%         cueResponseWindowOffset, cueOnsetWindow(1), cueOnsetSpikeTimes, cueOnsetSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelay = computeAverageFiringRateByCount(...
%         cueTargetDelayWindowOffset, arrayOnsetWindow(1), arrayOnsetSpikeTimes, arrayOnsetSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelayRel = computeAverageFiringRateByCount(...
%         cueTargetDelayWindowOffset, arrayOnsetWindow(1), arrayOnsetRelSpikeTimes, arrayOnsetRelSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelayHold = computeAverageFiringRateByCount(...
%         cueTargetDelayWindowOffset, arrayOnsetWindow(1), arrayOnsetHoldSpikeTimes, arrayOnsetHoldSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelayLong = computeAverageFiringRateByCount(...
%         cueTargetDelayLongWindowOffset, arrayOnsetWindow(1), arrayOnsetSpikeTimes, arrayOnsetSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelayLongRel = computeAverageFiringRateByCount(...
%         cueTargetDelayLongWindowOffset, arrayOnsetWindow(1), arrayOnsetRelSpikeTimes, arrayOnsetRelSpikeTimesByLoc);
% averageFiringRatesByCount.cueTargetDelayLongHold = computeAverageFiringRateByCount(...
%         cueTargetDelayLongWindowOffset, arrayOnsetWindow(1), arrayOnsetHoldSpikeTimes, arrayOnsetHoldSpikeTimesByLoc);
% averageFiringRatesByCount.arrayResponse = computeAverageFiringRateByCount(...
%         arrayResponseWindowOffset, arrayOnsetWindow(1), arrayOnsetSpikeTimes, arrayOnsetSpikeTimesByLoc);
% averageFiringRatesByCount.arrayRelResponse = computeAverageFiringRateByCount(...
%         arrayResponseWindowOffset, arrayOnsetWindow(1), arrayOnsetRelSpikeTimes, arrayOnsetRelSpikeTimesByLoc);
% averageFiringRatesByCount.arrayHoldResponse = computeAverageFiringRateByCount(...
%         arrayResponseWindowOffset, arrayOnsetWindow(1), arrayOnsetHoldSpikeTimes, arrayOnsetHoldSpikeTimesByLoc);
% averageFiringRatesByCount.targetDimDelay = computeAverageFiringRateByCount(...
%         targetDimDelayWindowOffset, targetDimWindow(1), targetDimSpikeTimes, targetDimSpikeTimesByLoc);
% averageFiringRatesByCount.targetDimDelayLong = computeAverageFiringRateByCount(...
%         targetDimDelayLongWindowOffset, targetDimWindow(1), targetDimSpikeTimes, targetDimSpikeTimesByLoc);
% averageFiringRatesByCount.targetDimResponse = computeAverageFiringRateByCount(...
%         targetDimResponseWindowOffset, targetDimWindow(1), targetDimSpikeTimes, targetDimSpikeTimesByLoc);
% averageFiringRatesByCount.targetDimShortHoldResponse = computeAverageFiringRateByCount(...
%         targetDimResponseWindowOffset, targetDimWindow(1), targetDimShortHoldSpikeTimes, targetDimShortHoldSpikeTimesByLoc);
% averageFiringRatesByCount.targetDimLongHoldResponse = computeAverageFiringRateByCount(...
%         targetDimResponseWindowOffset, targetDimWindow(1), targetDimLongHoldSpikeTimes, targetDimLongHoldSpikeTimesByLoc);
% averageFiringRatesByCount.postTargetDimMotorResponse = computeAverageFiringRateByCount(...
%         postTargetDimMotorResponseWindowOffset, targetDimWindow(1), targetDimSpikeTimes, targetDimSpikeTimesByLoc);
% averageFiringRatesByCount.preExitFixation = computeAverageFiringRateByCount(...
%         preExitFixationWindowOffset, exitFixationWindow(1), exitFixationSpikeTimes, exitFixationSpikeTimesByLoc);
% averageFiringRatesByCount.postExitFixation = computeAverageFiringRateByCount(...
%         postExitFixationWindowOffset, exitFixationWindow(1), exitFixationSpikeTimes, exitFixationSpikeTimesByLoc);
% 
% %% average spdf value in each window
% averageFiringRatesBySpdf = struct();
% averageFiringRatesBySpdf.kernelSigma = kernelSigma;
% averageFiringRatesBySpdf.preEnterFixation = computeAverageFiringRateBySpdf(...
%         preEnterFixationWindowOffset, enterFixationT, enterFixationWindow(1), enterFixationSpdf, enterFixationSpdfByLoc);
% averageFiringRatesBySpdf.postEnterFixation = computeAverageFiringRateBySpdf(...
%         postEnterFixationWindowOffset, enterFixationT, enterFixationWindow(1), enterFixationSpdf, enterFixationSpdfByLoc);
% averageFiringRatesBySpdf.postEnterFixationLate = computeAverageFiringRateBySpdf(...
%         postEnterFixationLateWindowOffset, enterFixationT, enterFixationWindow(1), enterFixationSpdf, enterFixationSpdfByLoc);
% averageFiringRatesBySpdf.preCueBaseline = computeAverageFiringRateBySpdf(...
%         preCueBaselineWindowOffset, cueOnsetT, cueOnsetWindow(1), cueOnsetSpdf, cueOnsetSpdfByLoc);
% averageFiringRatesBySpdf.cueResponse = computeAverageFiringRateBySpdf(...
%         cueResponseWindowOffset, cueOnsetT, cueOnsetWindow(1), cueOnsetSpdf, cueOnsetSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelay = computeAverageFiringRateBySpdf(...
%         cueTargetDelayWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetSpdf, arrayOnsetSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelayRel = computeAverageFiringRateBySpdf(...
%         cueTargetDelayWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetRelSpdf, arrayOnsetRelSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelayHold = computeAverageFiringRateBySpdf(...
%         cueTargetDelayWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetHoldSpdf, arrayOnsetHoldSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelayLong = computeAverageFiringRateBySpdf(...
%         cueTargetDelayLongWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetSpdf, arrayOnsetSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelayLongRel = computeAverageFiringRateBySpdf(...
%         cueTargetDelayLongWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetRelSpdf, arrayOnsetRelSpdfByLoc);
% averageFiringRatesBySpdf.cueTargetDelayLongHold = computeAverageFiringRateBySpdf(...
%         cueTargetDelayLongWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetHoldSpdf, arrayOnsetHoldSpdfByLoc);
% averageFiringRatesBySpdf.arrayResponse = computeAverageFiringRateBySpdf(...
%         arrayResponseWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetSpdf, arrayOnsetSpdfByLoc);
% averageFiringRatesBySpdf.arrayRelResponse = computeAverageFiringRateBySpdf(...
%         arrayResponseWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetRelSpdf, arrayOnsetRelSpdfByLoc);
% averageFiringRatesBySpdf.arrayHoldResponse = computeAverageFiringRateBySpdf(...
%         arrayResponseWindowOffset, arrayOnsetT, arrayOnsetWindow(1), arrayOnsetHoldSpdf, arrayOnsetHoldSpdfByLoc);
% averageFiringRatesBySpdf.targetDimDelay = computeAverageFiringRateBySpdf(...
%         targetDimDelayWindowOffset, targetDimT, targetDimWindow(1), targetDimSpdf, targetDimSpdfByLoc);
% averageFiringRatesBySpdf.targetDimDelayLong = computeAverageFiringRateBySpdf(...
%         targetDimDelayLongWindowOffset, targetDimT, targetDimWindow(1), targetDimSpdf, targetDimSpdfByLoc);
% averageFiringRatesBySpdf.targetDimResponse = computeAverageFiringRateBySpdf(...
%         targetDimResponseWindowOffset, targetDimT, targetDimWindow(1), targetDimSpdf, targetDimSpdfByLoc);
% averageFiringRatesBySpdf.targetDimShortHoldResponse = computeAverageFiringRateBySpdf(...
%         targetDimResponseWindowOffset, targetDimT, targetDimWindow(1), targetDimShortHoldSpdf, targetDimShortHoldSpdfByLoc);
% averageFiringRatesBySpdf.targetDimLongHoldResponse = computeAverageFiringRateBySpdf(...
%         targetDimResponseWindowOffset, targetDimT, targetDimWindow(1), targetDimLongHoldSpdf, targetDimLongHoldSpdfByLoc);
% averageFiringRatesBySpdf.postTargetDimMotorResponse = computeAverageFiringRateBySpdf(...
%         postTargetDimMotorResponseWindowOffset, targetDimT, targetDimWindow(1), targetDimSpdf, targetDimSpdfByLoc);
% averageFiringRatesBySpdf.preExitFixation = computeAverageFiringRateBySpdf(...
%         preExitFixationWindowOffset, exitFixationT, exitFixationWindow(1), exitFixationSpdf, exitFixationSpdfByLoc);
% averageFiringRatesBySpdf.postExitFixation = computeAverageFiringRateBySpdf(...
%         postExitFixationWindowOffset, exitFixationT, exitFixationWindow(1), exitFixationSpdf, exitFixationSpdfByLoc);
% 
% %% compute max firing rate across RF locations and conditions
% fn = fieldnames(averageFiringRatesByCount);
% maxFiringRateByCount = -Inf;
% minFiringRateByCount = Inf;
% for i = 1:numel(fn)
%     maxFiringRate = max(averageFiringRatesByCount.(fn{i}).byLoc);
%     minFiringRate = min(averageFiringRatesByCount.(fn{i}).byLoc);
%     if maxFiringRate > maxFiringRateByCount
%         maxFiringRateByCount = maxFiringRate;
%     end
%     if minFiringRate < minFiringRateByCount
%         minFiringRateByCount = minFiringRate;
%     end
% end
% clear fn maxFiringRate minFiringRate;
% 
% %% compute max firing rate across RF locations and conditions
% % use spdf time courses
% cueOnsetNormalizationWindowOffset = [-0.3 0.6];
% arrayOnsetHoldNormalizationWindowOffset = [-0.6 0.65];
% arrayOnsetRelNormalizationWindowOffset = [-0.6 0.25];
% targetDimNormalizationWindowOffset = [-0.65 0.25];
% 
% cueOnsetNormalizationLogical = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + cueOnsetNormalizationWindowOffset);
% arrayOnsetHoldNormalizationLogical = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + arrayOnsetHoldNormalizationWindowOffset);
% arrayOnsetRelNormalizationLogical = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + arrayOnsetRelNormalizationWindowOffset);
% targetDimNormalizationLogical = getTimeLogicalWithTolerance(targetDimT, targetDimWindow(1) + targetDimNormalizationWindowOffset);
% % exclude motor response here
% 
% % chain time courses into long vectors for each location
% allSpdfsByLoc = [cueOnsetSpdfByLoc(:,cueOnsetNormalizationLogical) ...
%         arrayOnsetSpdfByLoc(:,arrayOnsetRelNormalizationLogical) ... % be conservative and use rel here
%         arrayOnsetRelSpdfByLoc(:,arrayOnsetRelNormalizationLogical) ...
%         arrayOnsetHoldSpdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
%         targetDimSpdfByLoc(:,targetDimNormalizationLogical) ...
%         targetDimShortHoldSpdfByLoc(:,targetDimNormalizationLogical) ...
%         targetDimLongHoldSpdfByLoc(:,targetDimNormalizationLogical)];
% maxFiringRateBySpdf = max(max(allSpdfsByLoc));
% minFiringRateBySpdf = min(min(allSpdfsByLoc));
% 
% enterFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];
% cueOnsetInclMotorNormalizationWindowOffset = [-0.6 0.6];
% arrayOnsetRelInclMotorNormalizationWindowOffset = [-0.6 0.6];
% targetDimInclMotorNormalizationWindowOffset = [-0.65 0.6];
% exitFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];
% 
% enterFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(enterFixationT, enterFixationWindow(1) + enterFixationInclMotorNormalizationWindowOffset);
% cueOnsetInclMotorNormalizationLogical = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + cueOnsetInclMotorNormalizationWindowOffset);
% arrayOnsetRelInclMotorNormalizationLogical = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + arrayOnsetRelInclMotorNormalizationWindowOffset);
% targetDimInclMotorNormalizationLogical = getTimeLogicalWithTolerance(targetDimT, targetDimWindow(1) + targetDimInclMotorNormalizationWindowOffset);
% exitFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(exitFixationT, exitFixationWindow(1) + exitFixationInclMotorNormalizationWindowOffset);
% 
% allSpdfsByLoc = [enterFixationSpdfByLoc(:,enterFixationInclMotorNormalizationLogical) ...
%         cueOnsetSpdfByLoc(:,cueOnsetInclMotorNormalizationLogical) ...
%         arrayOnsetSpdfByLoc(:,arrayOnsetHoldNormalizationLogical) ... % be conservative and use hold here
%         arrayOnsetRelSpdfByLoc(:,arrayOnsetRelInclMotorNormalizationLogical) ...
%         arrayOnsetHoldSpdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
%         targetDimSpdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
%         targetDimShortHoldSpdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
%         targetDimLongHoldSpdfByLoc(:,targetDimInclMotorNormalizationLogical)...
%         exitFixationSpdfByLoc(:,exitFixationInclMotorNormalizationLogical)];
% maxFiringRateBySpdfInclMotor = max(max(allSpdfsByLoc));
% minFiringRateBySpdfInclMotor = min(min(allSpdfsByLoc));
% 
% clear allSpdfsByLoc;
% 
% %% compute RF by mean cue response
% meanCueResponseBaselineCorrByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     meanCueResponseBaselineCorrByLoc(i) = averageFiringRatesBySpdf.cueResponse.byLoc(i) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
% end
% [~,inRFLoc] = max(meanCueResponseBaselineCorrByLoc);
% assert(nLoc == 4); % next line based on nLoc == 4
% exRFLoc = mod(inRFLoc + 1, 4) + 1; % opposite location
% 
% %% compute RF by max cue response and max SD cue response
% cueResponseWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + cueResponseWindowOffset);
% maxCueResponseBaselineCorrByLoc = nan(nLoc, 1);
% maxSDCueResponseBaselineCorrByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     maxCueResponseBaselineCorrByLoc(i) = max(cueOnsetSpdfByLoc(i,cueResponseWindowIndices)) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
%     maxSDCueResponseBaselineCorrByLoc(i) = (max(cueOnsetSpdfByLoc(i,cueResponseWindowIndices)) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i)) / ...
%             averageFiringRatesBySpdf.preCueBaseline.byLocSDOverTime(i);
% end
% [~,inRFLocByMax] = max(maxCueResponseBaselineCorrByLoc);
% exRFLocByMax = mod(inRFLocByMax + 1, 4) + 1; % opposite location
% [~,inRFLocByMaxSD] = max(maxSDCueResponseBaselineCorrByLoc);
% exRFLocByMaxSD = mod(inRFLocByMaxSD + 1, 4) + 1; % opposite location
% 
% %% compute cue response at InRF > baseline
% % bootstrap on mean of baseline SPDF with 500 shuffles
% % compare actual mean SPDF cue response to distribution of bootstrapped
% % pre-cue baselines
% % TODO using count should be faster than psth
% numRandomizations = 1;
% bootstrappedMeanPreCueBaselines = zeros(numRandomizations, 1);
% preCueBaselineWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + preCueBaselineWindowOffset);
% for m = 1:numRandomizations
%     bootRandInd = randi(nTrials, 1, nTrials);
%     bootstrapPsth = fixedPsth(cueOnsetSpikeTimes(bootRandInd), kernelSigma, 0, cueOnsetT);
%     if ~isempty(bootstrapPsth)
%         bootstrappedMeanPreCueBaselines(m) = mean(bootstrapPsth(preCueBaselineWindowIndices));
%     end % else, then there were no spikes and leave as zeros
% end
% meanBootstrappedMeanPreCueBaselines = mean(bootstrappedMeanPreCueBaselines);
% 
% clear m bootRandInd bootstrapPsthResponse;
% 
% % inRFLoc is defined as location with largest mean > baseline
% % directional test
% cueResponsePValueByBootstrapBaselineSpdfInRF = sum(averageFiringRatesBySpdf.cueResponse.byLoc(inRFLoc) < bootstrappedMeanPreCueBaselines) / numRandomizations;
% 
% %% compute cue response at any location ~= baseline
% cueResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     if averageFiringRatesBySpdf.cueResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
%         cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     else
%         cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     end
% end
% 
% %% compute cue target delay at any location ~= baseline
% % TODO consider change from delay also as a response
% cueTargetDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     if averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
%         cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     else
%         cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     end
% end
% 
% %% compute array response (HOLD) at any location ~= baseline
% % TODO consider change from delay also as a response
% arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     if averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
%         arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     else
%         arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     end
% end
% 
% %% compute target dim delay at any location ~= baseline
% % TODO consider change from delay also as a response
% targetDimDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     if averageFiringRatesBySpdf.targetDimDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
%         targetDimDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     else
%         targetDimDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     end
% end
% 
% %% compute target dim response at any location ~= baseline
% % TODO consider change from delay also as a response
% targetDimResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
% for i = 1:nLoc
%     if averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
%         targetDimResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     else
%         targetDimResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
%     end
% end
% 
% %% rank-sum test on x vs baseline using individual trial spike counts
% % spike rates per trial are not normally distributed
% % they may also not have equal variance
% % use nonparametric or resampling methods
% % in fact, spike rates per trial are quantized, and because the length of
% % the time periods differ between the baseline and the other periods,
% % distributions may have different means (or maybe even medians) when that
% % is due to the low resolution of the data
% 
% % could change the analysis windows to be equal in duration -- which should
% % make the periods interchangeable, and make the rank sum test or
% % permutation tests work
% 
% % also, consider a paired test may be more powerful TODO test this
% 
% % does it matter much that the cue response by location is compared to all
% % the baseline trials together? shouldn't matter much. more data here, but 
% % lose statistical power in that it's not a paired test
% cueResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
%         averageFiringRatesByCount.cueResponse.trialRateByLoc);
% cueTargetDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
%         averageFiringRatesByCount.cueTargetDelay.trialRateByLoc);
% arrayHoldResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
%         averageFiringRatesByCount.arrayHoldResponse.trialRateByLoc);
% targetDimDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
%         averageFiringRatesByCount.targetDimDelay.trialRateByLoc);
% targetDimResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
%         averageFiringRatesByCount.targetDimResponse.trialRateByLoc);
% 
% %% shuffle test on cue onset time as another test for cue-evoked change in activity
% % concatenate the two time periods in a trial and circularly permute the
% % event time across time bins 
% % I think the event response time period has to be reasonably short (e.g.
% % doesn't have contamination from other evoked responses), but the
% % preceding time can be however long. 
% % This is similar to the permutation test of having group 1 being mean 
% % baseline activity over N1 trials and group 2 being mean cue response 
% % activity over N2 trials, pooling them and assigning N1 trials to group 1
% % and N2 trials ot group 2, and re-computing the mean activity across
% % trials P times. 
% 
% 
% %% permutation test on cue-target delay period
% numRandomizations = 1;
% shuffleCueTargetDelayDiff = zeros(numRandomizations, 1);
% shuffleCueTargetDelayAI = zeros(numRandomizations, 1);
% cueTargetDelayWindowIndices = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + cueTargetDelayWindowOffset);
% nTrialsInRF = numel(arrayOnsetSpikeTimesByLoc{inRFLoc});
% nTrialsExRF = numel(arrayOnsetSpikeTimesByLoc{exRFLoc});
% nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
% arrayOnsetSpikeTimesInRFExRF = [arrayOnsetSpikeTimesByLoc{inRFLoc} arrayOnsetSpikeTimesByLoc{exRFLoc}];
% for m = 1:numRandomizations
%     shuffleIndices = randperm(nTrialsInRFExRF);
%     randSpikeTimesInRF = arrayOnsetSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
%     randSpikeTimesExRF = arrayOnsetSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
%     
%     shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, arrayOnsetT);
%     shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, arrayOnsetT);
%     meanResponseShufflePsthInRF = mean(shufflePsthInRF(cueTargetDelayWindowIndices));
%     meanResponseShufflePsthExRF = mean(shufflePsthExRF(cueTargetDelayWindowIndices));
%     shuffleCueTargetDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
%     shuffleCueTargetDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
%             (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
% end
% cueTargetDelayDiff = averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc);
% cueTargetDelayAI = (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc)) / ...
%         (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc));
% % inRFLoc is defined as location with largest mean > baseline
% if cueTargetDelayDiff > mean(shuffleCueTargetDelayDiff)
%     cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff < shuffleCueTargetDelayDiff) / numRandomizations * 2;
% else
%     cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff > shuffleCueTargetDelayDiff) / numRandomizations * 2;
% end
% 
% if cueTargetDelayAI > mean(shuffleCueTargetDelayAI)
%     cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI < shuffleCueTargetDelayAI) / numRandomizations * 2;
% else
%     cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI > shuffleCueTargetDelayAI) / numRandomizations * 2;
% end
% 
% %% rank sum test on target-dim delay period using individual trial spike counts
% cueTargetDelayDiffPValueByRankSumTest = ranksum(...
%         averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{inRFLoc}, ...
%         averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{exRFLoc});
% 
% %% permutation test on target-dim delay period
% numRandomizations = 1;
% shuffleTargetDimDelayDiff = zeros(numRandomizations, 1);
% shuffleTargetDimDelayAI = zeros(numRandomizations, 1);
% targetDimDelayWindowIndices = getTimeLogicalWithTolerance(targetDimT, targetDimWindow(1) + targetDimDelayWindowOffset);
% nTrialsInRF = numel(targetDimSpikeTimesByLoc{inRFLoc});
% nTrialsExRF = numel(targetDimSpikeTimesByLoc{exRFLoc});
% nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
% targetDimSpikeTimesInRFExRF = [targetDimSpikeTimesByLoc{inRFLoc} targetDimSpikeTimesByLoc{exRFLoc}];
% for m = 1:numRandomizations    
%     shuffleIndices = randperm(nTrialsInRFExRF);
%     randSpikeTimesInRF = targetDimSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
%     randSpikeTimesExRF = targetDimSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
%     
%     shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, targetDimT);
%     shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, targetDimT);
%     meanResponseShufflePsthInRF = mean(shufflePsthInRF(targetDimDelayWindowIndices));
%     meanResponseShufflePsthExRF = mean(shufflePsthExRF(targetDimDelayWindowIndices));
%     shuffleTargetDimDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
%     shuffleTargetDimDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
%             (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
% end
% targetDimDelayDiff = averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc);
% targetDimDelayAI = (averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc)) / ...
%         (averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc));
% % inRFLoc is defined as location with largest mean > baseline
% if targetDimDelayDiff > mean(shuffleTargetDimDelayDiff)
%     targetDimDelayDiffPValueByShuffleSpdf = sum(targetDimDelayDiff < shuffleTargetDimDelayDiff) / numRandomizations * 2;
% else
%     targetDimDelayDiffPValueByShuffleSpdf = sum(targetDimDelayDiff > shuffleTargetDimDelayDiff) / numRandomizations * 2;
% end
% 
% if targetDimDelayAI > mean(shuffleTargetDimDelayAI)
%     targetDimDelayAIPValueByShuffleSpdf = sum(targetDimDelayAI < shuffleTargetDimDelayAI) / numRandomizations * 2;
% else
%     targetDimDelayAIPValueByShuffleSpdf = sum(targetDimDelayAI > shuffleTargetDimDelayAI) / numRandomizations * 2;
% end
% 
% %% rank-sum test on target-dim delay period using individual trial spike counts
% targetDimDelayDiffPValueByRankSumTest = ranksum(...
%         averageFiringRatesByCount.targetDimDelay.trialRateByLoc{inRFLoc}, ...
%         averageFiringRatesByCount.targetDimDelay.trialRateByLoc{exRFLoc});
% 
% %% spatial selectivity during pre-cue baseline using info rate (control)
% preCueBaselineAnalysisWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + preCueBaselineWindowOffset);
% [preCueBaselineInfoRate,shuffledPreCueBaselineInfoRates,preCueBaselineInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         cueOnsetSpikeTimes, cueOnsetSpikeTimesByLoc, averageFiringRatesBySpdf.preCueBaseline, ...
%         preCueBaselineAnalysisWindowIndices, kernelSigma, cueOnsetT);
%     
% %% spatial selectivity during cue response using info rate
% cueResponseAnalysisWindowIndices = getTimeLogicalWithTolerance(cueOnsetT, cueOnsetWindow(1) + cueResponseWindowOffset);
% [cueResponseInfoRate,shuffledCueResponseInfoRates,cueResponseInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         cueOnsetSpikeTimes, cueOnsetSpikeTimesByLoc, averageFiringRatesBySpdf.cueResponse, ...
%         cueResponseAnalysisWindowIndices, kernelSigma, cueOnsetT);
% 
% %% spatial selectivity during cue target delay period using info rate
% cueTargetDelayAnalysisWindowIndices = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + cueTargetDelayWindowOffset);
% [cueTargetDelayInfoRate,shuffledCueTargetDelayInfoRates,cueTargetDelayInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         arrayOnsetSpikeTimes, arrayOnsetSpikeTimesByLoc, averageFiringRatesBySpdf.cueTargetDelay, ...
%         cueTargetDelayAnalysisWindowIndices, kernelSigma, arrayOnsetT);
% 
% %% spatial selectivity during array response using info rate
% arrayResponseAnalysisWindowIndices = getTimeLogicalWithTolerance(arrayOnsetT, arrayOnsetWindow(1) + arrayResponseWindowOffset);
% [arrayHoldResponseInfoRate,shuffledArrayHoldResponseInfoRates,arrayHoldResponseInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         arrayOnsetHoldSpikeTimes, arrayOnsetHoldSpikeTimesByLoc, averageFiringRatesBySpdf.arrayHoldResponse, ...
%         arrayResponseAnalysisWindowIndices, kernelSigma, arrayOnsetT);
% 
% %% spatial selectivity during target dim delay period using info rate
% targetDimDelayAnalysisWindowIndices = getTimeLogicalWithTolerance(targetDimT, targetDimWindow(1) + targetDimDelayWindowOffset);
% [targetDimDelayInfoRate,shuffledTargetDimDelayInfoRates,targetDimDelayInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         targetDimSpikeTimes, targetDimSpikeTimesByLoc, averageFiringRatesBySpdf.targetDimDelay, ...
%         targetDimDelayAnalysisWindowIndices, kernelSigma, targetDimT);
% 
% %% spatial selectivity during target dimming using info rate
% targetDimResponseAnalysisWindowIndices = getTimeLogicalWithTolerance(targetDimT, targetDimWindow(1) + targetDimResponseWindowOffset);
% [targetDimResponseInfoRate,shuffledTargetDimResponseInfoRates,targetDimResponseInfoRatePValueByShuffleSpdf] = computeInfoRatePValueByShuffle(...
%         targetDimSpikeTimes, targetDimSpikeTimesByLoc, averageFiringRatesBySpdf.targetDimResponse, ...
%         targetDimResponseAnalysisWindowIndices, kernelSigma, targetDimT);
% 
% %% save
% clear i m spikeTs;
% clear nTrialsInRF nTrialsExRF nTrialsInRFExRF shuffleIndices;
% clear randSpikeTimesInRF randSpikeTimesExRF shufflePsthInRF shufflePsthExRF;
% clear meanResponseShufflePsthInRF meanResponseShufflePsthExRF;
% 
% save(saveFileName);
