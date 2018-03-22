function saveFileName = computeEvokedSpiking(saveFileName, spikeStruct, nLoc, UE)
% spikeStruct = struct of spiking data
% UE = useful events struct

spikeTs = spikeStruct.ts;
kernelSigma = 0.01;
numRandomizations = 2;

clear spikeStruct;

nTrials = numel(UE.cueOnset);
nTrialsHold = numel(UE.arrayOnsetHold);
nTrialsByLoc = cellfun(@numel, UE.cueOnsetByLoc);
isLocUsed = nTrialsByLoc > 0;
fprintf('Computing evoked spiking with SPDF sigma %0.3f seconds and %d randomizations over %d trials...\n', ...
        kernelSigma, numRandomizations, nTrials);
fprintf('\t');

%% align spikes to cue onset, compute spdf
cueOnset.window = [0.8 0.8]; % seconds before, after
cueOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
cueOnsetRel = cueOnset; % copy
cueOnsetHold = cueOnset; % copy

cueOnset = createTimeLockedSpdf(spikeTs, UE.cueOnset, UE.cueOnsetByLoc, cueOnset, kernelSigma);
cueOnsetRel = createTimeLockedSpdf(spikeTs, UE.cueOnsetRel, UE.cueOnsetRelByLoc, cueOnsetRel, kernelSigma);
cueOnsetHold = createTimeLockedSpdf(spikeTs, UE.cueOnsetHold, UE.cueOnsetHoldByLoc, cueOnsetHold, kernelSigma);

%% align spikes to array onset, compute spdf
arrayOnset.window = [0.8 0.8]; % seconds before, after
arrayOnset.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
arrayOnsetRel = arrayOnset; % copy
arrayOnsetHold = arrayOnset; % copy

arrayOnset = createTimeLockedSpdf(spikeTs, UE.arrayOnset, UE.arrayOnsetByLoc, arrayOnset, kernelSigma);
arrayOnsetRel = createTimeLockedSpdf(spikeTs, UE.arrayOnsetRel, UE.arrayOnsetRelByLoc, arrayOnsetRel, kernelSigma);
arrayOnsetHold = createTimeLockedSpdf(spikeTs, UE.arrayOnsetHold, UE.arrayOnsetHoldByLoc, arrayOnsetHold, kernelSigma);

%% align spikes to target dimming, compute spdf
targetDim.window = [0.8 0.8]; % seconds before, after
targetDim.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
targetDimShortHold = targetDim;
targetDimLongHold = targetDim;

targetDim = createTimeLockedSpdf(spikeTs, UE.targetDim, UE.targetDimByLoc, targetDim, kernelSigma);
targetDimShortHold = createTimeLockedSpdf(spikeTs, UE.targetDimShortHold, UE.targetDimShortHoldByLoc, targetDimShortHold, kernelSigma);
targetDimLongHold = createTimeLockedSpdf(spikeTs, UE.targetDimLongHold, UE.targetDimLongHoldByLoc, targetDimLongHold, kernelSigma);

%% look at lever-locked responses and saccade-locked responses
% note that this includes both hold and release trials
assert(~isempty(UE.fixationAndLeverTimes))

%% align spikes to enter/exit fixation, compute spdf
enterFixation.window = [0.8 0.8]; % seconds before, after
enterFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
enterFixationRel = enterFixation; % copy
enterFixationHold = enterFixation; % copy

enterFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCueByLoc, enterFixation, kernelSigma);
enterFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationRelTimesPreCueByLoc, enterFixationRel, kernelSigma);
enterFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCue, UE.fixationAndLeverTimes.firstEnterFixationHoldTimesPreCueByLoc, enterFixationHold, kernelSigma);

exitFixation.window = [0.8 0.8]; % seconds before, after
exitFixation.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects
exitFixationLeft = exitFixation; % copy
exitFixationRight = exitFixation; % copy
exitFixationRel = exitFixation; % copy
exitFixationHold = exitFixation; % copy

exitFixation = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuiceByLoc, exitFixation, kernelSigma);
exitFixationLeft = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationLeftTimesAroundJuiceByLoc, exitFixationLeft, kernelSigma);
exitFixationRight = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationRightTimesAroundJuiceByLoc, exitFixationRight, kernelSigma);
exitFixationRel = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationRelTimesAroundJuiceByLoc, exitFixationRel, kernelSigma);
exitFixationHold = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuice, UE.fixationAndLeverTimes.firstExitFixationHoldTimesAroundJuiceByLoc, exitFixationHold, kernelSigma);

%% align spikes to lever press and release, compute spdf
leverPress.window = [0.8 0.8]; % seconds before, after
leverPress.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverPress = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverPressTimesPreCue, UE.fixationAndLeverTimes.firstLeverPressTimesPreCueByLoc, leverPress, kernelSigma);

leverRelease.window = [0.8 0.8]; % seconds before, after
leverRelease.spdfWindowOffset = [-0.7 0.7]; % tighter window for spdf to avoid edge effects

leverRelease = createTimeLockedSpdf(spikeTs, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuiceByLoc, leverRelease, kernelSigma);

%% calc time between motor events
assert(numel(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue) == numel(UE.fixationAndLeverTimes.firstLeverPressTimesPreCue));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice));

arrayOnsetRelToJuiceEventTime = UE.firstJuiceEvent(~UE.isHoldTrial) - UE.arrayOnsetRel; % almost same as UE.rt(~isHoldTrial)?
targetDimToJuiceEventTime = UE.firstJuiceEvent(UE.isHoldTrial) - UE.targetDim; % almost same as UE.rt(isHoldTrial)?
enterFixationToLeverPressTime = UE.fixationAndLeverTimes.firstLeverPressTimesPreCue - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
enterFixationToCueOnsetTime = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
exitFixationToJuiceEventTime = UE.firstJuiceEvent - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
exitFixationToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
arrayOnsetRelToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;
arrayOnsetRelToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
targetDimToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;

%% analysis time windows
preLeverReleaseWindowOffset = [-0.175 0];
preCueBaselineWindowOffset = [-0.175 0];
cueResponseWindowOffset = [0.025 0.2];
cueTargetDelayWindowOffset = [-0.175 0];
cueTargetDelayLongWindowOffset = [-0.4 0];
arrayResponseWindowOffset = [0.025 0.2];
targetDimDelayWindowOffset = [-0.175 0];
targetDimDelayLongWindowOffset = [-0.4 0];
targetDimResponseWindowOffset = [0.025 0.2]; 
postTargetDimMotorResponseWindowOffset = [0.425 0.6];
preEnterFixationWindowOffset = [-0.175 0];
postEnterFixationWindowOffset = [0.025 0.2];
postEnterFixationLateWindowOffset = [0.15 0.325];
preExitFixationWindowOffset = [-0.05 0]; % note smaller window - TODO does this matter much?
preExitFixationEarlyWindowOffset = [-0.2 -0.05]; % note smaller window - TODO does this matter much?
postExitFixationWindowOffset = [0.0 0.1]; % note smaller window - TODO does this matter much?
postLeverReleaseWindowOffset = [0.025 0.2];

%% average number of spikes in each window
averageFiringRatesByCount = struct();

averageFiringRatesByCount.preEnterFixation = computeAverageFiringRateByCount(...
        preEnterFixationWindowOffset, enterFixation);
averageFiringRatesByCount.postEnterFixation = computeAverageFiringRateByCount(...
        postEnterFixationWindowOffset, enterFixation);
averageFiringRatesByCount.postEnterFixationLate = computeAverageFiringRateByCount(...
        postEnterFixationLateWindowOffset, enterFixation);
averageFiringRatesByCount.preCueBaseline = computeAverageFiringRateByCount(...
        preCueBaselineWindowOffset, cueOnset);
averageFiringRatesByCount.cueResponse = computeAverageFiringRateByCount(...
        cueResponseWindowOffset, cueOnset);
averageFiringRatesByCount.cueTargetDelay = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnset);
averageFiringRatesByCount.cueTargetDelayRel = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetRel);
averageFiringRatesByCount.cueTargetDelayHold = computeAverageFiringRateByCount(...
        cueTargetDelayWindowOffset, arrayOnsetHold);
averageFiringRatesByCount.cueTargetDelayLong = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnset);
averageFiringRatesByCount.cueTargetDelayLongRel = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetRel);
averageFiringRatesByCount.cueTargetDelayLongHold = computeAverageFiringRateByCount(...
        cueTargetDelayLongWindowOffset, arrayOnsetHold);
averageFiringRatesByCount.arrayResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnset);
averageFiringRatesByCount.arrayRelResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetRel);
averageFiringRatesByCount.arrayHoldResponse = computeAverageFiringRateByCount(...
        arrayResponseWindowOffset, arrayOnsetHold);
averageFiringRatesByCount.targetDimDelay = computeAverageFiringRateByCount(...
        targetDimDelayWindowOffset, targetDim);
averageFiringRatesByCount.targetDimDelayLong = computeAverageFiringRateByCount(...
        targetDimDelayLongWindowOffset, targetDim);
averageFiringRatesByCount.targetDimResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDim);
averageFiringRatesByCount.targetDimShortHoldResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimShortHold);
averageFiringRatesByCount.targetDimLongHoldResponse = computeAverageFiringRateByCount(...
        targetDimResponseWindowOffset, targetDimLongHold);
averageFiringRatesByCount.postTargetDimMotorResponse = computeAverageFiringRateByCount(...
        postTargetDimMotorResponseWindowOffset, targetDim);
averageFiringRatesByCount.preExitFixation = computeAverageFiringRateByCount(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesByCount.preExitFixationEarly = computeAverageFiringRateByCount(...
        preExitFixationEarlyWindowOffset, exitFixation);
averageFiringRatesByCount.postExitFixation = computeAverageFiringRateByCount(...
        postExitFixationWindowOffset, exitFixation);

%% average spdf value in each window
averageFiringRatesBySpdf = struct();
averageFiringRatesBySpdf.kernelSigma = kernelSigma;
averageFiringRatesBySpdf.preEnterFixation = computeAverageFiringRateBySpdf(...
        preEnterFixationWindowOffset, enterFixation);
averageFiringRatesBySpdf.postEnterFixation = computeAverageFiringRateBySpdf(...
        postEnterFixationWindowOffset, enterFixation);
averageFiringRatesBySpdf.postEnterFixationLate = computeAverageFiringRateBySpdf(...
        postEnterFixationLateWindowOffset, enterFixation);
averageFiringRatesBySpdf.preCueBaseline = computeAverageFiringRateBySpdf(...
        preCueBaselineWindowOffset, cueOnset);
averageFiringRatesBySpdf.cueResponse = computeAverageFiringRateBySpdf(...
        cueResponseWindowOffset, cueOnset);
averageFiringRatesBySpdf.cueTargetDelay = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnset);
averageFiringRatesBySpdf.cueTargetDelayRel = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetRel);
averageFiringRatesBySpdf.cueTargetDelayHold = computeAverageFiringRateBySpdf(...
        cueTargetDelayWindowOffset, arrayOnsetHold);
averageFiringRatesBySpdf.cueTargetDelayLong = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnset);
averageFiringRatesBySpdf.cueTargetDelayLongRel = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetRel);
averageFiringRatesBySpdf.cueTargetDelayLongHold = computeAverageFiringRateBySpdf(...
        cueTargetDelayLongWindowOffset, arrayOnsetHold);
averageFiringRatesBySpdf.arrayResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnset);
averageFiringRatesBySpdf.arrayRelResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetRel);
averageFiringRatesBySpdf.arrayHoldResponse = computeAverageFiringRateBySpdf(...
        arrayResponseWindowOffset, arrayOnsetHold);
averageFiringRatesBySpdf.targetDimDelay = computeAverageFiringRateBySpdf(...
        targetDimDelayWindowOffset, targetDim);
averageFiringRatesBySpdf.targetDimDelayLong = computeAverageFiringRateBySpdf(...
        targetDimDelayLongWindowOffset, targetDim);
averageFiringRatesBySpdf.targetDimResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDim);
averageFiringRatesBySpdf.targetDimShortHoldResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimShortHold);
averageFiringRatesBySpdf.targetDimLongHoldResponse = computeAverageFiringRateBySpdf(...
        targetDimResponseWindowOffset, targetDimLongHold);
averageFiringRatesBySpdf.postTargetDimMotorResponse = computeAverageFiringRateBySpdf(...
        postTargetDimMotorResponseWindowOffset, targetDim);
averageFiringRatesBySpdf.preExitFixation = computeAverageFiringRateBySpdf(...
        preExitFixationWindowOffset, exitFixation);
averageFiringRatesBySpdf.preExitFixationEarly = computeAverageFiringRateBySpdf(...
        preExitFixationEarlyWindowOffset, exitFixation);
averageFiringRatesBySpdf.postExitFixation = computeAverageFiringRateBySpdf(...
        postExitFixationWindowOffset, exitFixation);

%% compute max firing rate across RF locations and conditions
fn = fieldnames(averageFiringRatesByCount);
maxFiringRateByCount = -Inf;
minFiringRateByCount = Inf;
for i = 1:numel(fn)
    maxFiringRate = max(averageFiringRatesByCount.(fn{i}).byLoc);
    minFiringRate = min(averageFiringRatesByCount.(fn{i}).byLoc);
    if maxFiringRate > maxFiringRateByCount
        maxFiringRateByCount = maxFiringRate;
    end
    if minFiringRate < minFiringRateByCount
        minFiringRateByCount = minFiringRate;
    end
end
clear fn maxFiringRate minFiringRate;

%% compute max firing rate across RF locations and conditions
% use spdf time courses
cueOnsetNormalizationWindowOffset = [-0.3 0.6];
arrayOnsetHoldNormalizationWindowOffset = [-0.6 0.65];
arrayOnsetRelNormalizationWindowOffset = [-0.6 0.25];
targetDimNormalizationWindowOffset = [-0.65 0.25];

cueOnsetNormalizationLogical = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueOnsetNormalizationWindowOffset);
arrayOnsetHoldNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetHoldNormalizationWindowOffset);
arrayOnsetRelNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetRelNormalizationWindowOffset);
targetDimNormalizationLogical = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimNormalizationWindowOffset);
% exclude motor response here

% chain time courses into long vectors for each location
allSpdfsByLoc = [cueOnset.spdfByLoc(:,cueOnsetNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ... % be conservative and use rel here
        arrayOnsetRel.spdfByLoc(:,arrayOnsetRelNormalizationLogical) ...
        arrayOnsetHold.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDim.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimShortHold.spdfByLoc(:,targetDimNormalizationLogical) ...
        targetDimLongHold.spdfByLoc(:,targetDimNormalizationLogical)];
maxFiringRateBySpdf = max(max(allSpdfsByLoc));
minFiringRateBySpdf = min(min(allSpdfsByLoc));

enterFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];
cueOnsetInclMotorNormalizationWindowOffset = [-0.6 0.6];
arrayOnsetRelInclMotorNormalizationWindowOffset = [-0.6 0.6];
targetDimInclMotorNormalizationWindowOffset = [-0.65 0.6];
exitFixationInclMotorNormalizationWindowOffset = [-0.4 0.4];

enterFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(enterFixation.t, enterFixation.window(1) + enterFixationInclMotorNormalizationWindowOffset);
cueOnsetInclMotorNormalizationLogical = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueOnsetInclMotorNormalizationWindowOffset);
arrayOnsetRelInclMotorNormalizationLogical = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + arrayOnsetRelInclMotorNormalizationWindowOffset);
targetDimInclMotorNormalizationLogical = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimInclMotorNormalizationWindowOffset);
exitFixationInclMotorNormalizationLogical = getTimeLogicalWithTolerance(exitFixation.t, exitFixation.window(1) + exitFixationInclMotorNormalizationWindowOffset);

allSpdfsByLoc = [enterFixation.spdfByLoc(:,enterFixationInclMotorNormalizationLogical) ...
        cueOnset.spdfByLoc(:,cueOnsetInclMotorNormalizationLogical) ...
        arrayOnset.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ... % be conservative and use hold here
        arrayOnsetRel.spdfByLoc(:,arrayOnsetRelInclMotorNormalizationLogical) ...
        arrayOnsetHold.spdfByLoc(:,arrayOnsetHoldNormalizationLogical) ...
        targetDim.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimShortHold.spdfByLoc(:,targetDimInclMotorNormalizationLogical) ...
        targetDimLongHold.spdfByLoc(:,targetDimInclMotorNormalizationLogical)...
        exitFixation.spdfByLoc(:,exitFixationInclMotorNormalizationLogical)];
maxFiringRateBySpdfInclMotor = max(max(allSpdfsByLoc));
minFiringRateBySpdfInclMotor = min(min(allSpdfsByLoc));

clear allSpdfsByLoc;
clear cueOnsetNormalizationLogical arrayOnsetHoldNormalizationLogical ...
        arrayOnsetRelNormalizationLogical targetDimNormalizationLogical ...
        enterFixationInclMotorNormalizationLogical cueOnsetInclMotorNormalizationLogical ...
        arrayOnsetRelInclMotorNormalizationLogical targetDimInclMotorNormalizationLogical ...
        exitFixationInclMotorNormalizationLogical;

%% compute RF as largest mean cue response over baseline
meanCueResponseBaselineCorrByLoc = nan(nLoc, 1);
for i = 1:nLoc
    meanCueResponseBaselineCorrByLoc(i) = averageFiringRatesBySpdf.cueResponse.byLoc(i) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
end
[~,inRFLoc] = max(meanCueResponseBaselineCorrByLoc);
assert(nLoc == 4); % next line based on nLoc == 4
exRFLoc = mod(inRFLoc + 1, 4) + 1; % opposite location

%% compute RF by max cue response and max SD cue response
cueResponseWindowIndices = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + cueResponseWindowOffset);
maxCueResponseBaselineCorrByLoc = nan(nLoc, 1);
maxSDCueResponseBaselineCorrByLoc = nan(nLoc, 1);
for i = 1:nLoc
    maxCueResponseBaselineCorrByLoc(i) = max(cueOnset.spdfByLoc(i,cueResponseWindowIndices)) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i);
    maxSDCueResponseBaselineCorrByLoc(i) = (max(cueOnset.spdfByLoc(i,cueResponseWindowIndices)) - averageFiringRatesBySpdf.preCueBaseline.byLoc(i)) / ...
            averageFiringRatesBySpdf.preCueBaseline.byLocSDOverTime(i);
end
[~,inRFLocByMax] = max(maxCueResponseBaselineCorrByLoc);
exRFLocByMax = mod(inRFLocByMax + 1, 4) + 1; % opposite location
[~,inRFLocByMaxSD] = max(maxSDCueResponseBaselineCorrByLoc);
exRFLocByMaxSD = mod(inRFLocByMaxSD + 1, 4) + 1; % opposite location

clear cueResponseWindowIndices maxCueResponseBaselineCorrByLoc maxSDCueResponseBaselineCorrByLoc;

%% compute cue response at InRF > baseline
% bootstrap on mean of baseline SPDF with 500 shuffles
% compare actual mean SPDF cue response to distribution of bootstrapped
% pre-cue baselines
% TODO using count should be faster than psth
bootstrappedMeanPreCueBaselines = zeros(numRandomizations, 1);
preCueBaselineWindowIndices = getTimeLogicalWithTolerance(cueOnset.t, cueOnset.window(1) + preCueBaselineWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrials, 1, nTrials);
    bootstrapPsth = fixedPsth(cueOnset.spikeTimes(bootRandInd), kernelSigma, 0, cueOnset.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanPreCueBaselines(m) = mean(bootstrapPsth(preCueBaselineWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanPreCueBaselines = mean(bootstrappedMeanPreCueBaselines);
clear m bootRandInd bootstrapPsthResponse preCueBaselineWindowIndices;

% inRFLoc is defined as location with largest mean > baseline
% directional test
cueResponsePValueByBootstrapBaselineSpdfInRF = sum(averageFiringRatesBySpdf.cueResponse.byLoc(inRFLoc) < bootstrappedMeanPreCueBaselines) / numRandomizations;

%% compute cue response at any location ~= baseline
cueResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.cueResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        cueResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute cue target delay at any location ~= baseline
% TODO consider change from delay also as a response
cueTargetDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.cueTargetDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute array response (HOLD) at any location ~= baseline
% TODO consider change from delay also as a response
arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute target dim delay at any location ~= baseline
% TODO consider change from delay also as a response
targetDimDelayPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimDelay.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        targetDimDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimDelay.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        targetDimDelayPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimDelay.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute target dim response at any location ~= baseline
% TODO consider change from delay also as a response
targetDimResponsePValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        targetDimResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        targetDimResponsePValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute pre exit fixation response at any location ~= baseline
% TODO consider change from delay also as a response
preExitFixationPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.preExitFixation.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        preExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixation.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        preExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixation.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute post exit fixation response at any location ~= baseline
% TODO consider change from delay also as a response
postExitFixationPValueByBootstrapBaselineSpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.postExitFixation.byLoc(i) > meanBootstrappedMeanPreCueBaselines
        postExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.postExitFixation.byLoc(i) < bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    else
        postExitFixationPValueByBootstrapBaselineSpdfByLoc(i) = sum(averageFiringRatesBySpdf.postExitFixation.byLoc(i) > bootstrappedMeanPreCueBaselines) / numRandomizations * 2;
    end
end

%% compute array response at InRF > cue-target delay
% bootstrap on mean of cue-target delay SPDF with 500 shuffles
% compare actual mean SPDF array response to distribution of bootstrapped
% cue-target delay
% TODO using count should be faster than psth
bootstrappedMeanCueTargetDelays = zeros(numRandomizations, 1);
cueTargetDelayWindowIndices = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + cueTargetDelayWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrialsHold, 1, nTrialsHold);
    bootstrapPsth = fixedPsth(arrayOnset.spikeTimes(bootRandInd), kernelSigma, 0, arrayOnset.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanCueTargetDelays(m) = mean(bootstrapPsth(cueTargetDelayWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanCueTargetDelays = mean(bootstrappedMeanCueTargetDelays);
clear m bootRandInd bootstrapPsthResponse cueTargetDelayWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfInRF = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(inRFLoc) < bootstrappedMeanCueTargetDelays) / numRandomizations;

%% compute array response (HOLD) at any location ~= cue-target delay
arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > meanBootstrappedMeanCueTargetDelays
        arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) < bootstrappedMeanCueTargetDelays) / numRandomizations * 2;
    else
        arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.arrayHoldResponse.byLoc(i) > bootstrappedMeanCueTargetDelays) / numRandomizations * 2;
    end
end

%% compute target-dim response at InRF > target-dim delay
% bootstrap on mean of target-dim delay SPDF with 500 shuffles
% compare actual mean SPDF target-dim response to distribution of bootstrapped
% target-dim delay
% TODO using count should be faster than psth
bootstrappedMeanTargetDimDelays = zeros(numRandomizations, 1);
targetDimDelayWindowIndices = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimDelayWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrialsHold, 1, nTrialsHold);
    bootstrapPsth = fixedPsth(targetDim.spikeTimes(bootRandInd), kernelSigma, 0, targetDim.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanTargetDimDelays(m) = mean(bootstrapPsth(targetDimDelayWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanTargetDimDelays = mean(bootstrappedMeanTargetDimDelays);
clear m bootRandInd bootstrapPsthResponse targetDimDelayWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
targetDimResponsePValueByBootstrapTargetDimDelaySpdfInRF = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(inRFLoc) < bootstrappedMeanTargetDimDelays) / numRandomizations;

%% compute target dim response at any location ~= target-dim delay
targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > meanBootstrappedMeanTargetDimDelays
        targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) < bootstrappedMeanTargetDimDelays) / numRandomizations * 2;
    else
        targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > bootstrappedMeanTargetDimDelays) / numRandomizations * 2;
    end
end

%% compute target-dim response at InRF > target-dim delay
% bootstrap on mean of target-dim delay SPDF with 500 shuffles
% compare actual mean SPDF target-dim response to distribution of bootstrapped
% target-dim delay
% TODO using count should be faster than psth
bootstrappedMeanTargetDimDelays = zeros(numRandomizations, 1);
targetDimDelayWindowIndices = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimDelayWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrialsHold, 1, nTrialsHold);
    bootstrapPsth = fixedPsth(targetDim.spikeTimes(bootRandInd), kernelSigma, 0, targetDim.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanTargetDimDelays(m) = mean(bootstrapPsth(targetDimDelayWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanTargetDimDelays = mean(bootstrappedMeanTargetDimDelays);
clear m bootRandInd bootstrapPsthResponse targetDimDelayWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
targetDimResponsePValueByBootstrapTargetDimDelaySpdfInRF = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(inRFLoc) < bootstrappedMeanTargetDimDelays) / numRandomizations;

%% compute target dim response at any location ~= target-dim delay
targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > meanBootstrappedMeanTargetDimDelays
        targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) < bootstrappedMeanTargetDimDelays) / numRandomizations * 2;
    else
        targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc(i) = sum(averageFiringRatesBySpdf.targetDimResponse.byLoc(i) > bootstrappedMeanTargetDimDelays) / numRandomizations * 2;
    end
end

%% compute pre-exit fixation response at InRF > earlier pre-exit fixation window
% bootstrap on mean of pre-exit fixation early SPDF with 500 shuffles
% compare actual mean SPDF pre-exit fixation response to distribution of bootstrapped
% earlier pre-exit fixation responses
% TODO using count should be faster than psth
bootstrappedMeanPreExitFixationEarlys = zeros(numRandomizations, 1);
preExitFixationEarlyWindowIndices = getTimeLogicalWithTolerance(exitFixation.t, exitFixation.window(1) + preExitFixationEarlyWindowOffset);
for m = 1:numRandomizations
    bootRandInd = randi(nTrials, 1, nTrials);
    bootstrapPsth = fixedPsth(exitFixation.spikeTimes(bootRandInd), kernelSigma, 0, exitFixation.t);
    if ~isempty(bootstrapPsth)
        bootstrappedMeanPreExitFixationEarlys(m) = mean(bootstrapPsth(preExitFixationEarlyWindowIndices));
    end % else, then there were no spikes and leave as zeros
end
fprintf('.');
meanBootstrappedMeanPreExitFixationEarlys = mean(bootstrappedMeanPreExitFixationEarlys);
clear m bootRandInd bootstrapPsthResponse preExitFixationEarlyWindowIndices;

% inRFLoc is defined as location with largest mean cue response > baseline
% directional test
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfInRF = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(inRFLoc) < bootstrappedMeanPreExitFixationEarlys) / numRandomizations;

%% compute pre-exit fixation response at any location ~= earlier pre-exit fixation window
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc = nan(nLoc, 1);
for i = 1:nLoc
    if averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) > meanBootstrappedMeanPreExitFixationEarlys
        preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) < bootstrappedMeanPreExitFixationEarlys) / numRandomizations * 2;
    else
        preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(i) = sum(averageFiringRatesBySpdf.preExitFixationEarly.byLoc(i) > bootstrappedMeanPreExitFixationEarlys) / numRandomizations * 2;
    end
end

%% rank-sum test on x vs baseline using individual trial spike counts
% spike rates per trial are not normally distributed
% they may also not have equal variance
% use nonparametric or resampling methods
% in fact, spike rates per trial are quantized, and because the length of
% the time periods differ between the baseline and the other periods,
% distributions may have different means (or maybe even medians) when that
% is due to the low resolution of the data

% could change the analysis windows to be equal in duration -- which should
% make the periods interchangeable, and make the rank sum test or
% permutation tests work

% also, consider a paired test may be more powerful TODO test this

% does it matter much that the cue response by location is compared to all
% the baseline trials together? shouldn't matter much. more data here, but 
% lose statistical power in that it's not a paired test
cueResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueResponse.trialRateByLoc);
cueTargetDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc);
arrayHoldResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.arrayHoldResponse.trialRateByLoc);
targetDimDelayVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimDelay.trialRateByLoc);
targetDimResponseVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.targetDimResponse.trialRateByLoc);
preExitFixationVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);
postExitFixationVsBaselineRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preCueBaseline.trialRate, ...
        averageFiringRatesByCount.postExitFixation.trialRateByLoc);

arrayHoldResponseVsCueTargetDelayRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.cueTargetDelay.trialRate, ...
        averageFiringRatesByCount.arrayHoldResponse.trialRateByLoc);
targetDimResponseVsTargetDimDelayRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.targetDimDelay.trialRate, ...
        averageFiringRatesByCount.targetDimResponse.trialRateByLoc);
preExitFixationVsPreExitFixationEarlyRankSumTestStatsByLoc = computeRankSumTestByLoc(averageFiringRatesByCount.preExitFixationEarly.trialRate, ...
        averageFiringRatesByCount.preExitFixation.trialRateByLoc);

%% shuffle test on cue onset time as another test for cue-evoked change in activity
% concatenate the two time periods in a trial and circularly permute the
% event time across time bins 
% I think the event response time period has to be reasonably short (e.g.
% doesn't have contamination from other evoked responses), but the
% preceding time can be however long. 
% This is similar to the permutation test of having group 1 being mean 
% baseline activity over N1 trials and group 2 being mean cue response 
% activity over N2 trials, pooling them and assigning N1 trials to group 1
% and N2 trials ot group 2, and re-computing the mean activity across
% trials P times. 

%% permutation test on cue-target delay period
shuffleCueTargetDelayDiff = zeros(numRandomizations, 1);
shuffleCueTargetDelayAI = zeros(numRandomizations, 1);
cueTargetDelayWindowIndices = getTimeLogicalWithTolerance(arrayOnset.t, arrayOnset.window(1) + cueTargetDelayWindowOffset);
nTrialsInRF = numel(arrayOnset.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(arrayOnset.spikeTimesByLoc{exRFLoc});
nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
arrayOnsetSpikeTimesInRFExRF = [arrayOnset.spikeTimesByLoc{inRFLoc} arrayOnset.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations
    shuffleIndices = randperm(nTrialsInRFExRF);
    randSpikeTimesInRF = arrayOnsetSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = arrayOnsetSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, arrayOnset.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, arrayOnset.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(cueTargetDelayWindowIndices));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(cueTargetDelayWindowIndices));
    shuffleCueTargetDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleCueTargetDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end
fprintf('.');
cueTargetDelayDiff = averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc);
cueTargetDelayAI = (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc)) / ...
        (averageFiringRatesBySpdf.cueTargetDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.cueTargetDelay.byLoc(exRFLoc));
% inRFLoc is defined as location with largest mean > baseline
if cueTargetDelayDiff > mean(shuffleCueTargetDelayDiff)
    cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff < shuffleCueTargetDelayDiff) / numRandomizations * 2;
else
    cueTargetDelayDiffPValueByShuffleSpdf = sum(cueTargetDelayDiff > shuffleCueTargetDelayDiff) / numRandomizations * 2;
end

if cueTargetDelayAI > mean(shuffleCueTargetDelayAI)
    cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI < shuffleCueTargetDelayAI) / numRandomizations * 2;
else
    cueTargetDelayAIPValueByShuffleSpdf = sum(cueTargetDelayAI > shuffleCueTargetDelayAI) / numRandomizations * 2;
end
clear cueTargetDelayWindowIndices arrayOnsetSpikeTimesInRFExRF;

%% rank sum test on target-dim delay period using individual trial spike counts
cueTargetDelayDiffPValueByRankSumTest = ranksum(...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.cueTargetDelay.trialRateByLoc{exRFLoc});

%% permutation test on target-dim delay period
shuffleTargetDimDelayDiff = zeros(numRandomizations, 1);
shuffleTargetDimDelayAI = zeros(numRandomizations, 1);
targetDimDelayWindowIndices = getTimeLogicalWithTolerance(targetDim.t, targetDim.window(1) + targetDimDelayWindowOffset);
nTrialsInRF = numel(targetDim.spikeTimesByLoc{inRFLoc});
nTrialsExRF = numel(targetDim.spikeTimesByLoc{exRFLoc});
nTrialsInRFExRF = nTrialsInRF + nTrialsExRF;
targetDimSpikeTimesInRFExRF = [targetDim.spikeTimesByLoc{inRFLoc} targetDim.spikeTimesByLoc{exRFLoc}];
for m = 1:numRandomizations    
    shuffleIndices = randperm(nTrialsInRFExRF);
    randSpikeTimesInRF = targetDimSpikeTimesInRFExRF(shuffleIndices(1:nTrialsInRF));
    randSpikeTimesExRF = targetDimSpikeTimesInRFExRF(shuffleIndices(nTrialsInRF+1:end));
    
    shufflePsthInRF = fixedPsth(randSpikeTimesInRF, kernelSigma, 0, targetDim.t);
    shufflePsthExRF = fixedPsth(randSpikeTimesExRF, kernelSigma, 0, targetDim.t);
    meanResponseShufflePsthInRF = mean(shufflePsthInRF(targetDimDelayWindowIndices));
    meanResponseShufflePsthExRF = mean(shufflePsthExRF(targetDimDelayWindowIndices));
    shuffleTargetDimDelayDiff(m) = meanResponseShufflePsthInRF - meanResponseShufflePsthExRF;
    shuffleTargetDimDelayAI(m) = (meanResponseShufflePsthInRF - meanResponseShufflePsthExRF) / ...
            (meanResponseShufflePsthInRF + meanResponseShufflePsthExRF);
end
fprintf('.');
targetDimDelayDiff = averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc);
targetDimDelayAI = (averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) - averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc)) / ...
        (averageFiringRatesBySpdf.targetDimDelay.byLoc(inRFLoc) + averageFiringRatesBySpdf.targetDimDelay.byLoc(exRFLoc));
% inRFLoc is defined as location with largest mean > baseline
if targetDimDelayDiff > mean(shuffleTargetDimDelayDiff)
    targetDimDelayDiffPValueByShuffleSpdf = sum(targetDimDelayDiff < shuffleTargetDimDelayDiff) / numRandomizations * 2;
else
    targetDimDelayDiffPValueByShuffleSpdf = sum(targetDimDelayDiff > shuffleTargetDimDelayDiff) / numRandomizations * 2;
end

if targetDimDelayAI > mean(shuffleTargetDimDelayAI)
    targetDimDelayAIPValueByShuffleSpdf = sum(targetDimDelayAI < shuffleTargetDimDelayAI) / numRandomizations * 2;
else
    targetDimDelayAIPValueByShuffleSpdf = sum(targetDimDelayAI > shuffleTargetDimDelayAI) / numRandomizations * 2;
end
clear targetDimDelayWindowIndices targetDimSpikeTimesInRFExRF;

%% rank-sum test on target-dim delay period using individual trial spike counts
targetDimDelayDiffPValueByRankSumTest = ranksum(...
        averageFiringRatesByCount.targetDimDelay.trialRateByLoc{inRFLoc}, ...
        averageFiringRatesByCount.targetDimDelay.trialRateByLoc{exRFLoc});

%% spatial selectivity during pre-cue baseline using info rate (control)
preCueBaselineInfoRateStruct = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.preCueBaseline, preCueBaselineWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during cue response using info rate
cueResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        cueOnset, averageFiringRatesBySpdf.cueResponse, cueResponseWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during cue target delay period using info rate
cueTargetDelayInfoRateStruct = computeInfoRatePValueByShuffle(...
        arrayOnset, averageFiringRatesBySpdf.cueTargetDelay, cueTargetDelayWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during array response using info rate
arrayHoldResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        arrayOnsetHold, averageFiringRatesBySpdf.arrayHoldResponse, arrayResponseWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during target dim delay period using info rate
targetDimDelayInfoRateStruct = computeInfoRatePValueByShuffle(...
        targetDim, averageFiringRatesBySpdf.targetDimDelay, targetDimDelayWindowOffset, numRandomizations);
fprintf('.');

%% spatial selectivity during target dimming using info rate
targetDimResponseInfoRateStruct = computeInfoRatePValueByShuffle(...
        targetDim, averageFiringRatesBySpdf.targetDimResponse, targetDimResponseWindowOffset, numRandomizations);
fprintf('.');

%% latency of cue response by location
cueOnset = computeResponseLatencyByLoc(cueOnset, isLocUsed);
arrayOnset = computeResponseLatencyByLoc(arrayOnset, isLocUsed);
arrayOnsetRel = computeResponseLatencyByLoc(arrayOnsetRel, isLocUsed);
arrayOnsetHold = computeResponseLatencyByLoc(arrayOnsetHold, isLocUsed);
targetDim = computeResponseLatencyByLoc(targetDim, isLocUsed);

%% clean up unused locations
% temp until the change is made in computeEvokedSpiking
cueResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
cueTargetDelayPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
arrayHoldResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
targetDimDelayPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
targetDimResponsePValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
preExitFixationPValueByBootstrapBaselineSpdfByLoc(~isLocUsed) = NaN;
arrayHoldResponsePValueByBootstrapCueTargetDelaySpdfByLoc(~isLocUsed) = NaN;
targetDimResponsePValueByBootstrapTargetDimDelaySpdfByLoc(~isLocUsed) = NaN;
preExitFixationPValueByBootstrapPreExitFixationEarlySpdfByLoc(~isLocUsed) = NaN;

unusedLocs = find(~isLocUsed);
for i = 1:numel(unusedLocs)
    k = unusedLocs(i);
    if k <= nLoc
        cueResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        cueTargetDelayVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        arrayHoldResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        targetDimDelayVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        targetDimResponseVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        preExitFixationVsBaselineRankSumTestStatsByLoc(k).p = NaN;
        arrayHoldResponseVsCueTargetDelayRankSumTestStatsByLoc(k).p = NaN;
        targetDimResponseVsTargetDimDelayRankSumTestStatsByLoc(k).p = NaN;
        preExitFixationVsPreExitFixationEarlyRankSumTestStatsByLoc(k).p = NaN;
    end
end

%% save
clear i m k;
clear nTrialsInRF nTrialsExRF nTrialsInRFExRF shuffleIndices;
clear randSpikeTimesInRF randSpikeTimesExRF shufflePsthInRF shufflePsthExRF;
clear meanResponseShufflePsthInRF meanResponseShufflePsthExRF;

% save all workspace variables not cleared above into mat file
fprintf('\n');
fprintf('Saving to %s.\n', saveFileName);
save(saveFileName);

end

function timeLockStruct = createTimeLockedSpdf(spikeTs, eventTimes, eventTimesByLoc, timeLockStruct, kernelSigma)
% timeLockStruct.window and timeLockStruct.spdfWindowOffset must exist
nLoc = numel(eventTimesByLoc);

timeLockStruct.kernelSigma = kernelSigma;
timeLockStruct.t = computeTForSpdf(timeLockStruct.window(1), timeLockStruct.spdfWindowOffset, kernelSigma);

timeLockStruct.spikeTimes = createnonemptydatamatpt(spikeTs, eventTimes, timeLockStruct.window);

[timeLockStruct.spdf,~,timeLockStruct.spdfErr] = fixedPsth(timeLockStruct.spikeTimes, kernelSigma, 2, timeLockStruct.t);

timeLockStruct.spikeTimesByLoc = cell(nLoc, 1);
timeLockStruct.spdfByLoc = nan(nLoc, numel(timeLockStruct.t));
timeLockStruct.spdfErrByLoc = nan(nLoc, numel(timeLockStruct.t));
for i = 1:nLoc
    timeLockStruct.spikeTimesByLoc{i} = createnonemptydatamatpt(spikeTs, eventTimesByLoc{i}, timeLockStruct.window);
    [timeLockStruct.spdfByLoc(i,:),~,timeLockStruct.spdfErrByLoc(i,:)] = fixedPsth(timeLockStruct.spikeTimesByLoc{i}, kernelSigma, 2, timeLockStruct.t);
end

end
