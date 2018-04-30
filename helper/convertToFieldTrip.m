function [data,trl] = convertToFieldTrip(D, lfps, UE)
% D is struct of recording information
% lfps is nChannel x nSamples

timeFromCueOnset = -0.7;
timeFromFirstJuice = 0.7;

data = struct();
% shorten label
data.label = cellfun(@(x) x(strfind(x, 'FP'):end), D.lfpNames, 'UniformOutput', false);
data.fsample = D.lfpFs;
[data.trial,data.time,data.sampleinfo,offset] = sliceLfpAroundTrials(lfps, UE.cueOnset, UE.firstJuiceEvent, ...
        D.lfpFs, timeFromCueOnset, timeFromFirstJuice);

% set up data.trialinfo with condition information
nTrial = numel(data.time);
data.trialinfo = nan(nTrial, 6);
data.trialinfo(:,1) = UE.cueLoc;
data.trialinfo(:,2) = UE.isHoldTrial;
data.trialinfo(:,3) = UE.cueTargetDelayDur;
data.trialinfo(:,4) = UE.targetDimDelayDur;
data.trialinfo(:,5) = UE.rt;
data.trialinfo(:,6) = UE.fixationAndLeverTimes.firstExitFixationPolarAngle;

% set up trl if needed
trl = nan(nTrial, 3);
trl(:,1:2) = data.sampleinfo;
trl(:,3) = offset;