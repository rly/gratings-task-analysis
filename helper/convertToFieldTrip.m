function [data,trl] = convertToFieldTrip(D, lfps, UE)

timeFromCueOnset = -0.7;
timeFromFirstJuice = 0.7;

data = struct();
data.label = D.lfpNames;
data.fsample = D.lfpFs;
[aligned,t,sampleInfo,offset] = sliceLfpAroundTrials(lfps, UE.cueOnset, UE.firstJuiceEvent, ...
        D.lfpFs, timeFromCueOnset, timeFromFirstJuice);
nTrial = numel(t);
data.trial = aligned;
data.time = t;
data.sampleinfo = sampleInfo;

data.trialInfo = nan(nTrial, 3);
data.trialInfo(:,1) = UE.cueLoc;
data.trialInfo(:,2) = UE.isHoldTrial;
data.trialInfo(:,3) = UE.cueTargetDelayDur;
data.trialInfo(:,4) = UE.targetDimDelayDur;
data.trialInfo(:,5) = UE.rt;
data.trialInfo(:,6) = UE.fixationAndLeverTimes.firstExitFixationPolarAngle;

trl = nan(nTrial, 3);
trl(:,1:2) = data.sampleinfo;
trl(:,3) = offset;