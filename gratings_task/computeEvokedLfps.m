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