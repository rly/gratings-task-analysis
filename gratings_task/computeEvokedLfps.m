function saveFileName = computeEvokedLfps(saveFileName, lfp, Fs, nLoc, UE)
% lfp = 1-D vector
% UE = useful events struct
% assert(iscolumn(lfp));

nTrials = numel(UE.cueOnset);

%% align spikes to cue onset
cueOnsetLfp.windowOffset = [-0.7 0.7];
cueOnsetLfp = alignLfpToEvents(lfp, UE.cueOnset, Fs, cueOnsetLfp);

% don't split by location -- just index into cueOnsetLfp with logical check
% for location -- this saves disk space by 2x
% cueOnsetLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     cueOnsetLfpByLoc{i} = alignLfpToEvents(lfp, UE.cueOnsetByLoc{i}, Fs, cueOnsetWindowOffset);
% end

%% align spikes to array onset
arrayOnsetLfp.windowOffset = [-0.7 0.7];
arrayOnsetRelBalLfp.windowOffset = [-0.7 0.7];
arrayOnsetHoldBalLfp.windowOffset = [-0.7 0.7];
arrayOnsetLfp = alignLfpToEvents(lfp, UE.arrayOnset, Fs, arrayOnsetLfp);
arrayOnsetRelBalLfp = alignLfpToEvents(lfp, UE.arrayOnsetRelBal, Fs, arrayOnsetRelBalLfp);
arrayOnsetHoldBalLfp = alignLfpToEvents(lfp, UE.arrayOnsetHoldBal, Fs, arrayOnsetHoldBalLfp);

% arrayOnsetLfpByLoc = cell(nLoc, 1);
% arrayOnsetRelLfpByLoc = cell(nLoc, 1);
% arrayOnsetHoldLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     arrayOnsetLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetByLoc{i}, Fs, arrayOnsetWindowOffset);
%     arrayOnsetRelLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetRelByLoc{i}, Fs, arrayOnsetWindowOffset);
%     arrayOnsetHoldLfpByLoc{i} = alignLfpToEvents(lfp, UE.arrayOnsetHoldByLoc{i}, Fs, arrayOnsetWindowOffset);
% end

%% align spikes to target dimming
targetDimBalLfp.windowOffset = [-0.7 0.7];
targetDimBalLfp = alignLfpToEvents(lfp, UE.targetDimBal, Fs, targetDimBalLfp);

% targetDimLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     targetDimLfpByLoc{i} = alignLfpToEvents(lfp, UE.targetDimByLoc{i}, Fs, targetDimWindowOffset);
% end


%% look at lever-locked responses and saccade-locked responses
% note that this includes both hold and release trials
assert(~isempty(UE.fixationAndLeverTimes))

%% align spikes to enter/exit fixation
enterFixationLfp.windowOffset = [-0.7 0.7];
enterFixationLfp = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue, Fs, enterFixationLfp);

% enterFixationLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     enterFixationLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstEnterFixationTimesPreCueByLoc{i}, Fs, enterFixationWindowOffset);
% end

exitFixationLfp.windowOffset = [-0.7 0.7];
exitFixationLfp = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice, Fs, exitFixationLfp);

% exitFixationLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     exitFixationLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuiceByLoc{i}, Fs, exitFixationWindowOffset);
% end

%% align spikes to lever press and release
% leverPressLfp.windowOffset = [-0.7 0.7];
% leverPressLfp = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverPressTimesPreCue, Fs, leverPressLfp);

% leverPressLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     leverPressLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverPressTimesPreCueByLoc{i}, Fs, leverPressWindowOffset);
% end

% leverReleaseLfp.windowOffset = [-0.7 0.7];
% leverReleaseLfp = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice, Fs, leverReleaseLfp);

% leverReleaseLfpByLoc = cell(nLoc, 1);
% for i = 1:nLoc
%     leverReleaseLfpByLoc{i} = alignLfpToEvents(lfp, UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuiceByLoc{i}, Fs, leverReleaseWindowOffset);
% end

%% calc time between motor events
assert(numel(UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue) == numel(UE.fixationAndLeverTimes.firstLeverPressTimesPreCue));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice));
assert(numel(UE.firstJuiceEvent) == numel(UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice));

% arrayOnsetRelToJuiceEventTime = UE.firstJuiceEvent(~UE.isHoldTrial) - UE.arrayOnsetRel; % almost same as UE.rt(~isHoldTrial)?
% targetDimToJuiceEventTime = UE.firstJuiceEvent(UE.isHoldTrial) - UE.targetDim; % almost same as UE.rt(isHoldTrial)?
% enterFixationToLeverPressTime = UE.fixationAndLeverTimes.firstLeverPressTimesPreCue - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
% enterFixationToCueOnsetTime = UE.cueOnset - UE.fixationAndLeverTimes.firstEnterFixationTimesPreCue;
% exitFixationToJuiceEventTime = UE.firstJuiceEvent - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
% exitFixationToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice - UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice;
% arrayOnsetRelToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
% targetDimToExitFixationTime = UE.fixationAndLeverTimes.firstExitFixationTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;
% arrayOnsetRelToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(~UE.isHoldTrial) - UE.arrayOnsetRel;
% targetDimToLeverReleaseTime = UE.fixationAndLeverTimes.firstLeverReleaseTimesAroundJuice(UE.isHoldTrial) - UE.targetDim;

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

%%
clear lfp;
tic;
fprintf('Saving to file: %s...\n', saveFileName);
save(saveFileName);
fileInfo = dir(saveFileName);
fprintf('%d MB file created. Took %0.1f minutes.\n', round(fileInfo.bytes/1024^2), toc/60);

end