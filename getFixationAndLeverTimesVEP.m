function [fixationAndLeverTimes,isMissingData] = getFixationAndLeverTimesVEP(D, cueOnset, firstJuiceEvent)

%% find mean fixation and lever voltage
% look at voltages -1.5 seconds to -0.5 seconds from first juice event
% should not be contaminated with fixation breaks or lever releases
% for at least session M20170311, -1.5 seconds is the largest time before
% first juice event to look before there are saccades or lever movement
direct1LockedToJuice = createdatamatc(D.adjDirects(1,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);
direct2LockedToJuice = createdatamatc(D.adjDirects(2,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);
direct3LockedToJuice = createdatamatc(D.adjDirects(3,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);

direct1LockedToJuice(direct1LockedToJuice < 0) = NaN;
direct2LockedToJuice(direct2LockedToJuice < 0) = NaN;
direct3LockedToJuice(direct3LockedToJuice < 0) = NaN;

fixationXVoltage = nanmean(direct1LockedToJuice(:));
fixationYVoltage = nanmean(direct2LockedToJuice(:));
leverPressedVoltage = nanmean(direct3LockedToJuice(:));

% rangeFixationXVoltage = max(abs(fixationXVoltage - direct1LockedToJuice(:))) + 2; % could make smaller or bigger
% rangeFixationYVoltage = max(abs(fixationYVoltage - direct2LockedToJuice(:))) + 2; % could make smaller or bigger
% rangeLeverPressedVoltage = max(abs(leverPressedVoltage - direct3LockedToJuice(:))) + 1; % could make smaller or bigger
rangeFixationXVoltage = min(10, 5*nanstd(direct1LockedToJuice(:))); % could make smaller or bigger
rangeFixationYVoltage = min(10, 5*nanstd(direct2LockedToJuice(:))); % could make smaller or bigger
rangeLeverPressedVoltage = max(4, 15*nanstd(direct3LockedToJuice(:))); % could make smaller or bigger

fprintf('\tX, Y bounds: %0.2f +/- %0.2f, %0.2f +/- %0.2f\n', ...
        fixationXVoltage, rangeFixationXVoltage, fixationYVoltage, rangeFixationYVoltage);

fprintf('\tLever bounds: %0.2f +/- %0.2f\n', ...
        leverPressedVoltage, rangeLeverPressedVoltage);

directFs = D.directFs;
removeEdgeTime = 0.5;

%%
direct1AtPresentationEnterFixation = nan(numel(D.events{2}), 1);
direct2AtPresentationEnterFixation = nan(numel(D.events{2}), 1);
for i = 1:numel(D.events{2})
    direct1AtPresentationEnterFixation(i) = D.adjDirects(1,round(D.events{2}(i)*1000-325));
    direct2AtPresentationEnterFixation(i) = D.adjDirects(2,round(D.events{2}(i)*1000-325));
end
figure
histogram(direct1AtPresentationEnterFixation)
figure
histogram(direct2AtPresentationEnterFixation)

%% find times of start fixation and exit fixation
holdFixationMin = 0.3;
holdFixationMinInds = round(directFs * holdFixationMin);
notFixating = abs(D.adjDirects(1,:) - fixationXVoltage) > rangeFixationXVoltage | ...
        abs(D.adjDirects(2,:) - fixationYVoltage) > rangeFixationYVoltage;

diffNotFixating = diff(notFixating);
% zero out ends to remove edge effects
diffNotFixating(1:removeEdgeTime*directFs) = 0;
diffNotFixating(end-removeEdgeTime*directFs+1:end) = 0;
enterFixation = diffNotFixating == -1;
exitFixation = diffNotFixating == 1;

enterFixationInds = find(enterFixation)';
badEnterFixation = false(size(enterFixationInds));
for i = 1:numel(enterFixationInds)
    badEnterFixation(i) = any(notFixating(enterFixationInds(i)+1:enterFixationInds(i)+holdFixationMinInds) == 1);
end

exitFixationInds = find(exitFixation)';
badExitFixation = false(size(exitFixationInds));
for i = 1:numel(exitFixationInds)
    badExitFixation(i) = any(notFixating(exitFixationInds(i)-holdFixationMinInds+1:exitFixationInds(i)) == 1);
end

enterFixationTimes = enterFixationInds(~badEnterFixation) / directFs; % assume starting at 1
exitFixationTimes = exitFixationInds(~badExitFixation) / directFs; % assume starting at 1

% get direction of saccade centered on center of FP -- NOT actual fixation
% point. so this is just a rough estimate.
% exitFixationVector = [D.adjDirects(2,exitFixationInds(~badExitFixation))' - fixationYVoltage ...
%         D.adjDirects(1,exitFixationInds(~badExitFixation))' - fixationXVoltage];
% exitFixationPolarAngle = cart2pol(exitFixationVector(:,1), exitFixationVector(:,2));

%% find times of lever press and release
holdLeverMin = 0.3;
holdLeverMinInds = round(directFs * holdLeverMin);
notLeverHeld = abs(D.adjDirects(3,:) - leverPressedVoltage) > rangeLeverPressedVoltage;

diffNotLeverHeld = diff(notLeverHeld);
% zero out ends to remove edge effects
diffNotLeverHeld(1:removeEdgeTime*directFs) = 0;
diffNotLeverHeld(end-removeEdgeTime*directFs+1:end) = 0;
leverPress = diffNotLeverHeld == -1;
leverRelease = diffNotLeverHeld == 1;

leverPressInds = find(leverPress)';
badLeverPress = false(size(leverPressInds));
for i = 1:numel(leverPressInds)
    badLeverPress(i) = any(notLeverHeld(leverPressInds(i)+1:leverPressInds(i)+holdLeverMinInds) == 1);
end

leverReleaseInds = find(leverRelease)';
badLeverRelease = false(size(leverReleaseInds));
for i = 1:numel(leverReleaseInds)
    badLeverRelease(i) = any(notLeverHeld(leverReleaseInds(i)-holdLeverMinInds+1:leverReleaseInds(i)) == 1);
end

leverPressTimes = leverPressInds(~badLeverPress) / directFs; % assume starting at 1
leverReleaseTimes = leverReleaseInds(~badLeverRelease) / directFs; % assume starting at 1

assert(numel(cueOnset) == numel(firstJuiceEvent))

%% find enter/exit fixation times around trial start/juice
cueOnset = D.events{2}
maxEnterFixationTimesToCueOnset = 4.05;
minEnterFixationTimesToCueOnset = 0.325;
isEnterFixationTimesFirstPreCue = false(size(enterFixationTimes));
firstEnterFixationTimesPreCue = nan(size(cueOnset)); % use precue marker for flash mapping
for i = 1:numel(cueOnset)
    % enter fixation must precede trial start by 325+ ms
    enterFixationTimesToCueOnset = cueOnset(i) - enterFixationTimes;
    enterFixationTimesGoodInd = find(enterFixationTimesToCueOnset > minEnterFixationTimesToCueOnset & ...
            enterFixationTimesToCueOnset < maxEnterFixationTimesToCueOnset, 1, 'last');
    if ~isempty(enterFixationTimesGoodInd)
        if isEnterFixationTimesFirstPreCue(enterFixationTimesGoodInd)
            fprintf('\tEnter fixation event already assigned to cue onset event: %d, %d\n', i, enterFixationTimesGoodInd);
        end
        isEnterFixationTimesFirstPreCue(enterFixationTimesGoodInd) = 1;
        firstEnterFixationTimesPreCue(i) = enterFixationTimes(enterFixationTimesGoodInd);
    else
        fprintf('\tCue onset event without associated enter fixation event: %d\n', i);
    end
end
% firstEnterFixationTimesPreCue = enterFixationTimes(isEnterFixationTimesFirstPreCue);
otherEnterFixationTimes = enterFixationTimes(~isEnterFixationTimesFirstPreCue);

maxDiffExitFixationTimesToJuice = 1.0; % TODO deal with rare case where exit fixation is too long after juice
isExitFixationTimesFirstAroundJuice = false(size(exitFixationTimes));
firstExitFixationTimesAroundJuice = nan(size(firstJuiceEvent));
firstExitFixationVector = nan(numel(firstJuiceEvent), 2);
firstExitFixationPolarAngle = nan(size(firstJuiceEvent));
for i = 1:numel(firstJuiceEvent)
    % exit fixation may happen anytime within say 1 second around juice
    % event. find the earliest
    exitFixationTimesToJuice = firstJuiceEvent(i) - exitFixationTimes;
    exitFixationTimesGoodInd = find(abs(exitFixationTimesToJuice) < maxDiffExitFixationTimesToJuice, 1, 'first');
    if ~isempty(exitFixationTimesGoodInd)
        if isExitFixationTimesFirstAroundJuice(exitFixationTimesGoodInd)
            fprintf('\tExit fixation event already assigned to first juice event: %d, %d\n', i, exitFixationTimesGoodInd);
        end
        isExitFixationTimesFirstAroundJuice(exitFixationTimesGoodInd) = 1;
        firstExitFixationTimesAroundJuice(i) = exitFixationTimes(exitFixationTimesGoodInd);
%         firstExitFixationVector(i,:) = exitFixationVector(isExitFixationTimesFirstAroundJuice(exitFixationTimesGoodInd),:);
%         disp(firstExitFixationVector(i,:))
%         firstExitFixationPolarAngle(i) = exitFixationPolarAngle(isExitFixationTimesFirstAroundJuice(exitFixationTimesGoodInd));
    else
        fprintf('\tFirst juice event without associated exit fixation event: %d\n', i);
    end
end
% firstExitFixationTimesAroundJuice = exitFixationTimes(isExitFixationTimesFirstAroundJuice);
otherExitFixationTimes = exitFixationTimes(~isExitFixationTimesFirstAroundJuice);

% firstExitFixationVector = exitFixationVector(isExitFixationTimesFirstAroundJuice,:);
% firstExitFixationPolarAngle = exitFixationPolarAngle(isExitFixationTimesFirstAroundJuice);

%% find lever press/release around trial start/juice
maxLeverPressTimesToCueOnset = 4;
minLeverPressTimesToCueOnset = 0.35;
isLeverPressTimesFirstPreCue = false(size(leverPressTimes));
firstLeverPressTimesPreCue = nan(size(cueOnset));
for i = 1:numel(cueOnset)
    leverPressTimesToCueOnset = cueOnset(i) - leverPressTimes;
    leverPressTimesGoodInd = find(leverPressTimesToCueOnset > minLeverPressTimesToCueOnset & ...
            leverPressTimesToCueOnset < maxLeverPressTimesToCueOnset, 1, 'last');
    if ~isempty(leverPressTimesGoodInd)
        if isLeverPressTimesFirstPreCue(leverPressTimesGoodInd)
            % this can happen if I initiate a trial and he presses the
            % lever later or doesn't press the lever
            fprintf('\tLever press event already assigned to cue onset event: %d, %d\n', i, leverPressTimesGoodInd);
        end
        isLeverPressTimesFirstPreCue(leverPressTimesGoodInd) = 1;
        firstLeverPressTimesPreCue(i) = leverPressTimes(leverPressTimesGoodInd);
    else
        fprintf('\tCue onset event without associated lever press event: %d\n', i);
    end
end
% firstLeverPressTimesPreCue = leverPressTimes(isLeverPressTimesFirstPreCue);
otherLeverPressTimes = leverPressTimes(~isLeverPressTimesFirstPreCue);

maxDiffLeverReleaseTimesToJuice = 1;
isLeverReleaseTimesFirstAroundJuice = false(size(leverReleaseTimes));
firstLeverReleaseTimesAroundJuice = nan(size(firstJuiceEvent));
for i = 1:numel(firstJuiceEvent)
    % lever release must precede juice
    leverReleaseTimesToJuice = firstJuiceEvent(i) - leverReleaseTimes;
    leverReleaseTimesGoodInd = find(leverReleaseTimesToJuice >= 0 & ...
            leverReleaseTimesToJuice < maxDiffLeverReleaseTimesToJuice, 1, 'last');
    if ~isempty(leverReleaseTimesGoodInd)
        if isLeverReleaseTimesFirstAroundJuice(leverReleaseTimesGoodInd)
            fprintf('\tLever release event already assigned to first juice event: %d, %d\n', i, leverReleaseTimesGoodInd);
        end
        isLeverReleaseTimesFirstAroundJuice(leverReleaseTimesGoodInd) = 1;
        firstLeverReleaseTimesAroundJuice(i) = leverReleaseTimes(leverReleaseTimesGoodInd);
    else
        fprintf('\tFirst juice event without associated lever release event: %d\n', i);
    end
end
% firstLeverReleaseTimesAroundJuice = leverReleaseTimes(isLeverReleaseTimesFirstAroundJuice);
otherLeverReleaseTimes = leverReleaseTimes(~isLeverReleaseTimesFirstAroundJuice);

%%
isMissingData = isnan(firstEnterFixationTimesPreCue) | ...
        isnan(firstLeverPressTimesPreCue) | ...
        isnan(firstExitFixationTimesAroundJuice) | ...
        isnan(firstLeverReleaseTimesAroundJuice);
if sum(isMissingData) > 0
    fprintf('Missing matching lever/fixation data for %d trials.\n', sum(isMissingData));
end

assert(all([numel(firstEnterFixationTimesPreCue) ...
        numel(firstLeverPressTimesPreCue) ...
        numel(firstExitFixationTimesAroundJuice)]));
    
firstEnterFixationTimesPreCue = firstEnterFixationTimesPreCue(~isMissingData);
firstLeverPressTimesPreCue = firstLeverPressTimesPreCue(~isMissingData);
firstExitFixationTimesAroundJuice = firstExitFixationTimesAroundJuice(~isMissingData);
firstLeverReleaseTimesAroundJuice = firstLeverReleaseTimesAroundJuice(~isMissingData);

%%
voltages = var2struct(fixationXVoltage, fixationYVoltage, leverPressedVoltage, ...
        rangeFixationXVoltage, rangeFixationYVoltage, rangeLeverPressedVoltage);

fixationAndLeverTimes = var2struct(firstEnterFixationTimesPreCue, otherEnterFixationTimes, ...
        firstExitFixationTimesAroundJuice, otherExitFixationTimes, ...
        firstLeverPressTimesPreCue, otherLeverPressTimes, ...
        firstLeverReleaseTimesAroundJuice, otherLeverReleaseTimes, ...
        isMissingData, ...
        voltages);