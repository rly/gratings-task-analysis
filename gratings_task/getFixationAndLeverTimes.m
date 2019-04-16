function fixationAndLeverTimes = getFixationAndLeverTimes(D, cueOnset, firstJuiceEvent, cueLoc, isHoldTrial, nLoc)

%% find mean fixation and lever voltage
% look at voltages -1.5 seconds to -0.5 seconds from first juice event
% should not be contaminated with fixation breaks or lever releases
% for at least session M20170311, -1.5 seconds is the largest time before
% first juice event to look before there are saccades or lever movement
direct1LockedToJuice = createdatamatc(D.adjDirects(1,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);
direct2LockedToJuice = createdatamatc(D.adjDirects(2,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);
direct3LockedToJuice = createdatamatc(D.adjDirects(3,:)', firstJuiceEvent, D.directFs, [1.5 -0.5]);

fixationXVoltage = mean(direct1LockedToJuice(:));
fixationYVoltage = mean(direct2LockedToJuice(:));
leverPressedVoltage = mean(direct3LockedToJuice(:));

% rangeFixationXVoltage = max(abs(fixationXVoltage - direct1LockedToJuice(:))) + 2; % could make smaller or bigger
% rangeFixationYVoltage = max(abs(fixationYVoltage - direct2LockedToJuice(:))) + 2; % could make smaller or bigger
% rangeLeverPressedVoltage = max(abs(leverPressedVoltage - direct3LockedToJuice(:))) + 1; % could make smaller or bigger
rangeFixationXVoltage = min(6, 5*std(direct1LockedToJuice(:))); % could make smaller or bigger
rangeFixationYVoltage = min(6, 5*std(direct2LockedToJuice(:))); % could make smaller or bigger
rangeLeverPressedVoltage = max(4, 15*std(direct3LockedToJuice(:))); % could make smaller or bigger

fprintf('\tX, Y bounds: %0.2f +/- %0.2f, %0.2f +/- %0.2f\n', ...
        fixationXVoltage, rangeFixationXVoltage, fixationYVoltage, rangeFixationYVoltage);

fprintf('\tLever bounds: %0.2f +/- %0.2f\n', ...
        leverPressedVoltage, rangeLeverPressedVoltage);

directFs = D.directFs;
removeEdgeTime = 0.5;

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
exitFixationVector = [D.adjDirects(2,exitFixationInds(~badExitFixation))' - fixationYVoltage ...
        D.adjDirects(1,exitFixationInds(~badExitFixation))' - fixationXVoltage];
exitFixationPolarAngle = cart2pol(exitFixationVector(:,1), exitFixationVector(:,2));

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

%% find enter/exit fixation times around trial start/juice
isEnterFixationTimesFirstPreCue = false(size(enterFixationTimes));
maxEnterFixationTimesToCueOnset = 4.05;
minEnterFixationTimesToCueOnset = 0.35;
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
    else
        fprintf('\tCue onset event without associated enter fixation event: %d\n', i);
    end
end
firstEnterFixationTimesPreCue = enterFixationTimes(isEnterFixationTimesFirstPreCue);
otherEnterFixationTimes = enterFixationTimes(~isEnterFixationTimesFirstPreCue);

maxDiffExitFixationTimesToJuice = 3.0; % TODO deal with rare case where exit fixation is too long after juice
isExitFixationTimesFirstAroundJuice = false(size(exitFixationTimes));
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
    else
        fprintf('\tFirst juice event without associated exit fixation event: %d\n', i);
    end
end
firstExitFixationTimesAroundJuice = exitFixationTimes(isExitFixationTimesFirstAroundJuice);
otherExitFixationTimes = exitFixationTimes(~isExitFixationTimesFirstAroundJuice);
firstExitFixationVector = exitFixationVector(isExitFixationTimesFirstAroundJuice,:);
firstExitFixationPolarAngle = exitFixationPolarAngle(isExitFixationTimesFirstAroundJuice);
assert(~any(any(isnan(firstExitFixationVector))));
firstExitFixationVector = firstExitFixationVector ./ sqrt(sum(firstExitFixationVector.^2, 2)); % normalize

leftPolarAngleBounds = [-5*pi/8 5*pi/8];
rightPolarAngleBounds = [-3*pi/8 3*pi/8];
isExitFixationLeft = firstExitFixationPolarAngle < leftPolarAngleBounds(1) | firstExitFixationPolarAngle > leftPolarAngleBounds(2);
isExitFixationRight = firstExitFixationPolarAngle > rightPolarAngleBounds(1) & firstExitFixationPolarAngle < rightPolarAngleBounds(2);
firstExitFixationLeftTimesAroundJuice = firstExitFixationTimesAroundJuice(isExitFixationLeft);
firstExitFixationRightTimesAroundJuice = firstExitFixationTimesAroundJuice(isExitFixationRight);

%% find lever press/release around trial start/juice
isLeverPressTimesFirstPreCue = false(size(leverPressTimes));
maxLeverPressTimesToCueOnset = 6;
minLeverPressTimesToCueOnset = 0.35;
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
    else
        fprintf('\tCue onset event without associated lever press event: %d\n', i);
    end
end
firstLeverPressTimesPreCue = leverPressTimes(isLeverPressTimesFirstPreCue);
otherLeverPressTimes = leverPressTimes(~isLeverPressTimesFirstPreCue);

maxDiffLeverReleaseTimesToJuice = 1;
isLeverReleaseTimesFirstAroundJuice = false(size(leverReleaseTimes));
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
    else
        fprintf('\tFirst juice event without associated lever release event: %d\n', i);
    end
end
firstLeverReleaseTimesAroundJuice = leverReleaseTimes(isLeverReleaseTimesFirstAroundJuice);
otherLeverReleaseTimes = leverReleaseTimes(~isLeverReleaseTimesFirstAroundJuice);

%%
assert(all([numel(firstEnterFixationTimesPreCue) ...
        numel(firstLeverPressTimesPreCue) ...
        numel(firstExitFixationTimesAroundJuice) ...
        numel(firstLeverReleaseTimesAroundJuice)] == numel(cueLoc)));


firstEnterFixationHoldTimesPreCue = firstEnterFixationTimesPreCue(isHoldTrial);
firstEnterFixationRelTimesPreCue = firstEnterFixationTimesPreCue(~isHoldTrial);
firstExitFixationHoldTimesAroundJuice = firstExitFixationTimesAroundJuice(isHoldTrial);
firstExitFixationRelTimesAroundJuice = firstExitFixationTimesAroundJuice(~isHoldTrial);


firstEnterFixationTimesPreCueByLoc = cell(nLoc, 1);
firstLeverPressTimesPreCueByLoc = cell(nLoc, 1);
firstExitFixationTimesAroundJuiceByLoc = cell(nLoc, 1);
firstLeverReleaseTimesAroundJuiceByLoc = cell(nLoc, 1);
firstExitFixationVectorByLoc = cell(nLoc, 1);
firstExitFixationPolarAngleByLoc = cell(nLoc, 1);
firstExitFixationLeftTimesAroundJuiceByLoc = cell(nLoc, 1);
firstExitFixationRightTimesAroundJuiceByLoc = cell(nLoc, 1);
firstEnterFixationHoldTimesPreCueByLoc = cell(nLoc, 1);
firstEnterFixationRelTimesPreCueByLoc = cell(nLoc, 1);
firstExitFixationHoldTimesAroundJuiceByLoc = cell(nLoc, 1);
firstExitFixationRelTimesAroundJuiceByLoc = cell(nLoc, 1);
for i = 1:nLoc
    firstEnterFixationTimesPreCueByLoc{i} = firstEnterFixationTimesPreCue(cueLoc == i);
    firstLeverPressTimesPreCueByLoc{i} = firstLeverPressTimesPreCue(cueLoc == i);
    firstExitFixationTimesAroundJuiceByLoc{i} = firstExitFixationTimesAroundJuice(cueLoc == i);
    firstLeverReleaseTimesAroundJuiceByLoc{i} = firstLeverReleaseTimesAroundJuice(cueLoc == i);
    firstExitFixationVectorByLoc = firstExitFixationVector(cueLoc == i,:);
    firstExitFixationPolarAngleByLoc = firstExitFixationPolarAngle(cueLoc == i,:);
    firstExitFixationLeftTimesAroundJuiceByLoc{i} = firstExitFixationTimesAroundJuice(cueLoc == i & isExitFixationLeft);
    firstExitFixationRightTimesAroundJuiceByLoc{i} = firstExitFixationTimesAroundJuice(cueLoc == i  & isExitFixationRight);
    firstEnterFixationHoldTimesPreCueByLoc{i} = firstEnterFixationTimesPreCue(cueLoc == i & isHoldTrial);
    firstEnterFixationRelTimesPreCueByLoc{i} = firstEnterFixationTimesPreCue(cueLoc == i  & ~isHoldTrial);
    firstExitFixationHoldTimesAroundJuiceByLoc{i} = firstExitFixationTimesAroundJuice(cueLoc == i & isHoldTrial);
    firstExitFixationRelTimesAroundJuiceByLoc{i} = firstExitFixationTimesAroundJuice(cueLoc == i  & ~isHoldTrial);
end

%%
voltages = var2struct(fixationXVoltage, fixationYVoltage, leverPressedVoltage, ...
        rangeFixationXVoltage, rangeFixationYVoltage, rangeLeverPressedVoltage);

fixationAndLeverTimes = var2struct(firstEnterFixationTimesPreCue, otherEnterFixationTimes, ...
        firstExitFixationTimesAroundJuice, otherExitFixationTimes, ...
        firstLeverPressTimesPreCue, otherLeverPressTimes, ...
        firstLeverReleaseTimesAroundJuice, otherLeverReleaseTimes, ...
        firstEnterFixationHoldTimesPreCue, ...
        firstEnterFixationRelTimesPreCue, ...
        firstExitFixationHoldTimesAroundJuice, ...
        firstExitFixationRelTimesAroundJuice, ...
        firstEnterFixationTimesPreCueByLoc, ...
        firstLeverPressTimesPreCueByLoc, ...
        firstExitFixationTimesAroundJuiceByLoc, ...
        firstLeverReleaseTimesAroundJuiceByLoc, ...
        firstExitFixationVector, ...
        firstExitFixationVectorByLoc, ...
        firstExitFixationPolarAngle, ...
        firstExitFixationPolarAngleByLoc, ...
        firstExitFixationLeftTimesAroundJuice, ...
        firstExitFixationRightTimesAroundJuice, ...
        firstExitFixationLeftTimesAroundJuiceByLoc, ...
        firstExitFixationRightTimesAroundJuiceByLoc, ...
        firstEnterFixationHoldTimesPreCueByLoc, ...
        firstEnterFixationRelTimesPreCueByLoc, ...
        firstExitFixationHoldTimesAroundJuiceByLoc, ...
        firstExitFixationRelTimesAroundJuiceByLoc, ...
        voltages);