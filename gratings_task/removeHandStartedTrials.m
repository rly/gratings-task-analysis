function [cueOnsetClean,trialsToKeep] = removeHandStartedTrials(D, cueOnset, firstJuiceEvent)
% if I start the trial using the keyboard, monkey may not press the lever
% beforehand. remove those trials.

%% find lever voltage
% look at voltages -1.0 seconds to -0.5 seconds from first juice event
% should not be contaminated with fixation breaks or lever releases
direct3LockedToJuice = createdatamatc(D.adjDirects(3,:)', firstJuiceEvent, D.directFs, [1 -0.5]);

leverPressedVoltage = mean(direct3LockedToJuice(:));
rangeLeverPressedVoltage = max(4, 15*std(direct3LockedToJuice(:))); % could make smaller or bigger

directFs = D.directFs;
removeEdgeTime = 0.5;

%% find times of lever press
holdLeverMin = 0.3;
holdLeverMinInds = round(directFs * holdLeverMin);
notLeverHeld = abs(D.adjDirects(3,:) - leverPressedVoltage) > rangeLeverPressedVoltage;

diffNotLeverHeld = diff(notLeverHeld);
% zero out ends to remove edge effects
diffNotLeverHeld(1:removeEdgeTime*directFs) = 0;
diffNotLeverHeld(end-removeEdgeTime*directFs+1:end) = 0;
leverPress = diffNotLeverHeld == -1;

leverPressInds = find(leverPress);
badLeverPress = false(size(leverPressInds));
for i = 1:numel(leverPressInds)
    badLeverPress(i) = any(notLeverHeld(leverPressInds(i)+1:leverPressInds(i)+holdLeverMinInds) == 1);
end

leverPressTimes = leverPressInds(~badLeverPress) / directFs; % assume starting at 1
    
%%
isLeverPressTimesFirstPreCue = false(size(leverPressTimes));
maxLeverPressTimesToCueOnset = 5;
minLeverPressTimesToCueOnset = 0.35;
trialsToKeep = true(size(cueOnset));
for i = 1:numel(cueOnset)
    leverPressTimesToCueOnset = cueOnset(i) - leverPressTimes;
    leverPressTimesGoodInd = find(leverPressTimesToCueOnset > minLeverPressTimesToCueOnset & ...
            leverPressTimesToCueOnset < maxLeverPressTimesToCueOnset, 1, 'last');
    if ~isempty(leverPressTimesGoodInd)
        if isLeverPressTimesFirstPreCue(leverPressTimesGoodInd)
            fprintf('\tLever press event already assigned to cue onset event: %d, %d\n', i, leverPressTimesGoodInd);
            trialsToKeep(i) = false;
        end
    else
        fprintf('\tCue onset event without associated lever press event: %d\n', i);
        trialsToKeep(i) = false;
    end
end

cueOnsetClean = cueOnset(trialsToKeep);


