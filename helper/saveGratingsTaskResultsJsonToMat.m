function trialStructs = saveGratingsTaskResultsJsonToMat(logDir, logIndices, saveFileName)
% make sure that all JSON files are complete. if the block was aborted
% early, the JSON file may not be well formatted and will throw an error in
% loadjson.
% also, the original JSON results files had lines like:
% "shapeP1": H
% "shapeP2": R
% that had to be fixed to add quotes around H and R.
% run reformatGratingTaskResults.sh to fix that.

% should be sorted using my naming scheme
allLogFiles = dir(sprintf('%s/gratings_task_results_*.json', logDir));

trialStructs = {};
for i = 1:numel(allLogFiles)
    if ~any(i == logIndices)
        continue;
    end
    
    fileName = sprintf('%s/%s', allLogFiles(i).folder, allLogFiles(i).name);
    if ~exist(fileName, 'file')
        error('Cannot find file %s', fileName);
    end
    fprintf('Loading JSON %s...\n', fileName);
    trialStructsBlock = loadjson(fileName); % THIS IS SLOW!!!!
    trialStructs = [trialStructs trialStructsBlock];
end
clear i fileName trialStructsBlock; 

% flow:
% pre-array-on-early-release: trial_state == PRE_CUE (-105) || 
%               trial_state == CUE_ON (-106) || 
% 				trial_state == CUE_ARRAY_INTERVAL (-107)
%
% array-on-early-release: trial_state == ARRAY_ON (-108)
% 
% just-early-release: trial_state == TARGET_DIMMED (-110) || 
% 				trial_state == ARRAY_ON_RELEASE_SHAPE (-109)
%
% correct-response: trial_state == WAITING_FOR_RESPONSE (-112) || 
% 				trial_state == WAITING_FOR_RESPONSE_RELEASE_SHAPE (-117)
% 
% late-response: trial_state == PAST_RESPONSE_PERIOD (-113) || 
% 				trial_state == PAST_RESPONSE_PERIOD_RELEASE_SHAPE (-118) ||
% 				trial_state == RELEASE_SHAPE_DIMMED (-114)
%
% pre-array-eye-error: other eye error
%
% release-shape-on-eye-error: trial_state == ARRAY_ON_RELEASE_SHAPE (-109) ||
% 					trial_state == WAITING_FOR_RESPONSE_RELEASE_SHAPE (-117) ||
% 					trial_state == PAST_RESPONSE_PERIOD_RELEASE_SHAPE (-118) ||
% 					trial_state == RELEASE_SHAPE_DIMMED (-114)
%
% hold-shape-on-eye-error: trial_state == ARRAY_ON (-108) || 
% 					trial_state == TARGET_DIMMED (-110) ||
% 					trial_state == WAITING_FOR_RESPONSE (-112) ||
% 					trial_state == PAST_RESPONSE_PERIOD (-113)

trialResults = cellfun(@(x) x.response.type, trialStructs, 'UniformOutput', false);
nTrials = numel(trialResults);
isPreArrayOnEarlyReleaseError = strcmp(trialResults, 'pre-array-on-early-release');
isArrayOnEarlyReleaseError = strcmp(trialResults, 'array-on-early-release'); % hold trials only
isJustEarlyReleaseError = strcmp(trialResults, 'just-early-release'); % can be for hold or release
isCorrect = strcmp(trialResults, 'correct-response');
isLateResponseError = strcmp(trialResults, 'late-response');
isPreArrayEyeError = strcmp(trialResults, 'pre-array-eye-error');
isReleaseShapeOnEyeError = strcmp(trialResults, 'release-shape-on-eye-error');
isHoldShapeOnEyeError = strcmp(trialResults, 'hold-shape-on-eye-error');
assert(nTrials == (sum(isCorrect) + sum(isPreArrayEyeError) + sum(isHoldShapeOnEyeError) + ...
        sum(isReleaseShapeOnEyeError) + sum(isJustEarlyReleaseError) + sum(isLateResponseError) + ...
        sum(isPreArrayOnEarlyReleaseError) + sum(isArrayOnEarlyReleaseError)));

nTrialsArrayOnward = nTrials - sum(isPreArrayOnEarlyReleaseError) - sum(isPreArrayEyeError);

% verify two different methods of log files
correctTrialParams = readPresentationLogCorrectTrial(logDir, 'gratings_task_eye_tracking_', logIndices);
verifyGratingsTaskLogsAndJsonCorrect(correctTrialParams, trialStructs(isCorrect));

fprintf('%d/%d trials correct = %d%%\n', sum(isCorrect), nTrials, round(sum(isCorrect) / nTrials * 100));
fprintf('%d/%d trials from array onset onward correct = %d%%\n', sum(isCorrect), nTrialsArrayOnward, ...
        round(sum(isCorrect) / nTrialsArrayOnward * 100));

if isfield(trialStructs{1}, 'shapeP2')
    arrayShapes = cellfun(@(x) [x.shapeP1 x.shapeP2 x.shapeP3 x.shapeP4], trialStructs, 'UniformOutput', false);
else
    assert(isfield(trialStructs{1}, 'shapeP1'));
    arrayShapes = cellfun(@(x) x.shapeP1, trialStructs, 'UniformOutput', false);
end

fprintf('Saving to %s...\n', saveFileName);
save(saveFileName);