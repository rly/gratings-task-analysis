function returnStruct = computeAverageFiringRateBySpdf(analysisWindowOffset, timeLockedSpikesStruct, timeLockedSpikesErrorStruct)

returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;
analysisWindowIndices = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow);

returnStruct.all = mean(timeLockedSpikesStruct.spdf(analysisWindowIndices));
returnStruct.allSDOverTime = std(timeLockedSpikesStruct.spdf(analysisWindowIndices));
returnStruct.allMax = max(timeLockedSpikesStruct.spdf(analysisWindowIndices));

% TODO add bootstrap error

for i = 1:size(timeLockedSpikesStruct.spdfByLoc, 1)
    returnStruct.byLoc(i) = mean(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
    returnStruct.byLocSDOverTime(i) = std(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
    returnStruct.byLocMax(i) = max(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
end

if nargin == 3
    returnStruct.allError = mean(timeLockedSpikesErrorStruct.spdf(analysisWindowIndices));
    returnStruct.allSDOverTimeError = std(timeLockedSpikesErrorStruct.spdf(analysisWindowIndices));
    returnStruct.allMaxError = max(timeLockedSpikesErrorStruct.spdf(analysisWindowIndices));

    % TODO add bootstrap error

    for i = 1:size(timeLockedSpikesErrorStruct.spdfByLoc, 1)
        returnStruct.byLocError(i) = mean(timeLockedSpikesErrorStruct.spdfByLoc(i, analysisWindowIndices));
        returnStruct.byLocSDOverTimeError(i) = std(timeLockedSpikesErrorStruct.spdfByLoc(i, analysisWindowIndices));
        returnStruct.byLocMaxError(i) = max(timeLockedSpikesErrorStruct.spdfByLoc(i, analysisWindowIndices));
    end
end