function returnStruct = computeAverageFiringRateBySpdf(analysisWindowOffset, timeLockedSpikesStruct)

returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;
analysisWindowIndices = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow);

returnStruct.all = mean(timeLockedSpikesStruct.spdf(analysisWindowIndices));
returnStruct.allSDOverTime = std(timeLockedSpikesStruct.spdf(analysisWindowIndices));

% TODO add bootstrap error

for i = 1:size(timeLockedSpikesStruct.spdfByLoc, 1)
    returnStruct.byLoc(i) = mean(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
    returnStruct.byLocSDOverTime(i) = std(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
end