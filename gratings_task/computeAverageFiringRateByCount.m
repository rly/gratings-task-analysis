function returnStruct = computeAverageFiringRateByCount(analysisWindowOffset, timeLockedSpikesStruct, timeLockedSpikesErrorStruct)

returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;

[returnStruct.all,returnStruct.allSD,returnStruct.trialRate,returnStruct.trialCount] = computeAvgRateInWin(...
        timeLockedSpikesStruct.spikeTimes, analysisWindow);
    
for i = 1:size(timeLockedSpikesStruct.spikeTimesByLoc, 1)
    [returnStruct.byLoc(i),returnStruct.byLocSD(i),returnStruct.trialRateByLoc{i},returnStruct.trialCountByLoc{i}] = computeAvgRateInWin(...
            timeLockedSpikesStruct.spikeTimesByLoc{i}, analysisWindow);
end

if nargin == 3
    [returnStruct.allError,returnStruct.allSDError,returnStruct.trialRateError,returnStruct.trialCountError] = computeAvgRateInWin(...
            timeLockedSpikesErrorStruct.spikeTimes, analysisWindow);
    
    for i = 1:size(timeLockedSpikesErrorStruct.spikeTimesByLoc, 1)
        [returnStruct.byLocError(i),returnStruct.byLocSDError(i),returnStruct.trialRateByLocError{i},returnStruct.trialCountByLocError{i}] = computeAvgRateInWin(...
                timeLockedSpikesErrorStruct.spikeTimesByLoc{i}, analysisWindow);
    end
end