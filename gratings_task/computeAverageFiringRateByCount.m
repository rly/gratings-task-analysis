function returnStruct = computeAverageFiringRateByCount(analysisWindowOffset, timeLockedSpikesStruct)

returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;

[returnStruct.all,returnStruct.allSD,returnStruct.trialRate] = computeAvgRateInWin(...
        timeLockedSpikesStruct.spikeTimes, analysisWindow);
    
for i = 1:size(timeLockedSpikesStruct.spikeTimesByLoc, 1)
    [returnStruct.byLoc(i),returnStruct.byLocSD(i),returnStruct.trialRateByLoc{i}] = computeAvgRateInWin(...
            timeLockedSpikesStruct.spikeTimesByLoc{i}, analysisWindow);
end