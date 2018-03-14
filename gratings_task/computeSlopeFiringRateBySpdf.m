function returnStruct = computeSlopeFiringRateBySpdf(analysisWindowOffset, timeLockedSpikesStruct)
% compute slope and P-value of slope of linear fit of SPDF in the given
% window
returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;
analysisWindowIndices = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow);

returnStruct.allLM = fitlm(timeLockedSpikesStruct.t(analysisWindowIndices), timeLockedSpikesStruct.spdf(analysisWindowIndices));
returnStruct.all = table2array(returnStruct.allLM.Coefficients(2,1)); % slope
returnStruct.allPValue = coefTest(returnStruct.allLM);

for i = 1:size(timeLockedSpikesStruct.spdfByLoc, 1)
    if any(~isnan(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices)))
        returnStruct.byLocLM{i} = fitlm(timeLockedSpikesStruct.t(analysisWindowIndices), timeLockedSpikesStruct.spdfByLoc(i, analysisWindowIndices));
        returnStruct.byLoc(i) = table2array(returnStruct.byLocLM{i}.Coefficients(2,1)); % slope
        returnStruct.byLocPValue(i) = coefTest(returnStruct.byLocLM{i});
    else
        returnStruct.byLocLM{i} = NaN;
        returnStruct.byLoc(i) = NaN;
        returnStruct.byLocPValue(i) = NaN;
    end
end