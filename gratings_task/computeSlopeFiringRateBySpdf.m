function returnStruct = computeSlopeFiringRateBySpdf(analysisWindowOffset, timeLockedSpikesStruct)
% compute slope and P-value of slope of linear fit of SPDF in the given
% window
returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;
analysisWindowLogical = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow);

returnStruct.allLM = fitlm(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdf(analysisWindowLogical));
returnStruct.all = table2array(returnStruct.allLM.Coefficients(2,1)); % slope
returnStruct.allPValue = coefTest(returnStruct.allLM);

for i = 1:size(timeLockedSpikesStruct.spdfByLoc, 1)
    if any(~isnan(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowLogical)))
        returnStruct.byLocLM{i} = fitlm(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdfByLoc(i, analysisWindowLogical));
        returnStruct.byLoc(i) = table2array(returnStruct.byLocLM{i}.Coefficients(2,1)); % slope
        returnStruct.byLocPValue(i) = coefTest(returnStruct.byLocLM{i});
    else
        returnStruct.byLocLM{i} = NaN;
        returnStruct.byLoc(i) = NaN;
        returnStruct.byLocPValue(i) = NaN;
    end
end