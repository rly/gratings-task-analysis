function returnStruct = computeExpFiringRateBySpdf(analysisWindowOffset, timeLockedSpikesStruct)
% compute exponential decay / inverse decay fit for firing rate
returnStruct.windowOffset = analysisWindowOffset;

analysisWindow = timeLockedSpikesStruct.window(1) + analysisWindowOffset;
analysisWindowLogical = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow);

expFun = @(b,x) b(1) * exp(-abs(b(2)) * x) + b(3); % exponential decay or the inverse, depending on the sign of b(1)

% initialize b(1) with the slope so that the fitting can get the direction
% correct
allLM = fitlm(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdf(analysisWindowLogical));
allLMSlope = table2array(allLM.Coefficients(2,1)); % slope

expBeta0 = [-allLMSlope 1 mean(timeLockedSpikesStruct.spdf(analysisWindowLogical))];
returnStruct.allExpNLM = fitnlm(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdf(analysisWindowLogical), expFun, expBeta0);

% analysisWindow2 = timeLockedSpikesStruct.window(1) + [0 0.6];
% analysisWindow2Logical = getTimeLogicalWithTolerance(timeLockedSpikesStruct.t, analysisWindow2);
% analysisWindow3Logical = analysisWindowLogical | analysisWindow2Logical;
% figure;
% hold on;
% plot(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdf(analysisWindowLogical));
% plot(timeLockedSpikesStruct.t(analysisWindow2Logical), timeLockedSpikesStruct.spdf(analysisWindow2Logical));
% plot(timeLockedSpikesStruct.t(analysisWindow3Logical), predict(returnStruct.allExpNLM, timeLockedSpikesStruct.t(analysisWindow3Logical)'));

for i = 1:size(timeLockedSpikesStruct.spdfByLoc, 1)
    if any(~isnan(timeLockedSpikesStruct.spdfByLoc(i, analysisWindowLogical)))
        % initialize the per location fit with the grand fit coeffs
        expBeta0 = table2array(returnStruct.allExpNLM.Coefficients(1:3,1));
        returnStruct.byLocExpNLM{i} = fitnlm(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdfByLoc(i, analysisWindowLogical), expFun, expBeta0);
%         plot(timeLockedSpikesStruct.t(analysisWindowLogical), timeLockedSpikesStruct.spdfByLoc(i, analysisWindowLogical));
%         plot(timeLockedSpikesStruct.t(analysisWindow2Logical), timeLockedSpikesStruct.spdfByLoc(i, analysisWindow2Logical));
%         plot(timeLockedSpikesStruct.t(analysisWindow3Logical), predict(returnStruct.byLocExpNLM{i}, timeLockedSpikesStruct.t(analysisWindow3Logical)'));
    else
        returnStruct.byLocExpNLM{i} = NaN;
    end
end