S1 = load('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data//MUA_GRATINGS_SUMMARY//M20170311-sessionInd7-muaAnalysisSummaryData-v11.mat')
S2 = load('C:/Users/Ryan/Documents/MATLAB/gratings-task-analysis/processed_data//MUA_GRATINGS_SUMMARY//M20170311-sessionInd8-muaAnalysisSummaryData-v11.mat')

Data = struct();

t = S1.spdfInfo.concatCueOnsetT;
newTStep = 0.01;
newTStart = t(1):newTStep:(t(end) - newTStep);

R = [S1.spdfInfo.meanSpdfInRFConcatAll'; S2.spdfInfo.meanSpdfInRFConcatAll'];
newR = nan(numel(newTStart), size(R, 2));

% binning truncates the end
for i = 1:numel(newTStart)
    tInd = t >= newTStart(i) & t < (newTStart(i) + newTStep);
    newR(i,:) = mean(R(tInd,:), 1);
end


Data(1).A = newR;
Data(1).times = (newTStart + newTStep / 2) * 1000; % bin centers in ms


R = [S1.spdfInfo.meanSpdfExRFConcatAll'; S2.spdfInfo.meanSpdfExRFConcatAll'];
newR = nan(numel(newTStart), size(R, 2));

% binning truncates the end
for i = 1:numel(newTStart)
    tInd = t >= newTStart(i) & t < (newTStart(i) + newTStep);
    newR(i,:) = mean(R(tInd,:), 1);
end

Data(2).A = newR;
Data(2).times = (newTStart + newTStep / 2) * 1000; % bin centers in ms

save('test.mat', 'Data');
