function ax = computeFanoFactorAI(firingRates, inRFLocs, exRFLocs, col, isSigUnit)

nUnits = numel(firingRates);
fanoFactorInRF = nan(nUnits, 1);
fanoFactorExRF = nan(nUnits, 1);
for i = 1:nUnits
    fanoFactorInRF(i) = computeFanoFactor(firingRates(i), inRFLocs(i));
    fanoFactorExRF(i) = computeFanoFactor(firingRates(i), exRFLocs(i));
end

maxFanoFactor = 5;
toRemove = fanoFactorInRF > maxFanoFactor | fanoFactorExRF > maxFanoFactor;
fanoFactorInRF(toRemove) = [];
fanoFactorExRF(toRemove) = [];
isSigUnit(toRemove) = [];

fanoFactorAI = (fanoFactorInRF - fanoFactorExRF) ./ (fanoFactorInRF + fanoFactorExRF);

mInRF = fanoFactorAI;
mExRF = zeros(size(fanoFactorAI));

mDiff = mInRF - mExRF;

meanMDiff = mean(mDiff);
medianMDiff = median(mDiff);
p = signrank(mDiff);
fprintf('\tAll: Mean diff = %0.3f, median diff = %0.3f, sign rank test p = %0.5f, N = %d (%d%% units selective)\n', ...
        meanMDiff, medianMDiff, p, numel(mDiff), round(sum(isSigUnit) / numel(isSigUnit) * 100));

%% plot parameters
histXBounds = [-0.8 0.8];%[-ceil(maxAbsDiffFR / histBinStep) ceil(maxAbsDiffFR / histBinStep)] * histBinStep;
histBinEdges = histXBounds(1):0.1:histXBounds(2);%histXBounds(1):histBinStep:histXBounds(2);

%% plot
f1 = figure_tr_inch(5, 4);

%% histogram of differences ventral pulvinar
ax = subaxis(1, 1, 1, 'MB', 0.2, 'MT', 0.05, 'ML', 0.15); 
hold on;
histH = histogram(mDiff, histBinEdges);
histH.FaceColor = col;
histH = histogram(mDiff(isSigUnit), histBinEdges);
histH.FaceColor = zeros(3, 1);

origYLim = ylim();
plot([0 0], origYLim, 'k', 'LineWidth', 2); 
xlim(histXBounds);
ylim(origYLim);
ylabel('Number of Units');
box off;
set(gca, 'FontSize', 18);
set(gca, 'LineWidth', 1);
set(gca, 'TickDir', 'out');

