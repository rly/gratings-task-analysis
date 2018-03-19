function plotRTDistribution(rtRelInRF, rtRelExRF, rtHoldInRF, rtHoldExRF, ...
        checkRTStatAlpha, sessionName, isZeroDistractors, plotFileName)
%%
f = figure_tr_inch(13, 7.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('RT Distribution: %s', sessionName);
if isZeroDistractors
    modTitle = [modTitle ' (0 Distractors)'];
end
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
scatterW = 0.44;
scatterH = 0.385;

col1Left = 0.05;
col2Left = col1Left + scatterW + 0.05;

row2Btm = 0.07;
row1Btm = row2Btm + scatterH + 0.08;

%%
binEdges = 0.3:0.05:0.8;
xBounds = binEdges([1 end]);
textParams = {'FontSize', 10, 'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top'};
plotHs = nan(4, 1);
cols = lines(2);
inRFCol = cols(1,:);
exRFCol = cols(2,:);

% there should not be any difference between RTs on InRF and
% ExRF trials. otherwise there is significant spatial bias
% in this session.
p = ranksum(rtRelInRF, rtRelExRF);
if p < checkRTStatAlpha
    set(gcf, 'Color', 'red');
end
p = ranksum(rtHoldInRF, rtHoldExRF);
if p < checkRTStatAlpha
    set(gcf, 'Color', 'red');
end

%%
plotHs(1) = axes('Position', [col1Left row1Btm scatterW scatterH]); 
hold on;
hist1 = histogram(rtRelInRF, binEdges);
hist1.FaceColor = inRFCol;
xlim(xBounds);
title(sprintf('Release InRF (N=%d)', numel(rtRelInRF)));
text(1.0, 1.0, sprintf('Median = %0.3f s', median(rtRelInRF)), 'Color', inRFCol, textParams{:});
ylabel('Number of Trials');

%%
plotHs(2) = axes('Position', [col1Left row2Btm scatterW scatterH]); 
hold on;
hist2 = histogram(rtRelExRF, binEdges);
hist2.FaceColor = exRFCol;
xlim(xBounds);
title(sprintf('Release ExRF (N=%d)', numel(rtRelExRF)));
text(1.0, 1.0, sprintf('Median = %0.3f s', median(rtRelExRF)), 'Color', exRFCol, textParams{:});
xlabel('Response Time (s)');
ylabel('Number of Trials');

%%
plotHs(3) = axes('Position', [col2Left row1Btm scatterW scatterH]); 
hold on;
hist3 = histogram(rtHoldInRF, binEdges);
hist3.FaceColor = inRFCol;
xlim(xBounds);
title(sprintf('Hold InRF (N=%d)', numel(rtHoldInRF)));
text(1.0, 1.0, sprintf('Median = %0.3f s', median(rtHoldInRF)), 'Color', inRFCol, textParams{:});

%%
plotHs(4) = axes('Position', [col2Left row2Btm scatterW scatterH]); 
hold on;
hist4 = histogram(rtHoldExRF, binEdges);
hist4.FaceColor = exRFCol;
xlim(xBounds);
title(sprintf('Hold ExRF (N=%d)', numel(rtHoldExRF)));
text(1.0, 1.0, sprintf('Median = %0.3f s', median(rtHoldExRF)), 'Color', exRFCol, textParams{:});
xlabel('Response Time (s)');

%% set all y bounds the same
allYBounds = arrayfun(@(x) ylim(x), plotHs, 'UniformOutput', false);
allYBounds = [allYBounds{:}];
yBounds = [min(allYBounds) max(allYBounds)];
arrayfun(@(x) plot(x, [0 0], yBounds, 'Color', 0.3*ones(3, 1)), plotHs);
arrayfun(@(x) ylim(x, yBounds), plotHs);

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
