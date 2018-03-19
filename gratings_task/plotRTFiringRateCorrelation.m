function plotRTFiringRateCorrelation(cueTargetDelayLongRelInRFRate, ...
                    cueTargetDelayLongRelExRFRate, ...
                    cueTargetDelayLongHoldInRFRate, ...
                    targetDimDelayLongHoldInRFRate, ...
                    cueTargetDelayLongHoldExRFRate, ...
                    targetDimDelayLongHoldExRFRate, ...
                    rtRelInRF, rtRelExRF, rtHoldInRF, rtHoldExRF, ...
                    rtFiringRateInfo, unitName, isZeroDistractors, plotFileName)
%%
f = figure_tr_inch(13, 7.5); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

modTitle = sprintf('RT vs Firing Rate: %s', unitName);
if isZeroDistractors
    modTitle = [modTitle ' (0 Distractors)'];
end
titleParams = {'Interpreter', 'None', 'FontWeight', 'bold'};
title(modTitle, 'FontSize', 14, titleParams{:});

%% location params
scatterW = 0.275;
scatterH = 0.385;

col1Left = 0.05;
col2Left = col1Left + scatterW + 0.05;
col3Left = col2Left + scatterW + 0.05;

row2Btm = 0.07;
row1Btm = row2Btm + scatterH + 0.08;

%%
xBounds = [0 max([cueTargetDelayLongRelInRFRate; ...
                    cueTargetDelayLongRelExRFRate; ...
                    cueTargetDelayLongHoldInRFRate; ...
                    targetDimDelayLongHoldInRFRate; ...
                    cueTargetDelayLongHoldExRFRate; ...
                    targetDimDelayLongHoldExRFRate])];
                
yBounds = [0.3 0.8];
plotParams = {'MarkerSize', 15};
            
textParamsNormal = {'FontSize', 10, 'Units', 'normalized'};
textParamsBold = {'FontSize', 10, 'Units', 'normalized', 'FontWeight', 'bold'};
corrStatAlpha = 0.05;

%%
axes('Position', [col1Left row1Btm scatterW scatterH]); 
hold on;
plot(cueTargetDelayLongRelInRFRate, rtRelInRF, '.', plotParams{:});
LM = fitlm(cueTargetDelayLongRelInRFRate, rtRelInRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Release InRF: RT vs Cue-Target Delay Firing (N=%d)', numel(rtRelInRF)));
if rtFiringRateInfo.spearmanCorrCoefPValRelInRFCTDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefRelInRFCTDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValRelInRFCTDelayRT), textParams{:});
ylabel('Response Time (s)');

%%
axes('Position', [col1Left row2Btm scatterW scatterH]); 
hold on;
plot(cueTargetDelayLongRelExRFRate, rtRelExRF, '.', plotParams{:});
LM = fitlm(cueTargetDelayLongRelExRFRate, rtRelExRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Release ExRF: RT vs Cue-Target Delay Firing (N=%d)', numel(rtRelExRF)));
if rtFiringRateInfo.spearmanCorrCoefPValRelExRFCTDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefRelExRFCTDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValRelExRFCTDelayRT), textParams{:});
xlabel('Cue-Target Delay Firing Rate (Hz)');
ylabel('Response Time (s)');

%%
axes('Position', [col2Left row1Btm scatterW scatterH]); 
hold on;
plot(cueTargetDelayLongHoldInRFRate, rtHoldInRF, '.', plotParams{:});
LM = fitlm(cueTargetDelayLongHoldInRFRate, rtHoldInRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Hold InRF: RT vs Cue-Target Delay Firing (N=%d)', numel(rtHoldInRF)));
if rtFiringRateInfo.spearmanCorrCoefPValHoldInRFCTDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefHoldInRFCTDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValHoldInRFCTDelayRT), textParams{:});

%%
axes('Position', [col2Left row2Btm scatterW scatterH]); 
hold on;
plot(cueTargetDelayLongHoldExRFRate, rtHoldExRF, '.', plotParams{:});
LM = fitlm(cueTargetDelayLongHoldExRFRate, rtHoldExRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Hold ExRF: RT vs Cue-Target Delay Firing (N=%d)', numel(rtHoldExRF)));
if rtFiringRateInfo.spearmanCorrCoefPValHoldExRFCTDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefHoldExRFCTDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValHoldExRFCTDelayRT), textParams{:});
xlabel('Cue-Target Delay Firing Rate (Hz)');

%%
axes('Position', [col3Left row1Btm scatterW scatterH]); 
hold on;
plot(targetDimDelayLongHoldInRFRate, rtHoldInRF, '.', plotParams{:});
LM = fitlm(targetDimDelayLongHoldInRFRate, rtHoldInRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Hold InRF: RT vs Target-Dim Delay Firing (N=%d)', numel(rtHoldInRF)));
if rtFiringRateInfo.spearmanCorrCoefPValHoldInRFTDDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefHoldInRFTDDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValHoldInRFTDDelayRT), textParams{:});

%%
axes('Position', [col3Left row2Btm scatterW scatterH]); 
hold on;
plot(targetDimDelayLongHoldExRFRate, rtHoldExRF, '.', plotParams{:});
LM = fitlm(targetDimDelayLongHoldExRFRate, rtHoldExRF);
predY = predict(LM, xBounds');
plot(xBounds, predY, 'LineWidth', 2);
xlim(xBounds);
ylim(yBounds);
title(sprintf('Hold ExRF: RT vs Target-Dim Delay Firing (N=%d)', numel(rtHoldExRF)));
if rtFiringRateInfo.spearmanCorrCoefPValHoldExRFTDDelayRT < corrStatAlpha
    textParams = textParamsBold;
else
    textParams = textParamsNormal;
end
text(0.75, 0.1, sprintf('Spear R = %0.3f', rtFiringRateInfo.spearmanCorrCoefHoldExRFTDDelayRT), textParams{:});
text(0.75, 0.05, sprintf('p = %0.3f', rtFiringRateInfo.spearmanCorrCoefPValHoldExRFTDDelayRT), textParams{:});
xlabel('Target-Dim Delay Firing Rate (Hz)');

%% save
if ~isempty(plotFileName)
    fprintf('\tSaving figure to file %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end
