% call this script from grcMuaPlots.m
outputDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\GRC';

subdivisions = {'dPul2', 'vPul2'};
for j = 1:numel(subdivisions)
    subdivision = subdivisions{j};
    yBounds = [-0.2 0.6];
    isShowLabels = 1;
    
    if strcmp(subdivision, 'dPul')
        condition = isInDPulvinar & isSignificantCueResponseInc;
    elseif strcmp(subdivision, 'vPul')
        condition = isInVPulvinar  & isSignificantCueResponseInc;
    elseif strcmp(subdivision, 'dPul2')
        condition = isInDPulvinar & isSignificantCueResponseInc & isSignificantSelectivityCueResponse;
    elseif strcmp(subdivision, 'vPul2')
        condition = isInVPulvinar  & isSignificantCueResponseInc & isSignificantSelectivityCueResponse;
    end
    cueOnsetSpdfInRFNorm = (spdfInfo.cueOnsetSpdfInRFNorm(condition,:));
    cueOnsetSpdfExRFNorm = (spdfInfo.cueOnsetSpdfExRFNorm(condition,:));
    arrayOnsetHoldSpdfInRFNorm = (spdfInfo.arrayOnsetHoldSpdfInRFNorm(condition,:));
    arrayOnsetHoldSpdfExRFNorm = (spdfInfo.arrayOnsetHoldSpdfExRFNorm(condition,:));
    targetDimSpdfInRFNorm = (spdfInfo.targetDimSpdfInRFNorm(condition,:));
    targetDimSpdfExRFNorm = (spdfInfo.targetDimSpdfExRFNorm(condition,:));
    exitFixationSpdfInRFNorm = (spdfInfo.exitFixationSpdfInRFNorm(condition,:));
    exitFixationSpdfExRFNorm = (spdfInfo.exitFixationSpdfExRFNorm(condition,:));
    
    plotFileName = sprintf('%s/allSessions-%s-meanSpdfsCueArrayExitFix-v%d.png', outputDir, subdivision, v);

    fprintf('%s: %d units\n', subdivision, sum(condition));
    
    % check size assert
    nSessions = size(cueOnsetSpdfInRFNorm, 1);

    meanCueOnsetSpdfInRFNorm = mean(cueOnsetSpdfInRFNorm);
    meanCueOnsetSpdfExRFNorm = mean(cueOnsetSpdfExRFNorm);
    meanArrayOnsetHoldSpdfInRFNorm = mean(arrayOnsetHoldSpdfInRFNorm);
    meanArrayOnsetHoldSpdfExRFNorm = mean(arrayOnsetHoldSpdfExRFNorm);
    meanExitFixationSpdfInRFNorm = mean(exitFixationSpdfInRFNorm);
    meanExitFixationSpdfExRFNorm = mean(exitFixationSpdfExRFNorm);

    seCueOnsetSpdfInRFNorm = std(cueOnsetSpdfInRFNorm)/sqrt(nSessions);
    seCueOnsetSpdfExRFNorm = std(cueOnsetSpdfExRFNorm)/sqrt(nSessions);
    seArrayOnsetHoldSpdfInRFNorm = std(arrayOnsetHoldSpdfInRFNorm)/sqrt(nSessions);
    seArrayOnsetHoldSpdfExRFNorm = std(arrayOnsetHoldSpdfExRFNorm)/sqrt(nSessions);
    seExitFixationSpdfInRFNorm = std(exitFixationSpdfInRFNorm)/sqrt(nSessions);
    seExitFixationSpdfExRFNorm = std(exitFixationSpdfExRFNorm)/sqrt(nSessions);

    inRFCol = [0.9 0 0];
    exRFCol = [0 0 0.9];

    %%
    f = figure_tr_inch(17, 5); 
    clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% location params
    spdfW = 0.285;
    spdfH = 0.8;

    spdfCueOnsetLeft = 0.05;
    spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.04;
    spdfExitFixationLeft = spdfArrayOnsetLeft + spdfW + 0.04;

    btm = 0.15;

    %% spdf for cue onset
    axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

    xBounds = [-0.3 0.3];
    hold on;
    plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min enter fixation/lever press
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
    plot(cueOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

    jbfill(cueOnsetT, ...
            meanCueOnsetSpdfExRFNorm + seCueOnsetSpdfExRFNorm, ...
            meanCueOnsetSpdfExRFNorm - seCueOnsetSpdfExRFNorm, ...
            exRFCol, exRFCol, 0.5);
    plot(cueOnsetT, meanCueOnsetSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);
    jbfill(cueOnsetT, ...
            meanCueOnsetSpdfInRFNorm + seCueOnsetSpdfInRFNorm, ...
            meanCueOnsetSpdfInRFNorm - seCueOnsetSpdfInRFNorm, ...
            inRFCol, inRFCol, 0.5);
    plot(cueOnsetT, meanCueOnsetSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

    xlim(xBounds);
    if isShowLabels
        xlabel('Time from Cue Onset (s)');
        ylabel('Mean Normalized Firing Rate');
    end
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'LineWidth', 1);

    %% spdf for array onset
    axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

    xBounds = [-0.3 0.3];
    hold on;

    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
    plot(arrayOnsetT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

    jbfill(arrayOnsetT, ...
            meanArrayOnsetHoldSpdfExRFNorm + seArrayOnsetHoldSpdfExRFNorm, ...
            meanArrayOnsetHoldSpdfExRFNorm - seArrayOnsetHoldSpdfExRFNorm, ...
            exRFCol, exRFCol, 0.5);
    plot(arrayOnsetT, meanArrayOnsetHoldSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);

    jbfill(arrayOnsetT, ...
            meanArrayOnsetHoldSpdfInRFNorm + seArrayOnsetHoldSpdfInRFNorm, ...
            meanArrayOnsetHoldSpdfInRFNorm - seArrayOnsetHoldSpdfInRFNorm, ...
            inRFCol, inRFCol, 0.5);
    plot(arrayOnsetT, meanArrayOnsetHoldSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

    xlim(xBounds);
    if isShowLabels
        xlabel('Time from Array Onset (s)');
    end
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'YTickLabel', []);
    set(gca, 'LineWidth', 1);

    %% spdf for exit fixation
    axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

    xBounds = [-0.3 0.3];
    hold on;

    plot([-0.28 -0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min rt
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
    plot(exitFixationT([1 end]), [0 0], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

    jbfill(exitFixationT, ...
            meanExitFixationSpdfExRFNorm + seExitFixationSpdfExRFNorm, ...
            meanExitFixationSpdfExRFNorm - seExitFixationSpdfExRFNorm, ...
            exRFCol, exRFCol, 0.5);
    plot(exitFixationT, meanExitFixationSpdfExRFNorm, 'Color', exRFCol, 'LineWidth', 4);

    jbfill(exitFixationT, ...
            meanExitFixationSpdfInRFNorm + seExitFixationSpdfInRFNorm, ...
            meanExitFixationSpdfInRFNorm - seExitFixationSpdfInRFNorm, ...
            inRFCol, inRFCol, 0.5);
    plot(exitFixationT, meanExitFixationSpdfInRFNorm, 'Color', inRFCol, 'LineWidth', 4);

    xlim(xBounds);
    if isShowLabels
        xlabel('Time from End Fixation (s)');
    end
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'YTickLabel', []);
    set(gca, 'LineWidth', 1);

    %% adjust all ylim
    allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf) ylim(axExitFixationSpdf)];
    ylim(axCueOnsetSpdf, yBounds);
    ylim(axArrayOnsetSpdf, yBounds);
    ylim(axExitFixationSpdf, yBounds);

    %% save
    if ~isempty(plotFileName)
        export_fig(plotFileName, '-nocrop');
    end
end