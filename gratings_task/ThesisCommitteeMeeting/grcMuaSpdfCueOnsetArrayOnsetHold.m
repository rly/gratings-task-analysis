clear;

processedDataDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\PUL_MUA_GRATINGS_ALL';
fileNames = {'M20170311_PUL_48M-g3-g4-g6-g7-g8-g9-evokedSpiking-v12.mat', ...
        'M20170331_PUL_46M-g3-g4-g5-g6-evokedSpiking-v12.mat'};

saveFileDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\GRC';
saveFileNames = {'M20170311_PUL_48M-example.png', ...
        'M20170331_PUL_46M-example.png'};

inRFCol = [0.9 0 0];
exRFCol = [0 0 0.9];
nLoc = 4;

for fi = 1:numel(fileNames)
    ES = load(sprintf('%s/%s', processedDataDir, fileNames{fi}));

    %% create figure
    f = figure_tr_inch(10, 5); 
    clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% location params
    spdfW = 0.40;
    spdfH = 0.74;

    spdfCueOnsetLeft = 0.08;
    spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.05;

    btm = 0.15;

    %% spdf for cue onset
    axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

    t = ES.cueOnset.t - ES.cueOnset.window(1);
    xBounds = [-0.3 0.3];
    hold on;
    % TODO plot exRF first and inRF over it
    for i = 1:nLoc
        if any(isnan(ES.cueOnset.spdfByLoc(i,:)))
            continue;
        end
        if i == ES.inRFLoc
            col = inRFCol;
        elseif i == ES.exRFLoc
            col = exRFCol;
        else
            col = 0.3*ones(3, 1);
            continue; % just two lines
        end
        jbfill(t, ...
                ES.cueOnset.spdfByLoc(i,:) + ES.cueOnset.spdfErrByLoc(i,:), ...
                ES.cueOnset.spdfByLoc(i,:) - ES.cueOnset.spdfErrByLoc(i,:), ...
                col, col, 0.5);
        plot(t, ES.cueOnset.spdfByLoc(i,:), ...
                'Color', col, 'LineWidth', 4);
    end
    origYLim = ylim();
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
    plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

%     % shade in baseline and analysis windows
%     fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
%             [0.7 0.7 1], ones(3, 1), 0.1);
%     uistack(fillH, 'bottom');
%     fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
%             [0.3 1 0.3], ones(3, 1), 0.1);
%     uistack(fillH, 'bottom');
%     hold on;

    xlim(xBounds);
    ylim(origYLim);
    xlabel('Time from Cue Onset (s)');
    ylabel('Estimated Spike Rate (Hz)');
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'LineWidth', 1);

    %% spdf for array onset
    axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

    t = ES.arrayOnset.t - ES.arrayOnset.window(1);
    xBounds = [-0.3 0.3];
    hold on;
    for i = 1:nLoc
        if any(isnan(ES.arrayOnsetHoldBal.spdfByLoc(i,:)))
            continue;
        end
        if i == ES.inRFLoc
            col = inRFCol;
        elseif i == ES.exRFLoc
            col = exRFCol;
        else
            col = 0.3*ones(3, 1);
            continue; % just two lines
        end
        jbfill(t, ...
                ES.arrayOnsetHoldBal.spdfByLoc(i,:) + ES.arrayOnsetHoldBal.spdfErrByLoc(i,:), ...
                ES.arrayOnsetHoldBal.spdfByLoc(i,:) - ES.arrayOnsetHoldBal.spdfErrByLoc(i,:), ...
                col, col, 0.5);
        plot(t, ES.arrayOnsetHoldBal.spdfByLoc(i,:),  ...
                'Color', col, 'LineWidth', 4);
    end
    origYLim = ylim();
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

    % shade in baseline and analysis windows
%     fillH = jbfill(ES.cueTargetDelayWindowOffset, [0 0], [1000 1000], ...
%             [0.7 0.7 1], ones(3, 1), 0.1);
%     uistack(fillH, 'bottom');
%     fillH = jbfill(ES.arrayResponseWindowOffset, [0 0], [1000 1000], ...
%             [0.3 1 0.3], ones(3, 1), 0.1);
%     uistack(fillH, 'bottom');
%     hold on;

    xlim(xBounds);
    ylim(origYLim);
    xlabel('Time from Array Onset (s)');
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'YTickLabel', []);
    set(gca, 'LineWidth', 1);

    %% adjust all ylim
    allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf)];
    yBounds = [min(allYBounds) max(1, max(allYBounds))];
    ylim(axCueOnsetSpdf, yBounds);
    ylim(axArrayOnsetSpdf, yBounds);

    %% save
    plotFileName = sprintf('%s/%s', saveFileDir, saveFileNames{fi});
    if ~isempty(plotFileName)
        export_fig(plotFileName, '-nocrop');
    end
end