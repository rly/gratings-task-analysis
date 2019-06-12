clear;

processedDataDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\PUL_SUA_GRATINGS_ALL';
fileNames = {'M20170608_PUL_43b-g2-g3-g4-evokedSpiking-v13.mat'};

saveFileDir = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\demo_figures\';
saveFileNames = {'M20170608_PUL_43b-example'};

inRFCol = [0.9 0 0];
exRFCol = [0 0 0.9];
nLoc = 4;

for fi = 1:numel(fileNames)
    ES = load(sprintf('%s/%s', processedDataDir, fileNames{fi}));

    %% create figure
    f = figure_tr_inch(6.5, 5); 
    clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% location params
    rasterW = 0.8;
    spdfW = rasterW;
    rasterH = 0.35;
    spdfH = 0.35;
    
    rasterCueOnsetLeft = 0.15;
    spdfCueOnsetLeft = rasterCueOnsetLeft;

    btm = 0.17;
    spdfBtm = btm;
    rasterBtm = spdfBtm + spdfH + 0.09;
    
    %% plot raster aligned to cue
    axes('Position', [rasterCueOnsetLeft rasterBtm rasterW rasterH]); 
    hold on;

    xBounds = [-0.3 0.3];
    window = ES.cueOnset.window;
    lineParams = {'LineWidth', 2};
    lineHeight = 5;
    rasterY = 0;
    for i = nLoc:-1:1
        data = ES.cueOnset.spikeTimesByLoc{i};
        if isempty(data)
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
        for k = 1:numel(data)
            rasterY = rasterY + 1;
            if ~isempty(data(k).times)
                adjTimes = data(k).times - window(1);
                for l = 1:numel(data(k).times)
                    % plot mini lines at each spike time
                    plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:}, 'Color', col);
                end
            end
        end
    end

    % plot event line marker
    yBounds = [0 rasterY+1];
    plot([0 0], yBounds, '-', 'Color', 0.3*ones(3, 1), 'LineWidth', 1);

    xlim(xBounds);
    ylim(yBounds);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTickLabel', []);
    ylabel('Trial Number');
    set(gca, 'FontSize', 18);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out');

    %% spdf for cue onset
    axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft spdfBtm spdfW spdfH]); 

    t = ES.cueOnset.t - ES.cueOnset.window(1);
    xBounds = [-0.3 0.3];
    yBounds = [0 63];
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
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 1);
    plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

    xlim(xBounds);
    ylim(yBounds);
    xlabel('Time from Cue Onset (s)');
    ylabel('Spike Rate (Hz)');
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out');

    %% save
    plotFileName = sprintf('%s/%s-cue.png', saveFileDir, saveFileNames{fi});
    if ~isempty(plotFileName)
        export_fig(plotFileName, '-nocrop');
    end
    
    %% create figure
    f = figure_tr_inch(6.5, 5); 
    clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %% location params
    rasterW = 0.8;
    spdfW = rasterW;
    rasterH = 0.35;
    spdfH = 0.35;
    
    rasterArrayOnsetLeft = 0.15;
    spdfArrayOnsetLeft = rasterArrayOnsetLeft;

    btm = 0.17;
    spdfBtm = btm;
    rasterBtm = spdfBtm + spdfH + 0.09;
    
    %% plot raster aligned to array - NOT hold balanced
    axes('Position', [rasterArrayOnsetLeft rasterBtm rasterW rasterH]); 
    hold on;

    xBounds = [-0.3 0.3];
    window = ES.arrayOnset.window;
    lineParams = {'LineWidth', 2};
    lineHeight = 5;
    rasterY = 0;
    for i = nLoc:-1:1
        data = ES.arrayOnset.spikeTimesByLoc{i};
        if isempty(data)
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
        for k = 1:numel(data)
            rasterY = rasterY + 1;
            if ~isempty(data(k).times)
                adjTimes = data(k).times - window(1);
                for l = 1:numel(data(k).times)
                    % plot mini lines at each spike time
                    plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:}, 'Color', col);
                end
            end
        end
    end

    % plot event line marker
    yBounds = [0 rasterY+1];
    plot([0 0], yBounds, '-', 'Color', 0.3*ones(3, 1), 'LineWidth', 1);

    xlim(xBounds);
    ylim(yBounds);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTickLabel', []);
    ylabel('Trial Number');
    set(gca, 'FontSize', 18);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out');
    
    %% spdf for array onset - NOT hold balanced
    axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft spdfBtm spdfW spdfH]); 

    t = ES.arrayOnset.t - ES.arrayOnset.window(1);
    xBounds = [-0.3 0.3];
    hold on;
    for i = 1:nLoc
        if any(isnan(ES.arrayOnset.spdfByLoc(i,:)))
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
                ES.arrayOnset.spdfByLoc(i,:) + ES.arrayOnset.spdfErrByLoc(i,:), ...
                ES.arrayOnset.spdfByLoc(i,:) - ES.arrayOnset.spdfErrByLoc(i,:), ...
                col, col, 0.5);
        plot(t, ES.arrayOnset.spdfByLoc(i,:),  ...
                'Color', col, 'LineWidth', 4);
    end
    origYLim = ylim();
    plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 1);

    xlim(xBounds);
    ylim(origYLim);
    xlabel('Time from Array Onset (s)');
    ylabel('Spike Rate (Hz)');
    set(gca, 'FontSize', 18);
    set(gca, 'XTick', -0.3:0.1:0.3);
    set(gca, 'LineWidth', 1.5);
    set(gca, 'TickDir', 'out');
    
    %% save
    plotFileName = sprintf('%s/%s-arr.png', saveFileDir, saveFileNames{fi});
    if ~isempty(plotFileName)
        export_fig(plotFileName, '-nocrop');
    end
    
    
end