%% pul units w/ sig pos cue resp at contralateral (P3) 
precondition = isInPulvinar & isSignificantCueResponseInc & isInRFP3;
fprintf('Units in the pulvinar with significantly increased cue response compared to baseline and InRF P3: %d\n', sum(precondition));

suaAllPath = 'C:\Users\Ryan\Documents\MATLAB\gratings-task-analysis\processed_data\PUL_SUA_GRATINGS_ALL';

suaFiles = dir(suaAllPath);
esFilePathSub = esFileNames(precondition);

figLeft = 0.12;
figW = 0.85;
rasterBtm = 0.3;
rasterH = 0.68;
spdfBtm = 0.06;
spdfH = 0.2;

for i = 1:numel(esFilePathSub)
    fprintf('Processing file %d/%d = %d%%\n', i, numel(esFilePathSub), round(100*i/numel(esFilePathSub)));
    [~,esFileNameNoExt] = fileparts(esFilePathSub{i});
    suaPath = dir(sprintf('%s/%s.mat', suaAllPath, esFileNameNoExt));
    assert(numel(suaPath) == 1);
    ES = load([suaPath(1).folder '/' suaPath(1).name]);
    
    event = ES.enterFixation;
    eventName = 'Enter Fixation';
    xBounds = [-0.4 0.4];
    
    %%
    figure_tr_inch(5, 10);
    clf;
    set(gcf, 'Color', 'white');
    set(gcf, 'renderer', 'painters');

    %%
    axes('Position', [figLeft rasterBtm figW rasterH]); 
    hold on;

    window = event.window;
    data = event.spikeTimes;
    rasterY = 0;
    lineParams = {'Color', [0 0 0], 'LineWidth', 1};
    lineHeight = 5;
    for k = 1:numel(data)
        rasterY = rasterY + 1;
        if ~isempty(data(k).times)
            adjTimes = data(k).times - window(1);
            for l = 1:numel(data(k).times)
                % plot mini lines at each spike time
                plot([adjTimes(l) adjTimes(l)], rasterY + [-lineHeight lineHeight]/2, lineParams{:});
            end
        end
    end

    % plot event line marker
    yBounds = [0 rasterY+1];
    plot([0 0], yBounds, '-', 'Color', 0.5*ones(3, 1));

    xlim(xBounds);
    ylim(yBounds);
    set(gca, 'YDir', 'reverse');
    ylabel('Trial Number');
    
    %%
    axes('Position', [figLeft spdfBtm figW spdfH]); 

    t = event.t - event.window(1);
    hold on;
    legendEntry = cell(nLoc, 1);
    for j = 1:nLoc
        if any(isnan(event.spdfByLoc(j,:)))
            continue;
        end
        plot(t, event.spdfByLoc(j,:), ...
                'Color', cols(j,:), 'LineWidth', 2);

    end
    origYLim = ylim();
    plot([0 0], [0 1000], '-', 'Color', 0.3 * ones(3, 1));

    xlim(xBounds);
    ylim(origYLim);
    xlabel(['Time from ' eventName ' (s)']);
    ylabel('Estimated Spike Rate (Hz)');
    
    %% 
    drawnow;
end
