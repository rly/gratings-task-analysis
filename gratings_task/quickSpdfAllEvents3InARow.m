function quickSpdfAllEvents3InARow(evokedSpikingSaveFileName, ...
        nLoc, plotFileName)

ES = load(evokedSpikingSaveFileName);

cols = lines(4);
inRFCol = cols(1,:);
exRFCol = cols(2,:);

%%
f = figure_tr_inch(18, 4); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

%% location params
spdfW = 0.285;
spdfH = 0.74;

spdfCueOnsetLeft = 0.06;
spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.03;
spdfTargetDimLeft = spdfArrayOnsetLeft + spdfW + 0.03;

btm = 0.2;

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

t = ES.cueOnsetT - ES.cueOnsetWindow(1);
xBounds = [-0.4 0.4];
hold on;
% TODO plot exRF first and inRF over it
for i = 1:nLoc
    if any(isnan(ES.cueOnsetSpdfByLoc(i,:)))
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
            ES.cueOnsetSpdfByLoc(i,:) + ES.cueOnsetSpdfErrByLoc(i,:), ...
            ES.cueOnsetSpdfByLoc(i,:) - ES.cueOnsetSpdfErrByLoc(i,:), ...
            col, col, 0.5);
    plot(t, ES.cueOnsetSpdfByLoc(i,:), ...
            'Color', col, 'LineWidth', 4);
end
origYLim = ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press

% shade in baseline and analysis windows
fillH = jbfill(ES.preCueBaselineWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.cueResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Cue Onset (s)');
ylabel('Estimated Spike Rate (Hz)');
set(gca, 'FontSize', 18);
% set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'XTick', -0.4:0.2:0.4);
set(gca, 'FontName', 'Calibri');
set(gca, 'LineWidth', 1);

% text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

t = ES.arrayOnsetT - ES.arrayOnsetWindow(1);
xBounds = [-0.4 0.4];
hold on;
for i = 1:nLoc
    if any(isnan(ES.arrayOnsetHoldSpdfByLoc(i,:))) % TODO show release trials
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
            ES.arrayOnsetHoldSpdfByLoc(i,:) + ES.arrayOnsetHoldSpdfErrByLoc(i,:), ...
            ES.arrayOnsetHoldSpdfByLoc(i,:) - ES.arrayOnsetHoldSpdfErrByLoc(i,:), ...
            col, col, 0.5);
    plot(t, ES.arrayOnsetHoldSpdfByLoc(i,:),  ...
            'Color', col, 'LineWidth', 4);
end
origYLim = ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

% shade in baseline and analysis windows
fillH = jbfill(ES.cueTargetDelayWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.arrayResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Array Onset (s)');
% ylabel('Estimated spike rate (Hz)');
set(gca, 'FontSize', 18);
% set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'XTick', -0.4:0.2:0.4);
set(gca, 'FontName', 'Calibri');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 1);

% text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

t = ES.targetDimT - ES.targetDimWindow(1);
xBounds = [-0.4 0.4];
hold on;
for i = 1:nLoc
    if any(isnan(ES.targetDimSpdfByLoc(i,:)))
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
            ES.targetDimSpdfByLoc(i,:) + ES.targetDimSpdfErrByLoc(i,:), ...
            ES.targetDimSpdfByLoc(i,:) - ES.targetDimSpdfErrByLoc(i,:), ...
            col, col, 0.5);
    plot(t, ES.targetDimSpdfByLoc(i,:),  ...
            'Color', col, 'LineWidth', 4);
end
origYLim = ylim();
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);
plot([0.28 0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min rt

% shade in baseline and analysis windows
fillH = jbfill(ES.targetDimDelayWindowOffset, [0 0], [1000 1000], ...
        [0.7 0.7 1], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
fillH = jbfill(ES.targetDimResponseWindowOffset, [0 0], [1000 1000], ...
        [0.3 1 0.3], ones(3, 1), 0.1);
uistack(fillH, 'bottom');
hold on;

xlim(xBounds);
ylim(origYLim);
xlabel('Time from Target Dimming (s)');
% ylabel('Estimated spike rate (Hz)');
set(gca, 'FontSize', 18);
% set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'XTick', -0.4:0.2:0.4);
set(gca, 'FontName', 'Calibri');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 1);


%% adjust all ylim
allYBounds = [ylim(axCueOnsetSpdf) ylim(axArrayOnsetSpdf) ylim(axTargetDimSpdf)];
yBounds = [min(allYBounds) max(1, max(allYBounds))];
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);

%% save
if ~isempty(plotFileName)
    export_fig(plotFileName, '-nocrop');
end
