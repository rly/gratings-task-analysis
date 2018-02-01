function quickSpdfAllEvents5InARow(evokedSpikingSaveFileName, ...
        nLoc, plotFileName)

ES = load(evokedSpikingSaveFileName);

cols = lines(4);
inRFCol = [0.9 0 0];
exRFCol = [0 0 0.9];

%%
f = figure_tr_inch(19, 3.75); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% make main title
axBig = axes('Position', [0.04 0.045 0.92 0.91], 'Visible', 'off');
set(get(axBig, 'Title'), 'Visible', 'on')

%% location params
spdfW = 0.17;
spdfH = 0.74;

spdfEnterFixationLeft = 0.05;
spdfCueOnsetLeft = spdfEnterFixationLeft + spdfW + 0.02;
spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.02;
spdfTargetDimLeft = spdfArrayOnsetLeft + spdfW + 0.02;
spdfExitFixationLeft = spdfTargetDimLeft + spdfW + 0.02;

btm = 0.2;


%% spdf for enter fixation
axEnterFixationSpdf = axes('Position', [spdfEnterFixationLeft btm spdfW spdfH]); 

t = ES.enterFixationT - ES.enterFixationWindow(1);
xBounds = [-0.4 0.4];
hold on;

plot([0.35 0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

for i = [ES.exRFLoc ES.inRFLoc]
    if any(isnan(ES.enterFixationSpdfByLoc(i,:)))
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
            ES.enterFixationSpdfByLoc(i,:) + ES.enterFixationSpdfErrByLoc(i,:), ...
            ES.enterFixationSpdfByLoc(i,:) - ES.enterFixationSpdfErrByLoc(i,:), ...
            col, col, 0.5);
    plot(t, ES.enterFixationSpdfByLoc(i,:), ...
            'Color', col, 'LineWidth', 4);
end

xlim(xBounds);
% xlabel('Time from Begin Fixation (s)');
% ylabel('Estimated Spike Rate (Hz)');
set(gca, 'FontSize', 26);
set(gca, 'XTick', -0.25:0.25:0.25);
set(gca, 'FontName', 'Calibri');
set(gca, 'LineWidth', 2);

% text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

t = ES.cueOnsetT - ES.cueOnsetWindow(1);
xBounds = [-0.4 0.4];
hold on;

plot([-0.35 -0.35], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min enter fixation/lever press
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

for i = [ES.exRFLoc ES.inRFLoc]
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

xlim(xBounds);
% xlabel('Time from Cue Onset (s)');
% ylabel('Estimated Spike Rate (Hz)');
set(gca, 'FontSize', 26);
set(gca, 'XTick', -0.25:0.25:0.25);
set(gca, 'FontName', 'Calibri');
% set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);

% text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

t = ES.arrayOnsetT - ES.arrayOnsetWindow(1);
xBounds = [-0.4 0.4];
hold on;

plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

for i = [ES.exRFLoc ES.inRFLoc]
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

xlim(xBounds);
% xlabel('Time from Array Onset (s)');
% ylabel('Estimated spike rate (Hz)');
set(gca, 'FontSize', 26);
set(gca, 'XTick', -0.25:0.25:0.25);
set(gca, 'FontName', 'Calibri');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);

% text(0.9, 1, legendEntry, 'VerticalAlignment', 'top', 'Units', 'normalized', 'FontSize', 10);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

t = ES.targetDimT - ES.targetDimWindow(1);
xBounds = [-0.4 0.4];
hold on;

plot([0.28 0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1)); % line at min rt
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

for i = [ES.exRFLoc ES.inRFLoc]
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

xlim(xBounds);
% xlabel('Time from Target Dimming (s)');
% ylabel('Estimated spike rate (Hz)');
set(gca, 'FontSize', 26);
set(gca, 'XTick', -0.25:0.25:0.25);
set(gca, 'FontName', 'Calibri');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);


%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

t = ES.exitFixationT - ES.exitFixationWindow(1);
xBounds = [-0.4 0.4];
hold on;

plot([-0.28 -0.28], [-1000 1000], '--', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2); % line at min rt
plot([0 0], [-1000 1000], '-', 'Color', 0.3 * ones(3, 1), 'LineWidth', 2);

for i = [ES.exRFLoc ES.inRFLoc]
    if any(isnan(ES.exitFixationSpdfByLoc(i,:)))
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
            ES.exitFixationSpdfByLoc(i,:) + ES.exitFixationSpdfErrByLoc(i,:), ...
            ES.exitFixationSpdfByLoc(i,:) - ES.exitFixationSpdfErrByLoc(i,:), ...
            col, col, 0.5);
    plot(t, ES.exitFixationSpdfByLoc(i,:),  ...
            'Color', col, 'LineWidth', 4);
end

xlim(xBounds);
% xlabel('Time from End Fixation (s)');
% ylabel('Estimated spike rate (Hz)');
set(gca, 'FontSize', 26);
set(gca, 'XTick', -0.25:0.25:0.25);
set(gca, 'FontName', 'Calibri');
set(gca, 'YTickLabel', []);
set(gca, 'LineWidth', 2);



%% adjust all ylim
allYBounds = [ES.cueOnsetSpdfByLoc(i,:) + ES.cueOnsetSpdfErrByLoc(i,:), ...
        ES.arrayOnsetHoldSpdfByLoc(i,:) + ES.arrayOnsetHoldSpdfErrByLoc(i,:), ...
        ES.targetDimSpdfByLoc(i,:) + ES.targetDimSpdfErrByLoc(i,:), ...
        ES.exitFixationSpdfByLoc(i,:) + ES.exitFixationSpdfErrByLoc(i,:)];
yBounds = [0 max(allYBounds)];
ylim(axEnterFixationSpdf, yBounds);
ylim(axCueOnsetSpdf, yBounds);
ylim(axArrayOnsetSpdf, yBounds);
ylim(axTargetDimSpdf, yBounds);
ylim(axExitFixationSpdf, yBounds);

%% save
if ~isempty(plotFileName)
    export_fig(plotFileName, '-nocrop');
end
