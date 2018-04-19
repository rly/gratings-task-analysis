function quickImagePlotAllEvents5InARow(enterFixationR, ...
        cueOnsetR, arrayOnsetR, targetDimR, exitFixationR, ...
        enterFixationT, cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        groupName, isDiff, plotFileName)

% check size assert
nSessions = size(cueOnsetR, 1);
assert(all([size(arrayOnsetR, 1);
        size(targetDimR, 1);
        size(enterFixationR, 1); 
        size(exitFixationR, 1)] == nSessions))

%%
f = figure_tr_inch(19, 3.75); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% location params
spdfW = 0.17;
spdfH = 0.72;

spdfEnterFixationLeft = 0.04;
spdfCueOnsetLeft = spdfEnterFixationLeft + spdfW + 0.02;
spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW + 0.02;
spdfTargetDimLeft = spdfArrayOnsetLeft + spdfW + 0.02;
spdfExitFixationLeft = spdfTargetDimLeft + spdfW + 0.02;

btm = 0.16;

%% spdf for enter fixation
axEnterFixationSpdf = axes('Position', [spdfEnterFixationLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(enterFixationR, enterFixationT, 'Enter Fixation', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);
ylabel('Unit Number');

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(cueOnsetR, cueOnsetT, 'Cue Onset', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(arrayOnsetR, arrayOnsetT, 'Array Onset', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(targetDimR, targetDimT, 'Target Dimming', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);

%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
[hImage,hColorbar] = makeImagePlotOfPopulation(exitFixationR, exitFixationT, 'Exit Fixation', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);
set(hColorbar, 'Visible', 'on');

%% save
if ~isempty(plotFileName)
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end