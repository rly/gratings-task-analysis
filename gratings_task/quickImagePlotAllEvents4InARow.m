function quickImagePlotAllEvents4InARow(...
        cueOnsetR, arrayOnsetR, targetDimR, exitFixationR, ...
        cueOnsetT, arrayOnsetT, targetDimT, exitFixationT, ...
        groupName, isDiff, plotFileName)

% check size assert
nSessions = size(cueOnsetR, 1);
assert(all([size(arrayOnsetR, 1);
        size(targetDimR, 1);
        size(exitFixationR, 1)] == nSessions))

%%
f = figure_tr_inch(16, 3.75); clf;
set(gcf, 'Color', 'white');
set(gcf, 'renderer', 'painters');

%% location params
spdfW = 0.24;
spdfH = 0.72;

spdfCueOnsetLeft = 0.05;
spdfArrayOnsetLeft = spdfCueOnsetLeft + spdfW - 0.015;
spdfTargetDimLeft = spdfArrayOnsetLeft + spdfW- 0.015;
spdfExitFixationLeft = spdfTargetDimLeft + spdfW- 0.015;

btm = 0.16;

%% spdf for cue onset
axCueOnsetSpdf = axes('Position', [spdfCueOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(cueOnsetR, cueOnsetT, 'Cue Onset', isDiff);
ylabel('Unit Number');

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);

%% spdf for array onset
axArrayOnsetSpdf = axes('Position', [spdfArrayOnsetLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(arrayOnsetR, arrayOnsetT, 'Array Onset', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'YTick', []);

%% spdf for target dim
axTargetDimSpdf = axes('Position', [spdfTargetDimLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
makeImagePlotOfPopulation(targetDimR, targetDimT, 'Target Dimming', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);
set(gca, 'YTick', []);

%% spdf for exit fixation
axExitFixationSpdf = axes('Position', [spdfExitFixationLeft btm spdfW spdfH]); 

xBounds = [-0.4 0.4];
[hImage,hColorbar] = makeImagePlotOfPopulation(exitFixationR, exitFixationT, 'Break Fixation', isDiff);

xlim(xBounds);
set(gca, 'XTick', -0.5:0.25:0.5);
set(hColorbar, 'Visible', 'on');
set(gca, 'YTick', []);

%% save
if ~isempty(plotFileName)
    fprintf('Saving to %s...\n', plotFileName);
    export_fig(plotFileName, '-nocrop');
end