function plotRasterByTime(spikeStats, eventInfo, windowName, varargin)

%% setup
% figure details
makeFig = 1;
spikeFile = [];
overridedefaults(who, varargin);

if makeFig
    figure;
end
hold on;

window = eventInfo.eventWindows(windowName);
data = spikeStats.locData{:};

%% plotting
rasterY = 0;
plotParams = {'Color', [0 0 0], 'MarkerSize', 1.5, 'MarkerFaceColor', [0 0 0]};

% plot diamond at each spike, one row (y-coord) per trial
for j = 1:numel(data)
    rasterY = rasterY + 1;
    if ~isempty(data(j).times)
        plot(data(j).times - window(1), rasterY*ones(1, length(data(j))),...
                'd', plotParams{:});
    end
end

rasterY = rasterY + 1; % 1 past the last raster row

%% formatting
hold on; 
plot([0 0], [0 rasterY], 'k');
xlim([-window(1) window(2)]);
ylim([0 rasterY]);
set(gca, 'YDir', 'reverse');
xlabel('Time (s)'); ylabel('Trial #');
titleH = title([spikeFile ' Raster (' windowName '), Sorted by Trial Time']);
set(titleH, 'Interpreter', 'None');
set(gca, 'Layer', 'top')
