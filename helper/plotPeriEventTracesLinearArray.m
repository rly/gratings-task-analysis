function plotPeriEventTracesLinearArray(fPlot, t, averageResponses, varargin)

yOffset = 1;
yScale = 10;
overridedefaults(who, varargin);

nChannels = size(averageResponses, 1);

figure(fPlot);
hold on;
for j = 1:nChannels
    plot(t, yScale*averageResponses(j,:) - j*yOffset);
end
xlim(t([1 end]));
xlabel('Time from event onset (s)');
ylabel('Channel Number');
set(gca, 'YTickMode', 'manual');
set(gca, 'YTick', -1*nChannels:-1);
set(gca, 'YTickLabelMode', 'manual');
set(gca, 'YTickLabel', nChannels:-1:1);