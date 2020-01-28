clear all 
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally')
nUnits = dir('spikeTimes2use*');

binWidth4hist = 5;
plotOutput = 0;

percentageBurst  = nan(length(nUnits),1);
percBurstWMAttCT = nan(length(nUnits),1);
cueTargetAttDiff = nan(length(nUnits),1);
percBurstWMAttCA = nan(length(nUnits),1);
cueArrayAttDiff  = nan(length(nUnits),1);
for uniti = 1:length(nUnits)
    cd('/Users/labmanager/Documents/MATLAB/Data locally')
    load(nUnits(uniti).name)
    cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    [percentageBurst(uniti), percBurstWMAttCT(uniti), cueTargetAttDiff(uniti), percBurstWMAttCA(uniti), cueArrayAttDiff(uniti)] = histWDattCondi(spikeTimes2use, startTime, endTime, UE, binWidth4hist, firstSpikeTimes, nUnits(uniti).name, plotOutput);
    clear UE endTime firstSpikeTimes spikeTimes2use startTime
    close all
end

cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
data = [percentageBurst, percBurstWMAttCT, cueTargetAttDiff, percBurstWMAttCA, cueArrayAttDiff];
dataTbl = array2table(data);dataTbl.Properties.VariableNames = {'wdAll' 'wdCueTarg' 'cueTargDiff' 'wdCueArray' 'cueArrayDiff'};

save(['AttInMinusAttOutNw_' num2str(binWidth4hist) 'ms'],'data','dataTbl')

[~,sortedBC] = sort(data(:,1))
figure
plot(data(sortedBC,1))
hold on
plot(data(sortedBC,2))
plot(data(sortedBC,4))

nanmean(data(:,1))
sum(data(:,1)>1)

nanmean(data(:,2))
sum(data(:,2)>1)

nanmean(data(:,4))
sum(data(:,4)>1)

nanmean(data(:,3))
nanmean(data(:,5))

[H,P,CI,STATS] = ttest(data(data(:,1)>1,3))
[H,P,CI,STATS] = ttest(data(data(:,2)>1,3))

[H,P,CI,STATS] = ttest(data(data(:,1)>1,5))
[H,P,CI,STATS] = ttest(data(data(:,4)>1,5))





% eof