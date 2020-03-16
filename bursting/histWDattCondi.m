function [percentageBurst, percBurstWMAttCT, percBurstWMPoissonCT, realAttIn, realAttOut,...
    shuffledAttIn, shuffledAttOut,poissonAttIn, poissonAttOut] = histWDattCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut, binWidth4hist, plotOutput, saveFileName)

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;

% calculate percentage of spikes >= 200 Hz firing rate
percentageBurst =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);

data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

% calculate trial based spiking
% 200ms post-cue until TARGET DIM
data = binarySpikeTrain;
Fs = 1000;
NEAttIn=length(cueOnsetTrialAttIn); NEAttOut=length(cueOnsetTrialAttOut);
nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1; nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
datatmpAttIn=[];datatmpAttOut=[];
for n=1:NEAttIn
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttIn(n)*Fs);
    indx=nEAttIn(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
        nSpikesAttInPerTrial(n) = sum(data(indx));
        attIn(n).spikeTimes=find(data(indx));
        clear indx
    end
end
for n=1:NEAttOut;
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end
NEAttIn = length(attIn); NEAttOut = length(attOut);

% create Poisson data for comparison
[~,datatmpPoissonAttIn] = calcPoisson(NEAttIn,nEAttIn,endTimeAttIn,cueOnsetTrialAttIn,Fs,data);
[~,datatmpPoissonAttOut] = calcPoisson(NEAttOut,nEAttOut,endTimeAttOut,cueOnsetTrialAttOut,Fs,data);
spikeTimesAllTrials = [];
for triali = 1:NEAttOut
    nSpikesPerTrial(triali) = length(attOut(triali).spikeTimes);
    spikeTimesAllTrials = [spikeTimesAllTrials attOut(triali).spikeTimes];
end
nSpikesPerTrialAttOut = nSpikesPerTrial;
spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
% make sure that there are no neighboring spikes with the same
% spiketimes
while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
    spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
end
% reconstruct trials and extract interspike intervals
for triali = 1:NEAttOut
    if triali == 1
        if nSpikesPerTrial(triali) ~= 0
            nspike2use = spikeTimesAllTrials2use(1:nSpikesPerTrial(1));
            idx = nSpikesPerTrial(1) + 1;
        else
            idx = 1;
            continue
        end
    else
        if nSpikesPerTrial(triali) ~= 0
            nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesPerTrial(triali)-1);
            idx = nSpikesPerTrial(triali) + idx;
        else
            attOutRandomized(triali).spikeTimes = [];
            continue
        end
    end
    attOutRandomized(triali).spikeTimes = sort(nspike2use);
    clear nspike2use
end
% finally get the interspike intervals of the shuffled data
datatmpShuffledAttOut = [];
for n = 1:length(attOutRandomized)
    datatmpShuffledAttOut = [datatmpShuffledAttOut diff(attOutRandomized(n).spikeTimes)];
end

clear spikeTimesAllTrials2use nSpikesPerTrial
spikeTimesAllTrials = [];
for triali = 1:NEAttIn
    nSpikesPerTrial(triali) = length(attIn(triali).spikeTimes);
    spikeTimesAllTrials = [spikeTimesAllTrials attIn(triali).spikeTimes];
end
nSpikesPerTrialAttIn = nSpikesPerTrial;
spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
% make sure that there are no neighboring spikes with the same
% spiketimes
while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
    spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
end
% reconstruct trials and extract interspike intervals
for triali = 1:NEAttIn
    if triali == 1
        if nSpikesPerTrial(triali) ~= 0
            nspike2use = spikeTimesAllTrials2use(1:nSpikesPerTrial(1));
            idx = nSpikesPerTrial(1) + 1;
        else
            idx = 1;
            continue
        end
    else
        if nSpikesPerTrial(triali) ~= 0
            nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesPerTrial(triali)-1);
            idx = nSpikesPerTrial(triali) + idx;
        else
            attInRandomized(triali).spikeTimes = [];
            continue
        end
    end
    attInRandomized(triali).spikeTimes = sort(nspike2use);
    clear nspike2use
end
% finally get the interspike intervals of the shuffled data
datatmpShuffledAttIn = [];
for n = 1:length(attInRandomized)
    datatmpShuffledAttIn = [datatmpShuffledAttIn diff(attInRandomized(n).spikeTimes)];
end

totalAttConds = [datatmpAttIn datatmpAttOut];
percBurstWMAttCT = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
totalAttShuffledConds = [datatmpShuffledAttIn datatmpShuffledAttOut];
percBurstWMAttShuffledCT = sum(totalAttShuffledConds <= 5) * 100 / size(totalAttShuffledConds,2);
datatmpPoisson = [datatmpPoissonAttIn datatmpPoissonAttOut];
percBurstWMPoissonCT = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);

figure
subplot(141)
histogram(datatmpAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
hold on
plot(0,0)
attIn4statsCTPoisson = histogram(datatmpPoissonAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
ylimsub1 = ylim;
xlim([-10 500])
title('real data Att In')
subplot(142)
plot(0,0)
hold on
histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
attOut4statsCTPoisson = histogram(datatmpPoissonAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
ylimsub2 = ylim;
xlim([-10 500])
title('real data Att Out')
subplot(143)
attIn4statsCT = histogram(datatmpAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
hold on
attOut4statsCT = histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
xlim([-10 500])
title('overlay both condi + Poisson')
ylimsub3 = ylim;
subplot(144)
attIn4statsSHUFCT = histogram(datatmpShuffledAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
hold on
attOut4statsSHUFCT = histogram(datatmpShuffledAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
title('overlay shuffled data')
ylimsub4 = ylim;
xlim([-10 500])
subplot(141)
ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.9, ...
    {sprintf('First bin %d ms \n AttIn %0.4f', binWidth4hist, attIn4statsCT.Values(1)),...
    sprintf('AttIn Poisson %0.4f', attIn4statsCTPoisson.Values(1))}, 'FontSize',6);
subplot(142)
ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.9, ...
    {sprintf('First bin %d ms \n AttOut %0.4f', binWidth4hist, attOut4statsCT.Values(1)),...
    sprintf('AttOut Poisson %0.4f', attOut4statsCTPoisson.Values(1))}, 'FontSize',6);
subplot(143)
ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.9, ...
    {sprintf('Womelsdorf method \n (200Hz) burst %%: %0.2f', percentageBurst),...
    sprintf('Cue target (200Hz) \n burst %%: %0.2f', percBurstWMAttCT),...
    sprintf('First bin %d ms \n Real difference = %0.4f', ...
    binWidth4hist, attIn4statsCT.Values(1) - attOut4statsCT.Values(1)),...
    sprintf('Poisson difference = %0.4f', attIn4statsCTPoisson.Values(1) - attOut4statsCTPoisson.Values(1))}, 'FontSize',6);
subplot(144)
ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.9, ...
    {sprintf('Cue target \n (200Hz) burst %%: %0.2f', percBurstWMAttShuffledCT),...
    sprintf('First bin %d ms \n AttIn %0.4f - \n AttOut %0.4f \n = %0.4f', ...
    binWidth4hist, attIn4statsSHUFCT.Values(1), attOut4statsSHUFCT.Values(1), attIn4statsSHUFCT.Values(1) - attOut4statsSHUFCT.Values(1))}, 'FontSize',6);
if plotOutput
    cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
    saveas(gcf,saveFileName)
    cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
end

realAttIn = attIn4statsCT.Values(1);
realAttOut = attOut4statsCT.Values(1);
shuffledAttIn = attIn4statsSHUFCT.Values(1);
shuffledAttOut = attOut4statsSHUFCT.Values(1);
poissonAttIn = attIn4statsCTPoisson.Values(1);
poissonAttOut = attOut4statsCTPoisson.Values(1);

cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
save(saveFileName(1:end-4),'attIn4statsCT','attOut4statsCT','attIn4statsSHUFCT',...
    'attOut4statsSHUFCT','attIn4statsCTPoisson','attOut4statsCTPoisson')
cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')

end



% attIn4statsAll(idx).Values = attIn4stats.Values;
% attOut4statsAll(idx).Values = attOut4stats.Values;
% attIn4statsCTAll(idx).Values = attIn4statsCT.Values;
% attOut4statsCTAll(idx).Values = attOut4statsCT.Values;
% wdPercentageAll(idx) = percentageBurst;
%
% idx = idx + 1;
% close all
%
% clear binarySpikeTrain
%
% save(['AttInMinusAttOut_' num2str(binWidth4hist) 'ms'],'cueArrayAttDiff','cueTargetAttDiff',...
%     'attIn4statsAll', 'attOut4statsAll', 'attIn4statsCTAll', 'attOut4statsCTAll',...
%     'wdPercentageAll')
% [H,P,CI,STATS] = ttest(cueArrayAttDiff)
% [H,P,CI,STATS] = ttest(cueTargetAttDiff)
%
% % [~,indices] = sort(attIn4statsAll(end).Values)
% % [~,indicesCT] = sort(attIn4statsCTAll(end).Values)
% % figure
% % subplot(121)
% % plot(attIn4statsAll(end).Values(indices),'.')%(attIn4statsAll(end).Values>0))
% % hold on
% % plot(attOut4statsAll(end).Values(indices),'.')%(attOut4statsAll(end).Values>0))
% % subplot(122)
% % plot(attIn4statsCTAll(end).Values(indicesCT),'.')%(attIn4statsAll(end).Values>0))
% % hold on
% % plot(attOut4statsCTAll(end).Values(indicesCT),'.')%(attOut4statsAll(end).Values>0))
%
% figure
% histogram(cueArrayAttDiff)
% hold on
% histogram(cueTargetAttDiff)






% eof