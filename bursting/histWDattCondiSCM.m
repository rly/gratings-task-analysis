function [percentageBurst, percBurstWMAttCT, percBurstWMPoissonCT, realAttIn, realAttOut,...
    shuffledAttIn, shuffledAttOut,poissonAttIn, poissonAttOut] = histWDattCondiSCM(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut, binWidth4hist, plotOutput, saveFileName)

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;
%endTimeAttIn = targetDimTrialAttIn; endTimeAttOut =targetDimTrialAttOut;
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

% Perform spike count matching
clear datatmpAttOut2use datatmpAttIn2use
nSpikesAttIn = length(datatmpAttIn);
nSpikesAttOut = length(datatmpAttOut);
counter = 0;
if nSpikesAttIn > 10 && nSpikesAttOut > 10
    if nSpikesAttIn > nSpikesAttOut
        datatmpAttInOld = datatmpAttIn;
        trialSel = randperm(NEAttIn);
        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
            counter = counter + 1;
            NEAttInNw = NEAttIn - counter;
            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); endTimeAttInNw = endTimeAttIn(trialSel(1:NEAttInNw));
            for n=1:NEAttInNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttInNw(n)*Fs);
                indx=nEAttInNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                    attInNw(n).spikeTimes=find(data(indx));
                end
                clear indx
            end
            datatmpAttIn2use = datatmpAttIn;
            nSpikesAttIn = length(datatmpAttIn2use);
            datatmpAttIn = [];
        end
    elseif nSpikesAttOut > nSpikesAttIn
        datatmpAttOutOld = datatmpAttOut; datatmpAttOut = [];
        trialSel = randperm(NEAttOut);
        while nSpikesAttOut + mean(nSpikesAttOutPerTrial) / 2 > nSpikesAttIn
            counter = counter + 1;
            NEAttOutNw = NEAttOut - counter;
            nEAttOutNw = nEAttOut(trialSel(1:NEAttOutNw)); endTimeAttOutNw = endTimeAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(endTimeAttOutNw(n)*Fs);
                indx=nEAttOutNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                    attOutNw(n).spikeTimes=find(data(indx));
                end
                clear indx
            end
            datatmpAttOut2use = datatmpAttOut;
            nSpikesAttOut = length(datatmpAttOut2use);
            datatmpAttOut = [];
        end
    end
    if nSpikesAttOut > 10 && nSpikesAttIn > 10
        if exist('datatmpAttIn2use','var') == 0 && exist('datatmpAttOut2use','var') == 1
            datatmpAttIn2use = datatmpAttIn;
            endTimeAttInNw = endTimeAttIn;
            nEAttInNw = nEAttIn; NEAttInNw = NEAttIn;
            if exist('trialSel','var') == 0
                cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(1:NEAttOutNw);
            else
                cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(trialSel(1:NEAttOutNw));
            end
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn;
            attInNw = attIn;
        elseif exist('datatmpAttIn2use','var') == 1 && exist('datatmpAttOut2use','var') == 0
            datatmpAttOut2use = datatmpAttOut;
            endTimeAttOutNw = endTimeAttOut;
            nEAttOutNw = nEAttOut; NEAttOutNw = NEAttOut;
            if exist('trialSel','var') == 0
                cueOnsetTrialAttInNw = cueOnsetTrialAttIn(1:NEAttInNw);
            else
                cueOnsetTrialAttInNw = cueOnsetTrialAttIn(trialSel(1:NEAttInNw));
            end
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut;
            attOutNw = attOut;
        else
            datatmpAttIn2use = datatmpAttIn;
            endTimeAttInNw = endTimeAttIn;
            nEAttInNw = nEAttIn; NEAttInNw = NEAttIn;
            datatmpAttOut2use = datatmpAttOut;
            endTimeAttOutNw = endTimeAttOut;
            nEAttOutNw = nEAttOut; NEAttOutNw = NEAttOut;
            if exist('trialSel','var') == 0
                cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(1:NEAttOutNw);
            else
                cueOnsetTrialAttOutNw = cueOnsetTrialAttOut(trialSel(1:NEAttOutNw));
            end
            cueOnsetTrialAttInNw = cueOnsetTrialAttIn;
            attInNw = attIn;
            if exist('trialSel','var') == 0
                cueOnsetTrialAttInNw = cueOnsetTrialAttIn(1:NEAttInNw);
            else
                cueOnsetTrialAttInNw = cueOnsetTrialAttIn(trialSel(1:NEAttInNw));
            end
            cueOnsetTrialAttOutNw = cueOnsetTrialAttOut;
            attOutNw = attOut;
        end

        % create Poisson data for comparison
        [~,datatmpPoissonAttIn] = calcPoisson(NEAttInNw,endTimeAttInNw,cueOnsetTrialAttInNw,Fs,data);
        [~,datatmpPoissonAttOut] = calcPoisson(NEAttOutNw,endTimeAttOutNw,cueOnsetTrialAttOutNw,Fs,data);

        spikeTimesAllTrials = [];
        for triali = 1:NEAttOutNw
            nSpikesPerTrial(triali) = length(attOutNw(triali).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attOutNw(triali).spikeTimes];
        end
        nSpikesPerTrialAttOut = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:NEAttOutNw
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
        for triali = 1:NEAttInNw
            nSpikesPerTrial(triali) = length(attInNw(triali).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attInNw(triali).spikeTimes];
        end
        nSpikesPerTrialAttIn = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:NEAttInNw
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

        totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
        percBurstWMAttCT = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
        totalAttShuffledConds = [datatmpShuffledAttIn datatmpShuffledAttOut];
        percBurstWMAttShuffledCT = sum(totalAttShuffledConds <= 5) * 100 / size(totalAttShuffledConds,2);
        datatmpPoisson = [datatmpPoissonAttIn datatmpPoissonAttOut];
        percBurstWMPoissonCT = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);

        figure
        subplot(141)
        histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
        hold on
        plot(0,0)
        attIn4statsCTPoisson = histogram(datatmpPoissonAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
        ylimsub1 = ylim;
        xlim([-10 500])
        title('real data Att In')
        subplot(142)
        plot(0,0)
        hold on
        histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
        attOut4statsCTPoisson = histogram(datatmpPoissonAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
        ylimsub2 = ylim;
        xlim([-10 500])
        title('real data Att Out')
        subplot(143)
        attIn4statsCT = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
        hold on
        attOut4statsCT = histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
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
            %cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
            cd('/Volumes/HDD1/BurstSep4allHDD/')
            saveas(gcf,saveFileName)
            cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
        end

        realAttIn = attIn4statsCT.Values(1);
        realAttOut = attOut4statsCT.Values(1);
        shuffledAttIn = attIn4statsSHUFCT.Values(1);
        shuffledAttOut = attOut4statsSHUFCT.Values(1);
        poissonAttIn = attIn4statsCTPoisson.Values(1);
        poissonAttOut = attOut4statsCTPoisson.Values(1);
    else
        attIn4statsCT = NaN;
        attOut4statsCT = NaN;
        attIn4statsSHUFCT = NaN;
        attOut4statsSHUFCT = NaN;
        attIn4statsCTPoisson = NaN;
        attOut4statsCTPoisson = NaN;
        percBurstWMAttCT = NaN;
        percBurstWMPoissonCT = NaN;
        realAttIn = NaN;
        realAttOut = NaN;
        shuffledAttIn = NaN;
        shuffledAttOut = NaN;
        poissonAttIn = NaN;
        poissonAttOut = NaN;
    end
else
    attIn4statsCT = NaN;
    attOut4statsCT = NaN;
    attIn4statsSHUFCT = NaN;
    attOut4statsSHUFCT = NaN;
    attIn4statsCTPoisson = NaN;
    attOut4statsCTPoisson = NaN;
    percBurstWMAttCT = NaN;
    percBurstWMPoissonCT = NaN;
    realAttIn = NaN;
    realAttOut = NaN;
    shuffledAttIn = NaN;
    shuffledAttOut = NaN;
    poissonAttIn = NaN;
    poissonAttOut = NaN;
end

%cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
cd('/Volumes/HDD1/BurstSep4allHDD/')
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