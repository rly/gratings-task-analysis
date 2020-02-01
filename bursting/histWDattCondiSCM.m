function [percentageBurst, percBurstWMAttCT, cueTargetAttDiff, percBurstWMAttCA,...
    cueArrayAttDiff, cueTargetPoisson, cueArrayPoisson, percBurstWMPoissonCT,...
    percBurstWMPoissonCA, cueTargetShuffledDiff, cueArrayShuffledDiff] = histWDattCondiSCM(spikeTimes2use, startTime, endTime, UE, binWidth4hist,firstSpikeTimes, structName, plotOutput)
%structName = nUnits(uniti).name
unitStruct.name = structName(19:end-4);

% calculate percentage of spikes >= 200 Hz firing rate
percentageBurst =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);

targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
targetDimTrialAttOut = targetDimTrialAttOut - firstSpikeTimes;
targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut>0);

targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
targetDimTrialAttIn = targetDimTrialAttIn - firstSpikeTimes;
targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn>0);

arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - firstSpikeTimes;
arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut>0);

arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - firstSpikeTimes;
arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn>0);

cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
cueOnsetTrialAttOut = cueOnsetTrialAttOut - firstSpikeTimes;
cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut>0);

cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
cueOnsetTrialAttIn = cueOnsetTrialAttIn - firstSpikeTimes;
cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn>0);

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
for n=1:NEAttIn;
    %     nwinl=round(0.200*Fs);
    nwinr=round(targetDimTrialAttIn(n)*Fs);
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
    nwinr=round(targetDimTrialAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end

if ~isempty(datatmpAttIn) & ~isempty(datatmpAttOut)
    NEAttOut = size(attOut,2);
    NEAttIn = size(attIn,2);

    %     if plotOutput
    %         figure
    %         subplot(131)
    %         histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
    %         ylimsub1 = ylim;
    %         xlim([-10 500])
    %         title('real data Att In')
    %         subplot(132)
    %         plot(0,0)
    %         hold on
    %         histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
    %         ylimsub2 = ylim;
    %         xlim([-10 500])
    %         title('real data Att Out')
    %         subplot(133)
    %         histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
    %         hold on
    %         histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
    %         title('overlay both conditions')
    %         xlim([-10 500])
    %         subplot(131)
    %         ylim([0 max([ylimsub1 ylimsub2])])
    %         subplot(132)
    %         ylim([0 max([ylimsub1 ylimsub2])])
    %         subplot(133)
    %        ylim([0 max([ylimsub1 ylimsub2])])
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Repeat but with trial based SCM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear datatmpAttOut2use datatmpAttIn2use
    nSpikesAttIn = length(datatmpAttIn);
    nSpikesAttOut = length(datatmpAttOut);
    counter = 0;
    if nSpikesAttIn > nSpikesAttOut
        trialSel = randperm(NEAttIn);
        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
            counter = counter + 1;
            NEAttInNw = NEAttIn - counter;
            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); targetDimTrialAttInNw = targetDimTrialAttIn(trialSel(1:NEAttInNw));
            for n=1:NEAttInNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(targetDimTrialAttInNw(n)*Fs);
                indx=nEAttInNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                end
                clear indx
            end
            datatmpAttIn2use = datatmpAttIn;
            nSpikesAttIn = length(datatmpAttIn2use);
            datatmpAttIn = [];
        end
    elseif nSpikesAttOut > nSpikesAttIn
        trialSel = randperm(NEAttOut);
        while nSpikesAttOut + mean(nSpikesAttOutPerTrial) / 2 > nSpikesAttIn
            counter = counter + 1;
            NEAttOutNw = NEAttOut - counter;
            nEAttOutNw = nEAttOut(trialSel(1:NEAttOutNw)); targetDimTrialAttOutNw = targetDimTrialAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(targetDimTrialAttOutNw(n)*Fs);
                indx=nEAttOutNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                end
                clear indx
            end
            datatmpAttOut2use = datatmpAttOut;
            nSpikesAttOut = length(datatmpAttOut2use);
            datatmpAttOut = [];
        end
%     elseif nSpikesAttOut == nSpikesAttIn
%         AttInMore = 2;
    end

    if nSpikesAttIn > nSpikesAttOut
        AttInMore = 1;
    elseif nSpikesAttIn < nSpikesAttOut
        AttInMore = 0;
    end

    if exist('datatmpAttIn2use','var') == 0
        datatmpAttIn2use = datatmpAttIn;
    end
    if exist('datatmpAttOut2use','var') == 0
        datatmpAttOut2use = datatmpAttOut;
    end
    
    % create Poisson data for comparison
    if AttInMore == 1; nTrials = NEAttOut; else nTrials = NEAttIn; end
    datatmpPoisson = [];
    if ~AttInMore
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttIn,targetDimTrialAttIn,cueOnsetTrialAttIn,Fs,data);
        spikeTimesAllTrials = [];
        for triali = 1:length(trialSel)-counter
            nSpikesPerTrial(triali) = length(attOut(trialSel(triali)).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attOut(trialSel(triali)).spikeTimes];
        end
        nSpikesPerTrialAttOut = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:length(trialSel)-counter
            if triali == 1
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(1:nSpikesPerTrial(1));
                    idx = nSpikesPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(idx:idx+nSpikesPerTrial(triali)-1);
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
        for triali = 1:length(attIn)
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
        for triali = 1:length(attIn)
            if triali == 1
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(1:nSpikesPerTrial(1));
                    idx = nSpikesPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(idx:idx+nSpikesPerTrial(triali)-1);
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
    elseif AttInMore == 1
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttOut,targetDimTrialAttOut,cueOnsetTrialAttOut,Fs,data);
        spikeTimesAllTrials = [];
        for triali = 1:length(trialSel)-counter
            nSpikesPerTrial(triali) = length(attIn(trialSel(triali)).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attIn(trialSel(triali)).spikeTimes];
        end
        nSpikesPerTrialAttIn = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:length(trialSel)-counter
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
        
        clear spikeTimesAllTrials2use nSpikesPerTrial
        spikeTimesAllTrials = [];
        for triali = 1:length(attOut)
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
        for triali = 1:length(attOut)
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
%     else
%         [~,datatmpPoisson] = calcPoisson(NEAttOut,nEAttOut,targetDimTrialAttOut,cueOnsetTrialAttOut,Fs,data);
%         spikeTimesAllTrials = [];
%         for triali = 1:NEAttIn
%             nSpikesPerTrial(triali) = length(attIn(triali).spikeTimes);
%             spikeTimesAllTrials = [spikeTimesAllTrials attIn(triali).spikeTimes];
%         end
%         nSpikesPerTrialAttIn = nSpikesPerTrial;
%         spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
%         % make sure that there are no neighboring spikes with the same
%         % spiketimes
%         while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
%             spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
%         end
%         % reconstruct trials and extract interspike intervals
%         for triali = 1:NEAttIn
%             if triali == 1
%                 if nSpikesPerTrial(triali) ~= 0
%                     nspike2use = spikeTimesAllTrials2use(1:nSpikesPerTrial(1));
%                     idx = nSpikesPerTrial(1) + 1;
%                 else
%                     idx = 1;
%                     continue
%                 end
%             else
%                 if nSpikesPerTrial(triali) ~= 0
%                     nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesPerTrial(triali)-1);
%                     idx = nSpikesPerTrial(triali) + idx;
%                 else
%                     attInRandomized(triali).spikeTimes = [];
%                     continue
%                 end
%             end
%             attInRandomized(triali).spikeTimes = sort(nspike2use);
%             clear nspike2use
%         end
%         % finally get the interspike intervals of the shuffled data
%         datatmpShuffledAttIn = [];
%         for n = 1:length(attInRandomized)
%             datatmpShuffledAttIn = [datatmpShuffledAttIn diff(attInRandomized(n).spikeTimes)];
%         end
%         
%         clear spikeTimesAllTrials2use nSpikesPerTrial
%         spikeTimesAllTrials = [];
%         for triali = 1:NEAttOut
%             nSpikesPerTrial(triali) = length(attOut(triali).spikeTimes);
%             spikeTimesAllTrials = [spikeTimesAllTrials attOut(triali).spikeTimes];
%         end
%         nSpikesPerTrialAttOut = nSpikesPerTrial;
%         spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
%         % make sure that there are no neighboring spikes with the same
%         % spiketimes
%         while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
%             spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
%         end
%         % reconstruct trials and extract interspike intervals
%         for triali = 1:NEAttOut
%             if triali == 1
%                 if nSpikesPerTrial(triali) ~= 0
%                     nspike2use = spikeTimesAllTrials2use(1:nSpikesPerTrial(1));
%                     idx = nSpikesPerTrial(1) + 1;
%                 else
%                     idx = 1;
%                     continue
%                 end
%             else
%                 if nSpikesPerTrial(triali) ~= 0
%                     nspike2use = spikeTimesAllTrials2use(idx:idx+nSpikesPerTrial(triali)-1);
%                     idx = nSpikesPerTrial(triali) + idx;
%                 else
%                     attOutRandomized(triali).spikeTimes = [];
%                     continue
%                 end
%             end
%             attOutRandomized(triali).spikeTimes = sort(nspike2use);
%             clear nspike2use
%         end
%         % finally get the interspike intervals of the shuffled data
%         datatmpShuffledAttOut = [];
%         for n = 1:length(attOutRandomized)
%             datatmpShuffledAttOut = [datatmpShuffledAttOut diff(attOutRandomized(n).spikeTimes)];
%         end
    end
    
    totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
    percBurstWMAttCT = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
    percBurstWMPoissonCT = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);
    
    figure
    subplot(141)
    histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub1 = ylim;
    xlim([-10 500])
    title('real data Att In')
    subplot(142)
    plot(0,0)
    hold on
    histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub2 = ylim;
    xlim([-10 500])
    title('real data Att Out')
    subplot(143)
    attIn4statsCT = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    hold on
    attOut4statsCT = histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    poisson4statsCT = histogram(datatmpPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
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
    subplot(142)
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.95, ...
        {sprintf('Womelsdorf method \n (200Hz) burst %%: %0.2f', percentageBurst)});
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    subplot(143)
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.85, ...
        {sprintf('Cue target (200Hz) \n burst %%: %0.2f', percBurstWMAttCT),...
        sprintf('First bin %d ms \n AttIn %0.4f - \n AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsCT.Values(1), attOut4statsCT.Values(1), attIn4statsCT.Values(1) - attOut4statsCT.Values(1)),...
        sprintf('Poisson %0.4f', poisson4statsCT.Values(1))});
    subplot(144)
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.85, ...
        {sprintf('SHUFFLED DATA \n First bin %d ms \n AttIn %0.4f - \n AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsSHUFCT.Values(1), attOut4statsSHUFCT.Values(1), attIn4statsSHUFCT.Values(1) - attOut4statsSHUFCT.Values(1))});
    if plotOutput
        cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
        saveas(gcf,['HistogramISIrealCueTargetDimSCM' unitStruct.name '_' num2str(binWidth4hist) 'ms.png'])
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    end
    cueTargetAttDiff = attIn4statsCT.Values(1) - attOut4statsCT.Values(1);
    cueTargetShuffledDiff = attIn4statsSHUFCT.Values(1) - attOut4statsSHUFCT.Values(1);
    cueTargetPoisson = poisson4statsCT.Values(1);
    clear datatmpPoisson
else
    percentageBurst       = nan;
    percBurstWMAttCT      = nan;
    cueTargetAttDiff      = nan;
    percBurstWMAttCA      = nan;
    cueArrayAttDiff       = nan;
    cueTargetPoisson      = nan;
    cueArrayPoisson       = nan;
    percBurstWMPoissonCT  = nan;
    percBurstWMPoissonCA  = nan;
    cueTargetShuffledDiff = nan;
    cueArrayShuffledDiff  = nan;
end

%% repeat for array onset
clear datatmp*

% 200ms post-cue until ARRAY ONSET
data = binarySpikeTrain;
Fs = 1000;
NEAttIn=length(cueOnsetTrialAttIn);NEAttOut=length(cueOnsetTrialAttOut);
nEAttIn=floor(cueOnsetTrialAttIn*Fs)+1;nEAttOut=floor(cueOnsetTrialAttOut*Fs)+1;
datatmpAttIn=[];datatmpAttOut=[];
for n=1:NEAttIn;
    %     nwinl=round(0.200*Fs);
    nwinr=round(arrayOnsetTrialAttIn(n)*Fs);
    indx=nEAttIn(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
        nSpikesAttInPerTrialao(n) = sum(data(indx));
    end
end
for n=1:NEAttOut;
    %     nwinl=round(0.200*Fs);
    nwinr=round(arrayOnsetTrialAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
    end
end
if ~isempty(datatmpAttIn) & ~isempty(datatmpAttOut)
    %     if plotOutput
    %         figure
    %         subplot(131)
    %         histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
    %         ylimsub1 = ylim;
    %         xlim([-10 500])
    %         title('real data Att In')
    %         subplot(132)
    %         plot(0,0)
    %         hold on
    %         histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
    %         ylimsub2 = ylim;
    %         xlim([-10 500])
    %         title('real data Att Out')
    %         subplot(133)
    %         histogram(datatmpAttIn,'BinWidth',binWidth4hist,'Normalization','probability')
    %         hold on
    %         histogram(datatmpAttOut,'BinWidth',binWidth4hist,'Normalization','probability')
    %         title('overlay both conditions')
    %         xlim([-10 500])
    %         subplot(131)
    %         ylim([0 max([ylimsub1 ylimsub2])])
    %         subplot(132)
    %         ylim([0 max([ylimsub1 ylimsub2])])
    %         subplot(133)
    %         ylim([0 max([ylimsub1 ylimsub2])])
    %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Repeat but with trial based SCM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nSpikesAttIn = length(datatmpAttIn);
    nSpikesAttOut = length(datatmpAttOut);
    counter = 0;
    if nSpikesAttIn > nSpikesAttOut
        AttInMore = 1;
        trialSel = randperm(NEAttIn);
        while nSpikesAttIn + mean(nSpikesAttInPerTrial) / 2 > nSpikesAttOut
            counter = counter + 1;
            NEAttInNw = NEAttIn - counter;
            nEAttInNw = nEAttIn(trialSel(1:NEAttInNw)); arrayOnsetTrialAttInNw = arrayOnsetTrialAttIn(trialSel(1:NEAttInNw));
            for n=1:NEAttInNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(arrayOnsetTrialAttInNw(n)*Fs);
                indx=nEAttInNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttIn=[datatmpAttIn diff(find(data(indx)))];
                end
            end
            datatmpAttIn2use = datatmpAttIn;
            nSpikesAttIn = length(datatmpAttIn2use);
            datatmpAttIn = [];
        end
    elseif nSpikesAttOut > nSpikesAttIn
        AttInMore = 0;
        trialSel = randperm(NEAttOut);
        while nSpikesAttOut + mean(nSpikesAttOutPerTrial) / 2 > nSpikesAttIn
            counter = counter + 1;
            NEAttOutNw = NEAttOut - counter;
            nEAttOutNw = nEAttOut(trialSel(1:NEAttOutNw)); arrayOnsetTrialAttOutNw = arrayOnsetTrialAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(arrayOnsetTrialAttOutNw(n)*Fs);
                indx=nEAttOutNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                end
            end
            datatmpAttOut2use = datatmpAttOut;
            nSpikesAttOut = length(datatmpAttOut2use);
            datatmpAttOut = [];
        end
    elseif nSpikesAttOut == nSpikesAttIn
        if targetDimTrialAttIn(end) > targetDimTrialAttOut(end)
            AttInMore = 1;
        elseif targetDimTrialAttIn(end) < targetDimTrialAttOut(end)
            AttInMore = 0;
        end
    end
    if exist('datatmpAttIn2use','var') == 0
        datatmpAttIn2use = datatmpAttIn;
    end
    if exist('datatmpAttOut2use','var') == 0
        datatmpAttOut2use = datatmpAttOut;
    end
    
    % create Poisson data for comparison
    if AttInMore; nTrials = NEAttOut; else nTrials = NEAttIn; end
    datatmpPoisson = [];
    if ~AttInMore
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttIn,targetDimTrialAttIn,cueOnsetTrialAttIn,Fs,data);
        spikeTimesAllTrials = [];
        for triali = 1:length(trialSel)-counter
            nSpikesPerTrial(triali) = length(attOut(trialSel(triali)).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attOut(trialSel(triali)).spikeTimes];
        end
        nSpikesPerTrialAttOut = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:length(trialSel)-counter
            if triali == 1
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(1:nSpikesPerTrial(1));
                    idx = nSpikesPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(idx:idx+nSpikesPerTrial(triali)-1);
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
        for triali = 1:length(attIn)
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
        for triali = 1:length(attIn)
            if triali == 1
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(1:nSpikesPerTrial(1));
                    idx = nSpikesPerTrial(1) + 1;
                else
                    idx = 1;
                    continue
                end
            else
                if nSpikesPerTrial(triali) ~= 0
                    nspike2use = spikeTimesAllTrials(idx:idx+nSpikesPerTrial(triali)-1);
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
    else
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttOut,targetDimTrialAttOut,cueOnsetTrialAttOut,Fs,data);
        spikeTimesAllTrials = [];
        for triali = 1:length(trialSel)-counter
            nSpikesPerTrial(triali) = length(attIn(trialSel(triali)).spikeTimes);
            spikeTimesAllTrials = [spikeTimesAllTrials attIn(trialSel(triali)).spikeTimes];
        end
        nSpikesPerTrialAttIn = nSpikesPerTrial;
        spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        % make sure that there are no neighboring spikes with the same
        % spiketimes
        while sum(diff(spikeTimesAllTrials2use)==0) ~= 0
            spikeTimesAllTrials2use = spikeTimesAllTrials(randperm(length(spikeTimesAllTrials)));
        end
        % reconstruct trials and extract interspike intervals
        for triali = 1:length(trialSel)-counter
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
        
        clear spikeTimesAllTrials2use nSpikesPerTrial
        spikeTimesAllTrials = [];
        for triali = 1:length(attOut)
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
        for triali = 1:length(attOut)
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
    end
    
    totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
    percBurstWMAttCA = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
    percBurstWMPoissonCA = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);
    
    figure
    subplot(141)
    histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub1 = ylim;
    xlim([-10 500])
    title('real data Att In')
    subplot(142)
    plot(0,0)
    hold on
    histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub2 = ylim;
    xlim([-10 500])
    title('real data Att Out')
    subplot(143)
    attIn4statsCA = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    hold on
    attOut4statsCA = histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    poisson4statsCA = histogram(datatmpPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    xlim([-10 500])
    title('overlay both condi + Poisson')
    ylimsub3 = ylim;
    subplot(144)
    attIn4statsSHUFCA = histogram(datatmpShuffledAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    hold on
    attOut4statsSHUFCA = histogram(datatmpShuffledAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    title('overlay shuffled data')
    ylimsub4 = ylim;
    xlim([-10 500])
    subplot(141)
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    subplot(142)
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.95, ...
        {sprintf('Womelsdorf method \n (200Hz) burst %%: %0.2f', percentageBurst)});
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    subplot(143)
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.85, ...
        {sprintf('Cue array (200Hz) \n burst %%: %0.2f', percBurstWMAttCT),...
        sprintf('First bin %d ms \n AttIn %0.4f - \n AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsCA.Values(1), attOut4statsCA.Values(1), attIn4statsCA.Values(1) - attOut4statsCA.Values(1)),...
        sprintf('Poisson %0.4f', poisson4statsCT.Values(1))});
    subplot(144)
    ylim([0 max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])])
    text(40,max([ylimsub1 ylimsub2 ylimsub3 ylimsub4])*0.85, ...
        {sprintf('SHUFFLED DATA \n First bin %d ms \n AttIn %0.4f - \n AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsSHUFCA.Values(1), attOut4statsSHUFCA.Values(1), attIn4statsSHUFCA.Values(1) - attOut4statsSHUFCA.Values(1))});
    if plotOutput
        cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
        saveas(gcf,['HistogramISIrealCueArrayOnsetSCM' unitStruct.name '_' num2str(binWidth4hist) 'ms.png'])
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    end
    
    cueArrayAttDiff = attIn4statsCA.Values(1) - attOut4statsCA.Values(1);
    cueArrayShuffledDiff = attIn4statsSHUFCA.Values(1) - attOut4statsSHUFCA.Values(1);
    cueArrayPoisson = poisson4statsCA.Values(1);
    clear datatmpPoisson
else
    percentageBurst       = nan;
    percBurstWMAttCT      = nan;
    cueTargetAttDiff      = nan;
    percBurstWMAttCA      = nan;
    cueArrayAttDiff       = nan;
    cueTargetPoisson      = nan;
    cueArrayPoisson       = nan;
    percBurstWMPoissonCT  = nan;
    percBurstWMPoissonCA  = nan;
    cueTargetShuffledDiff = nan;
    cueArrayShuffledDiff  = nan;
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