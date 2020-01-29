function [percentageBurst, percBurstWMAttCT, cueTargetAttDiff, percBurstWMAttCA,...
 cueArrayAttDiff, cueTargetPoisson, cueArrayPoisson, percBurstWMPoissonCT,...
 percBurstWMPoissonCA] = histWDattCondi(spikeTimes2use, startTime, endTime, UE, binWidth4hist,firstSpikeTimes, structName, plotOutput)
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
        clear indx
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
        AttInMore = 1;
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
            nEAttOutNw = nEAttOut(trialSel(1:NEAttOutNw)); targetDimTrialAttOutNw = targetDimTrialAttOut(trialSel(1:NEAttOutNw));
            for n=1:NEAttOutNw;
                %     nwinl=round(0.200*Fs);
                nwinr=round(targetDimTrialAttOutNw(n)*Fs);
                indx=nEAttOutNw(n):nwinr-1;
                if length(indx) >1 && indx(end) < size(data,2)
                    datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
                end
            end
            datatmpAttOut2use = datatmpAttOut;
            nSpikesAttOut = length(datatmpAttOut2use);
            datatmpAttOut = [];
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
    else
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttOut,targetDimTrialAttOut,cueOnsetTrialAttOut,Fs,data);
    end
     
    totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
    percBurstWMAttCT = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
    percBurstWMPoissonCT = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);

    figure
    subplot(131)
    histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub1 = ylim;
    xlim([-10 500])
    title('real data Att In')
    subplot(132)
    plot(0,0)
    hold on
    histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub2 = ylim;
    xlim([-10 500])
    title('real data Att Out')
    subplot(133)
    attIn4statsCT = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    hold on
    attOut4statsCT = histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    poisson4statsCT = histogram(datatmpPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    title('overlay both conditions')
    xlim([-10 500])
    subplot(131)
    ylim([0 max([ylimsub1 ylimsub2])])
    subplot(132)
    ylim([0 max([ylimsub1 ylimsub2])])
    subplot(133)
    ylim([0 max([ylimsub1 ylimsub2])])
    text(20,max([ylimsub1 ylimsub2])*0.9, ...
        {sprintf('Womelsdorf method (200Hz) \n burst %%: %0.2f', percentageBurst),...
        sprintf('Cue target (200Hz) \n burst %%: %0.2f', percBurstWMAttCT),...
        sprintf('First bin %d ms \n AttIn %0.4f - AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsCT.Values(1), attOut4statsCT.Values(1), attIn4statsCT.Values(1) - attOut4statsCT.Values(1)),...
        sprintf('Poisson %0.4f', poisson4statsCT.Values(1))});
    if plotOutput
        cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
        saveas(gcf,['HistogramISIrealCueTargetDimSCM' unitStruct.name '_' num2str(binWidth4hist) 'ms.png'])
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    end
    cueTargetAttDiff = attIn4statsCT.Values(1) - attOut4statsCT.Values(1);
    cueTargetPoisson = poisson4statsCT.Values(1);
    clear datatmpPoisson
else
    percentageBurst  = nan;
    percBurstWMAttCT = nan;
    cueTargetAttDiff = nan;
    percBurstWMAttCA = nan;
    cueArrayAttDiff  = nan;
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
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttIn,arrayOnsetTrialAttIn,cueOnsetTrialAttIn,Fs,data);
    else
        [~,datatmpPoisson] = calcPoisson(nTrials,nEAttOut,arrayOnsetTrialAttOut,cueOnsetTrialAttOut,Fs,data);
    end

    totalAttConds = [datatmpAttIn2use datatmpAttOut2use];
    percBurstWMAttCA = sum(totalAttConds <= 5) * 100 / size(totalAttConds,2);
    percBurstWMPoissonCA = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);

    figure
    subplot(131)
    histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub1 = ylim;
    xlim([-10 500])
    title('real data Att In')
    subplot(132)
    plot(0,0)
    hold on
    histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability')
    ylimsub2 = ylim;
    xlim([-10 500])
    title('real data Att Out')
    subplot(133)
    attIn4statsCA = histogram(datatmpAttIn2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    hold on
    attOut4statsCA = histogram(datatmpAttOut2use,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    poisson4statsCA = histogram(datatmpPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    title('overlay both conditions')
    xlim([-10 500])
    subplot(131)
    ylim([0 max([ylimsub1 ylimsub2])])
    subplot(132)
    ylim([0 max([ylimsub1 ylimsub2])])
    subplot(133)
    ylim([0 max([ylimsub1 ylimsub2])])
    text(20,max([ylimsub1 ylimsub2])*0.9, ...
        {sprintf('Womelsdorf method (200Hz) \n burst %%: %0.2f', percentageBurst),...
        sprintf('Cue array (200Hz) \n burst %%: %0.2f', percBurstWMAttCA),...
        sprintf('First bin %d ms \n AttIn %0.4f - AttOut %0.4f \n = %0.4f', ...
        binWidth4hist, attIn4statsCA.Values(1), attOut4statsCA.Values(1), attIn4statsCA.Values(1) - attOut4statsCA.Values(1)),...
        sprintf('Poisson %0.4f', poisson4statsCA.Values(1))});
    if plotOutput
        cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
        saveas(gcf,['HistogramISIrealCueArrayOnsetSCM' unitStruct.name '_' num2str(binWidth4hist) 'ms.png'])
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    end

    cueArrayAttDiff = attIn4statsCA.Values(1) - attOut4statsCA.Values(1);
    cueArrayPoisson = poisson4statsCA.Values(1);
    clear datatmpPoisson
else
    percentageBurst  = nan;
    percBurstWMAttCT = nan;
    cueTargetAttDiff = nan;
    percBurstWMAttCA = nan;
    cueArrayAttDiff  = nan;
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