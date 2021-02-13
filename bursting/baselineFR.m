clc
clear all
close all
cd('E:\grating-data-diffbranch')
nUnitsDir = dir('*.mat');

nUnits = length(nUnitsDir);
cd('E:\grating-analysis-output')

for uniti = 1:nUnits
    tic
    cd('E:\grating-data-diffbranch')
    load(nUnitsDir(uniti).name)
    if ismember(unitStruct.channelID,R.dPulChannels)
        pulLocV(uniti) = 0;
    elseif ismember(unitStruct.channelID,R.vPulChannels)
        pulLocV(uniti) = 1;
    end
    
    chanNum(uniti) = unitStruct.channelID;
    startTime = unitStruct.ts(1); endTime = unitStruct.ts(end);
    cueOnset = UE.cueOnset;
    baselineOnset = UE.cueOnset - 0.300;
    tsInRange = baselineOnset >= startTime & cueOnset <= endTime;
    cueOnset = cueOnset(tsInRange); baselineOnset = baselineOnset(tsInRange);    
    
    spikeTimes2use = unitStruct.ts; %spikeStruct.tsTask;
    firingRate(uniti) = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

    % compute FR count
    for cuei = 1:numel(baselineOnset)
        nwinr = cueOnset(cuei);
        nwinl = baselineOnset(cuei);
        baseline.spikeCount(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        baseline.frCount(cuei) = baseline.spikeCount(cuei) / round((nwinr - nwinl),3);
    end
     
    allBaseline(uniti).meanSpikeCount = mean(baseline.spikeCount);
    allBaseline(uniti).medianSpikeCount = median(baseline.spikeCount);
    allBaseline(uniti).meanFRCount = mean(baseline.frCount);
    allBaseline(uniti).medianFRCount = median(baseline.frCount);
    
    cueOnsetTrialInRF = UE.cueOnset(UE.cueLoc == 3);%UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3);
    cueOnsetTrialOutRF = UE.cueOnset(UE.cueLoc == 4);%UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 4);
    
    % compute visually responsiveness
    for cuei = 1:numel(cueOnsetTrialInRF)
        nwinr = cueOnsetTrialInRF(cuei);
        nwinl = cueOnsetTrialInRF(cuei) - 0.175;
        spikeCount.cueInRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.cueInRFBaseline(cuei) = spikeCount.cueInRFBaseline(cuei) / round((nwinr - nwinl),3);
    end
    for cuei = 1:numel(cueOnsetTrialInRF)
        nwinr = cueOnsetTrialInRF(cuei) + 0.200;
        nwinl = cueOnsetTrialInRF(cuei) + 0.025;
        spikeCount.cueInRF(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.cueInRF(cuei) = spikeCount.cueInRF(cuei) / round((nwinr - nwinl),3);
    end
    for cuei = 1:numel(cueOnsetTrialOutRF)
        nwinr = cueOnsetTrialOutRF(cuei);
        nwinl = cueOnsetTrialOutRF(cuei) - 0.175;
        spikeCount.cueOutRFBaseline(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.cueOutRFBaseline(cuei) = spikeCount.cueOutRFBaseline(cuei) / round((nwinr - nwinl),3);
    end
    for cuei = 1:numel(cueOnsetTrialOutRF)
        nwinr = cueOnsetTrialOutRF(cuei) + 0.200;
        nwinl = cueOnsetTrialOutRF(cuei) + 0.025;
        spikeCount.cueOutRF(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.cueOutRF(cuei) = spikeCount.cueOutRF(cuei) / round((nwinr - nwinl),3);
    end
    meanSCInRFCount(uniti) = mean(spikeCount.cueInRF);
    meanSCInRFBlCount(uniti) = mean(spikeCount.cueInRFBaseline);
    meanSCOutRFCount(uniti) = mean(spikeCount.cueOutRF);
    meanSCOutRFBlCount(uniti) = mean(spikeCount.cueOutRFBaseline);
    [pSCCIsum(uniti)] = ranksum(frCount.cueInRFBaseline,frCount.cueInRF);
    [pSCCOsum(uniti)] = ranksum(frCount.cueOutRFBaseline,frCount.cueOutRF);

    targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
    targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
    arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
    arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
    cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
    cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;

    tsInRangeAttOut = cueOnsetTrialAttOut >= startTime & targetDimTrialAttOut <= endTime;
    targetDimTrialAttOut = targetDimTrialAttOut(tsInRangeAttOut); 
    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(tsInRangeAttOut); 
    cueOnsetTrialAttOut = cueOnsetTrialAttOut(tsInRangeAttOut); 
    tsInRangeAttIn = cueOnsetTrialAttIn >= startTime & targetDimTrialAttIn <= endTime;
    targetDimTrialAttIn = targetDimTrialAttIn(tsInRangeAttIn); 
    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(tsInRangeAttIn); 
    cueOnsetTrialAttIn = cueOnsetTrialAttIn(tsInRangeAttIn); 
    
    % compute attentional modulation
    for cuei = 1:numel(cueOnsetTrialAttIn)
        nwinr = arrayOnsetTrialAttIn(cuei);
        nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
        spikeCount.attIn(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.attIn(cuei) = spikeCount.attIn(cuei) / round((nwinr - nwinl),3);
    end
    for cuei = 1:numel(cueOnsetTrialAttOut)
        nwinr = arrayOnsetTrialAttOut(cuei);
        nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
        spikeCount.attOut(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.attOut(cuei) = spikeCount.attOut(cuei) / round((nwinr - nwinl),3);
    end   
    meanSCAttInCount(uniti) = mean(spikeCount.attIn);
    meanSCAttOutCount(uniti) = mean(spikeCount.attOut);
    [pSCATTsum(uniti)] = ranksum(frCount.attIn,frCount.attOut);
    
    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
    toc
end
cd('E:\grating-analysis-output')
visResp = (pSCCIsum <= .05 & meanSCInRFCount > meanSCInRFBlCount) | (pSCCOsum <= .05 & meanSCOutRFCount > meanSCOutRFBlCount);
attMod = pSCATTsum <= .05;
pulLocV = boolean(pulLocV);
save('infoUnits','visResp','attMod','p*','mean*')

%%

chanNumAdj = chanNum; chanNumAdj(chanNumAdj>32) = chanNumAdj(chanNumAdj>32)-32;
figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV).meanFRCount],chanNumAdj(~pulLocV)+32,'k*')
hold on
plot([allBaseline(pulLocV).meanFRCount],chanNumAdj(pulLocV),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV).meanFRCount]) median([allBaseline(~pulLocV).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV).meanFRCount]) median([allBaseline(pulLocV).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

figure
subplot(211)
plot(firingRate(~pulLocV),chanNumAdj(~pulLocV)+32,'k*')
hold on
plot(firingRate(pulLocV),chanNumAdj(pulLocV),'go')
title('Average FR across task')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
abHistD = histogram(firingRate(~pulLocV),linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram(firingRate(pulLocV),linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
plot([median(firingRate(~pulLocV)) median(firingRate(~pulLocV))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median(firingRate(pulLocV)) median(firingRate(pulLocV))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

%%
figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV & visResp).meanFRCount],chanNumAdj(~pulLocV & visResp)+32,'k*')
hold on
plot([allBaseline(pulLocV & visResp).meanFRCount],chanNumAdj(pulLocV & visResp),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV & visResp).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV & visResp).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV & visResp).meanFRCount]) median([allBaseline(~pulLocV & visResp).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV & visResp).meanFRCount]) median([allBaseline(pulLocV & visResp).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

figure
subplot(211)
plot(firingRate(~pulLocV & visResp),chanNumAdj(~pulLocV & visResp)+32,'k*')
hold on
plot(firingRate(pulLocV & visResp),chanNumAdj(pulLocV & visResp),'go')
title('Average FR across task')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
abHistD = histogram(firingRate(~pulLocV & visResp),linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram(firingRate(pulLocV & visResp),linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
plot([median(firingRate(~pulLocV & visResp)) median(firingRate(~pulLocV & visResp))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median(firingRate(pulLocV & visResp)) median(firingRate(pulLocV & visResp))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

%%
figure
subplot(211)
% plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
plot([allBaseline(~pulLocV & attMod).meanFRCount],chanNumAdj(~pulLocV & attMod)+32,'k*')
hold on
plot([allBaseline(pulLocV & attMod).meanFRCount],chanNumAdj(pulLocV & attMod),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
% abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
abHistD = histogram([allBaseline(~pulLocV & attMod).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram([allBaseline(pulLocV & attMod).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
% plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
%     [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(~pulLocV & attMod).meanFRCount]) median([allBaseline(~pulLocV & attMod).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median([allBaseline(pulLocV & attMod).meanFRCount]) median([allBaseline(pulLocV & attMod).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')

figure
subplot(211)
plot(firingRate(~pulLocV & attMod),chanNumAdj(~pulLocV & attMod)+32,'k*')
hold on
plot(firingRate(pulLocV & attMod),chanNumAdj(pulLocV & attMod),'go')
title('Average FR across task')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
abHistD = histogram(firingRate(~pulLocV & attMod),linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3,'normalization','Probability');
hold on
abHistV = histogram(firingRate(pulLocV & attMod),linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6,'normalization','Probability');
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
plot([median(firingRate(~pulLocV & attMod)) median(firingRate(~pulLocV & attMod))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'k')
plot([median(firingRate(pulLocV & attMod)) median(firingRate(pulLocV & attMod))],...
    [0 max([abHistD.Values abHistV.Values])+0.02],'g')
ylabel('Probability')
xlabel('Firing rate in Hz')





% eof