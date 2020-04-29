clc
clear all
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
nUnitsDir = dir('spikeTimes2use_UE_v15*');

nUnits = length(nUnitsDir);
%cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
cd('/Volumes/HDD1/BurstSep4allHDD/')

for uniti = 1:nUnits
    tic
    %         if uniti == 21
    %             continue
    %         end
    cd('/Users/labmanager/Documents/MATLAB/Data locally SQ3')
    load(nUnitsDir(uniti).name)
    load(['spikeTimes2use_UE_' nUnitsDir(uniti).name(23:end)],'pulLoc')
    %     if pulLoc ~= 'u'

    if pulLoc == 'u'
        pulLocAll(uniti) = 0;
    else
        pulLocAll(uniti) = 1;
        if pulLoc == 'v'
            pulLocV(uniti) = 1;
        else
            pulLocV(uniti) = 0;
        end
    end
    
    chanNum(uniti) = spikeStruct.channelID;
    
    cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')
    startTime = spikeStruct.unitStartTime; endTime = spikeStruct.unitEndTime;
    firstSpikeTimes = spikeStruct.tsTask(1);
    
    cueOnset = UE.cueOnset;
    baselineOnset = UE.cueOnset - 0.300;
    cueOnset = cueOnset(cueOnset >= startTime & cueOnset <= endTime);
    cueOnset = cueOnset - firstSpikeTimes;
    baselineOnset = baselineOnset(baselineOnset >= startTime & baselineOnset <= endTime);
    baselineOnset = baselineOnset - firstSpikeTimes;
    
    cueOnset = cueOnset(baselineOnset>0);
    baselineOnset = baselineOnset(baselineOnset>0);
    
    spikeTimes2use = spikeStruct.tsTask;
    firingRate(uniti) = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

    spikeTimes2use = spikeTimes2use - firstSpikeTimes;
    % compute FR count
    for cuei = 1:numel(cueOnset)
        nwinr = cueOnset(cuei);
        nwinl = baselineOnset(cuei);
        baseline.spikeCount(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        baseline.frCount(cuei) = baseline.spikeCount(cuei) / round((nwinr - nwinl),3);
    end
     
    allBaseline(uniti).meanSpikeCount = mean(baseline.spikeCount);
    allBaseline(uniti).medianSpikeCount = median(baseline.spikeCount);
    allBaseline(uniti).meanFRCount = mean(baseline.frCount);
    allBaseline(uniti).medianFRCount = median(baseline.frCount);
    
    targetDimTrialAttOut = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 1);
    targetDimTrialAttOut = targetDimTrialAttOut(targetDimTrialAttOut >= startTime & targetDimTrialAttOut <= endTime);
    targetDimTrialAttOut = targetDimTrialAttOut - firstSpikeTimes;
    targetDimTrialAttIn = UE.targetDim(UE.cueLoc(UE.isHoldTrial) == 3);
    targetDimTrialAttIn = targetDimTrialAttIn(targetDimTrialAttIn >= startTime & targetDimTrialAttIn <= endTime);
    targetDimTrialAttIn = targetDimTrialAttIn - firstSpikeTimes;
    arrayOnsetTrialAttOut = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 1);
    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(arrayOnsetTrialAttOut >= startTime & arrayOnsetTrialAttOut <= endTime);
    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut - firstSpikeTimes;
    arrayOnsetTrialAttIn = UE.arrayOnset(UE.isHoldTrial & UE.cueLoc == 3);
    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(arrayOnsetTrialAttIn >= startTime & arrayOnsetTrialAttIn <= endTime);
    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn - firstSpikeTimes;
    cueOnsetTrialAttOut = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 1) + 0.200;
    cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut >= startTime & cueOnsetTrialAttOut <= endTime);
    cueOnsetTrialAttOut = cueOnsetTrialAttOut - firstSpikeTimes;
    cueOnsetTrialAttIn = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3) + 0.200;
    cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn >= startTime & cueOnsetTrialAttIn <= endTime);
    cueOnsetTrialAttIn = cueOnsetTrialAttIn - firstSpikeTimes;

    cueOnsetTrialInRF = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 3);
    cueOnsetTrialInRF = cueOnsetTrialInRF(cueOnsetTrialInRF >= startTime & cueOnsetTrialInRF <= endTime);
    cueOnsetTrialInRF = cueOnsetTrialInRF - firstSpikeTimes;
    cueOnsetTrialOutRF = UE.cueOnset(UE.isHoldTrial & UE.cueLoc == 4);
    cueOnsetTrialOutRF = cueOnsetTrialOutRF(cueOnsetTrialOutRF >= startTime & cueOnsetTrialOutRF <= endTime);
    cueOnsetTrialOutRF = cueOnsetTrialOutRF - firstSpikeTimes;

    targetDimTrialAttOut = targetDimTrialAttOut(cueOnsetTrialAttOut>0);
    targetDimTrialAttIn = targetDimTrialAttIn(cueOnsetTrialAttIn>0);
    arrayOnsetTrialAttOut = arrayOnsetTrialAttOut(cueOnsetTrialAttOut>0);
    arrayOnsetTrialAttIn = arrayOnsetTrialAttIn(cueOnsetTrialAttIn>0);
    cueOnsetTrialAttOut = cueOnsetTrialAttOut(cueOnsetTrialAttOut>0);
    cueOnsetTrialAttIn = cueOnsetTrialAttIn(cueOnsetTrialAttIn>0);
    cueOnsetTrialInRF = cueOnsetTrialInRF(cueOnsetTrialInRF>0);
    cueOnsetTrialOutRF = cueOnsetTrialOutRF(cueOnsetTrialOutRF>0);

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
    [pSCCIsum(uniti)] = ranksum(spikeCount.cueInRFBaseline,spikeCount.cueInRF);
    [pSCCOsum(uniti)] = ranksum(spikeCount.cueOutRFBaseline,spikeCount.cueOutRF);

    % compute attentional modulation
    for cuei = 1:numel(cueOnsetTrialAttIn)
        nwinr = arrayOnsetTrialAttIn(cuei) - 0.175;
        nwinl = cueOnsetTrialAttIn(cuei) + 0.200;
        spikeCount.attIn(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.attIn(cuei) = spikeCount.attIn(cuei) / round((nwinr - nwinl),3);
    end
    for cuei = 1:numel(cueOnsetTrialAttOut)
        nwinr = arrayOnsetTrialAttOut(cuei) - 0.175;
        nwinl = cueOnsetTrialAttOut(cuei) + 0.200;
        spikeCount.attOut(cuei) = sum(spikeTimes2use>nwinl & spikeTimes2use<nwinr);
        frCount.attOut(cuei) = spikeCount.attOut(cuei) / round((nwinr - nwinl),3);
    end   
    meanSCAttInCount(uniti) = mean(spikeCount.attIn);
    meanSCAttOutCount(uniti) = mean(spikeCount.attOut);
    [pSCATTsum(uniti)] = ranksum(spikeCount.attIn,spikeCount.attOut);
    
    fprintf(['Unit: ' num2str(uniti) ' from ' num2str(nUnits) ' total Units\n'])
    toc
end

visResp = (pSCCIsum <= .05 & meanSCInRFCount > meanSCInRFBlCount) | (pSCCOsum <= .05 & meanSCOutRFCount > meanSCOutRFBlCount);
attMod = (pSCATTsum <= .05 & meanSCAttInCount > meanSCAttOutCount);

pulLocAll = boolean(pulLocAll);
if numel(pulLocV) ~= numel(pulLocAll)
    pulLocV = [pulLocV zeros(1,length(pulLocAll)-length(pulLocV))];
end
pulLocV = boolean(pulLocV);
figure
histogram([allBaseline(pulLocAll).meanFRCount],linspace(0,20,50))

figure
histogram([allBaseline(pulLocV).meanFRCount],linspace(0,20,50))
hold on
histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50))

chanNumAdj = chanNum; chanNumAdj(chanNumAdj>32) = chanNumAdj(chanNumAdj>32)-32;
figure
subplot(211)
plot([allBaseline(pulLocAll & ~pulLocV).meanFRCount],chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
hold on
plot([allBaseline(pulLocV).meanFRCount],chanNumAdj(pulLocV),'go')
title('Baseline FR')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
abHistD = histogram([allBaseline(pulLocAll & ~pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
hold on
abHistV = histogram([allBaseline(pulLocV).meanFRCount],linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6);
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
plot([median([allBaseline(pulLocAll & ~pulLocV).meanFRCount]) median([allBaseline(pulLocAll & ~pulLocV).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median([allBaseline(pulLocV).meanFRCount]) median([allBaseline(pulLocV).meanFRCount])],...
    [0 max([abHistD.Values abHistV.Values])+2],'g')
ylabel('Count')
xlabel('Firing rate in Hz')

figure
subplot(211)
plot(firingRate(pulLocAll & ~pulLocV),chanNumAdj(pulLocAll & ~pulLocV)+32,'k*')
hold on
plot(firingRate(pulLocV),chanNumAdj(pulLocV),'go')
title('Average FR across task')
ylabel('Electrode depth')
xlabel('Firing rate in Hz')
subplot(212)
abHistD = histogram(firingRate(pulLocAll & ~pulLocV),linspace(0,20,50),'FaceColor',[0 0 0],'FaceAlpha',0.3);
hold on
abHistV = histogram(firingRate(pulLocV),linspace(0,20,50),'FaceColor',[0 1 0],'FaceAlpha',0.6);
legend({'dorsal Pul' 'ventral Pul'},'AutoUpdate','off')
plot([median(firingRate(pulLocAll & ~pulLocV)) median(firingRate(pulLocAll & ~pulLocV))],...
    [0 max([abHistD.Values abHistV.Values])+2],'k')
plot([median(firingRate(pulLocV)) median(firingRate(pulLocV))],...
    [0 max([abHistD.Values abHistV.Values])+2],'g')
ylabel('Count')
xlabel('Firing rate in Hz')







% eof