clc
clear all 
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally')
nUnits = dir('spikeTimes2use*');
javaaddpath(which('MatlabGarbageCollector.jar'))
binWidth4hist = 5;
plotOutput = 0;

for stati = 1:50%:100
    jheapcl
    tic
    percentageBurst= nan(length(nUnits),1);
    percBurstWMAttCA= nan(length(nUnits),1);
    percBurstWMPoissonCA= nan(length(nUnits),1);
    cueArrayAttIn= nan(length(nUnits),1);
    cueArrayAttOut= nan(length(nUnits),1);
    cueArrayshuffledAttIn= nan(length(nUnits),1);
    cueArrayshuffledAttOut= nan(length(nUnits),1);
    cueArraypoissonAttIn= nan(length(nUnits),1);
    cueArraypoissonAttOut= nan(length(nUnits),1);

    percBurstWMAttCT= nan(length(nUnits),1);
    percBurstWMPoissonCT= nan(length(nUnits),1);
    cueTargetAttIn= nan(length(nUnits),1);
    cueTargetAttOut= nan(length(nUnits),1);
    cueTargetshuffledAttIn= nan(length(nUnits),1);
    cueTargetshuffledAttOut= nan(length(nUnits),1);
    cueTargetpoissonAttIn= nan(length(nUnits),1);
    cueTargetpoissonAttOut= nan(length(nUnits),1);

    percentageBurstSCM= nan(length(nUnits),1);
    percBurstWMAttCASCM= nan(length(nUnits),1);
    percBurstWMPoissonCASCM= nan(length(nUnits),1);
    cueArrayAttInSCM= nan(length(nUnits),1);
    cueArrayAttOutSCM= nan(length(nUnits),1);
    cueArrayshuffledAttInSCM= nan(length(nUnits),1);
    cueArrayshuffledAttOutSCM= nan(length(nUnits),1);
    cueArraypoissonAttInSCM= nan(length(nUnits),1);
    cueArraypoissonAttOutSCM= nan(length(nUnits),1);

    percBurstWMAttCTSCM= nan(length(nUnits),1);
    percBurstWMPoissonCTSCM= nan(length(nUnits),1);
    cueTargetAttInSCM= nan(length(nUnits),1);
    cueTargetAttOutSCM= nan(length(nUnits),1);
    cueTargetshuffledAttInSCM= nan(length(nUnits),1);
    cueTargetshuffledAttOutSCM= nan(length(nUnits),1);
    cueTargetpoissonAttInSCM= nan(length(nUnits),1);
    cueTargetpoissonAttOutSCM= nan(length(nUnits),1);
    for uniti = 16:length(nUnits)
        if uniti == 21
            continue
        end
        cd('/Users/labmanager/Documents/MATLAB/Data locally')
        load(nUnits(uniti).name)
        cd('/Users/labmanager/Documents/MATLAB/gratings-task-analysis/bursting')

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

        saveFileName = ['HistogramISIrealCueArrayDim' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
        [percentageBurst(uniti), percBurstWMAttCA(uniti), percBurstWMPoissonCA(uniti), cueArrayAttIn(uniti), cueArrayAttOut(uniti),...
        cueArrayshuffledAttIn(uniti), cueArrayshuffledAttOut(uniti),cueArraypoissonAttIn(uniti), cueArraypoissonAttOut(uniti)] = histWDattCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
            arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, binWidth4hist, plotOutput,saveFileName);
        saveFileName = ['HistogramISIrealCueTargetOnset' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
        [~, percBurstWMAttCT(uniti), percBurstWMPoissonCT(uniti), cueTargetAttIn(uniti), cueTargetAttOut(uniti),...
        cueTargetshuffledAttIn(uniti), cueTargetshuffledAttOut(uniti),cueTargetpoissonAttIn(uniti), cueTargetpoissonAttOut(uniti)] = histWDattCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
            targetDimTrialAttIn, targetDimTrialAttOut, binWidth4hist, plotOutput,saveFileName);

        saveFileName = ['HistogramISIrealCueArrayDimSCM' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
        [percentageBurstSCM(uniti), percBurstWMAttCASCM(uniti), percBurstWMPoissonCASCM(uniti), cueArrayAttInSCM(uniti), cueArrayAttOutSCM(uniti),...
        cueArrayshuffledAttInSCM(uniti), cueArrayshuffledAttOutSCM(uniti),cueArraypoissonAttInSCM(uniti), cueArraypoissonAttOutSCM(uniti)] = histWDattCondiSCM(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
            arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, binWidth4hist, plotOutput,saveFileName);
        saveFileName = ['HistogramISIrealCueTargetOnsetSCM' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
        [~, percBurstWMAttCTSCM(uniti), percBurstWMPoissonCTSCM(uniti), cueTargetAttInSCM(uniti), cueTargetAttOutSCM(uniti),...
        cueTargetshuffledAttInSCM(uniti), cueTargetshuffledAttOutSCM(uniti),cueTargetpoissonAttInSCM(uniti), cueTargetpoissonAttOutSCM(uniti)] = histWDattCondiSCM(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
            targetDimTrialAttIn, targetDimTrialAttOut, binWidth4hist, plotOutput,saveFileName);
        clear UE endTime firstSpikeTimes spikeTimes2use startTime ...
            cueOnsetTrialAttIn cueOnsetTrialAttOut arrayOnsetTrialAttIn ...
            arrayOnsetTrialAttOut targetDimTrialAttIn targetDimTrialAttOut
        close all
    end


cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
save(['StatRunMcCartney' num2str(stati) '_AttNw3_' num2str(binWidth4hist) 'ms'],'cueArray*',...
    'cueTarget*','perc*')

% look at first bin of 200ms post-cue until array onset
cueArrayDiff = cueArrayAttIn(percBurstWMAttCA>1) - cueArrayAttOut(percBurstWMAttCA>1);
cueArrayDiffPoisson = cueArraypoissonAttIn(percBurstWMAttCA>1) - cueArraypoissonAttOut(percBurstWMAttCA>1);
cueArrayDiffShuffled = cueArrayshuffledAttIn(percBurstWMAttCA>1) - cueArrayshuffledAttOut(percBurstWMAttCA>1);
[~,PCAP(stati),~,~] = ttest(cueArrayDiff,cueArrayDiffPoisson);
[~,PCAS(stati),~,~] = ttest(cueArrayDiff,cueArrayDiffShuffled);
[~,PCA(stati),~,~] = ttest(cueArrayDiff);

% repeat for SCM
cueArrayDiffSCM = cueArrayAttInSCM(percBurstWMAttCASCM>1) - cueArrayAttOutSCM(percBurstWMAttCASCM>1);
cueArrayDiffPoissonSCM = cueArraypoissonAttInSCM(percBurstWMAttCASCM>1) - cueArraypoissonAttOutSCM(percBurstWMAttCASCM>1);
cueArrayDiffShuffledSCM = cueArrayshuffledAttInSCM(percBurstWMAttCASCM>1) - cueArrayshuffledAttOutSCM(percBurstWMAttCASCM>1);
[~,PCAPSCM(stati),~,~] = ttest(cueArrayDiffSCM,cueArrayDiffPoissonSCM);
[~,PCASSCM(stati),~,~] = ttest(cueArrayDiffSCM,cueArrayDiffShuffledSCM);
[~,PCASCM(stati),~,~] = ttest(cueArrayDiffSCM);

% look at first bin of 200ms post-cue until array onset
cueTargetDiff = cueTargetAttIn(percBurstWMAttCT>1) - cueTargetAttOut(percBurstWMAttCT>1);
cueTargetDiffPoisson = cueTargetpoissonAttIn(percBurstWMAttCT>1) - cueTargetpoissonAttOut(percBurstWMAttCT>1);
cueTargetDiffShuffled = cueTargetshuffledAttIn(percBurstWMAttCT>1) - cueTargetshuffledAttOut(percBurstWMAttCT>1);
[~,PCTP(stati),~,~] = ttest(cueTargetDiff,cueTargetDiffPoisson);
[~,PCTS(stati),~,~] = ttest(cueTargetDiff,cueTargetDiffShuffled);
[~,PCT(stati),~,~] = ttest(cueTargetDiff);

% repeat for SCM
cueTargetDiffSCM = cueTargetAttInSCM(percBurstWMAttCTSCM>1) - cueTargetAttOutSCM(percBurstWMAttCTSCM>1);
cueTargetDiffPoissonSCM = cueTargetpoissonAttInSCM(percBurstWMAttCTSCM>1) - cueTargetpoissonAttOutSCM(percBurstWMAttCTSCM>1);
cueTargetDiffShuffledSCM = cueTargetshuffledAttInSCM(percBurstWMAttCTSCM>1) - cueTargetshuffledAttOutSCM(percBurstWMAttCTSCM>1);
[~,PCTPSCM(stati),~,~] = ttest(cueTargetDiffSCM,cueTargetDiffPoissonSCM);
[~,PCTSSCM(stati),~,~] = ttest(cueTargetDiffSCM,cueTargetDiffShuffledSCM);
[~,PCTSCM(stati),~,~] = ttest(cueTargetDiffSCM);
toc
sprintf(['statRun: ' num2str(stati)])
save(['pValues_McCartney_' num2str(binWidth4hist) 'ms'],'PCAP','PCAS','PCA','PCAPSCM','PCASSCM','PCASCM','PCTP',...
    'PCTS','PCT','PCTPSCM','PCTSSCM','PCTSCM')
end

%%
binWidth4hist = 5;
cd('/Users/labmanager/Documents/MATLAB/BurstSep4all/data')
load(['pValues_McCartney_' num2str(binWidth4hist) 'ms'])


% cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
% data = [percentageBurst, percBurstWMAttCT, cueTargetAttDiff, percBurstWMAttCA,...
%         cueArrayAttDiff,cueTargetP, cueArrayP, percBurstWMPoissonCT,percBurstWMPoissonCA,...
%         cueTargetShuffledDiff, cueArrayShuffledDiff];
% dataTbl = array2table(data);dataTbl.Properties.VariableNames = {'wdAll' 'wdCueTarg'...
%  'cueTargDiff' 'wdCueArray' 'cueArrayDiff' 'cueTargetP' 'cueArrayP' 'wdPoissonCT' 'wdPoissonCA'};
% 
% save(['AttInMinusAttOutNw2_' num2str(binWidth4hist) 'ms'],'data','dataTbl')
% 
% figure
% plot(percentageBurst,'.')
% hold on
% plot(percBurstWMAttCT,'.')
% plot(percBurstWMAttCA,'.')
% plot(percBurstWMPoissonCT,'*')
% plot(percBurstWMPoissonCA,'*')


% [~,sortedBC] = sort(data(:,1))
% figure
% plot(data(sortedBC,1))
% hold on
% plot(data(sortedBC,2))
% plot(data(sortedBC,4))
% 
% nanmean(data(:,1))
% sum(data(:,1)>1)
% 
% nanmean(data(:,2))
% sum(data(:,2)>1)
% 
% nanmean(data(:,4))
% sum(data(:,4)>1)
% 
% nanmean(data(:,3))
% nanmean(data(:,5))

% [H,P,CI,STATS] = ttest(data(data(:,1)>1,3))
% [H,P,CI,STATS] = ttest(data(data(:,2)>1,3))
% 
% [H,P,CI,STATS] = ttest(data(data(:,1)>1,5))
% [H,P,CI,STATS] = ttest(data(data(:,4)>1,5))

%%
binWidth4hist = 10;
load(['AttNw3_' num2str(binWidth4hist) 'ms'])

% look at first bin of 200ms post-cue until array onset
cueArrayDiff = cueArrayAttIn(percBurstWMAttCA>1) - cueArrayAttOut(percBurstWMAttCA>1);
cueArrayDiffPoisson = cueArraypoissonAttIn(percBurstWMAttCA>1) - cueArraypoissonAttOut(percBurstWMAttCA>1);
cueArrayDiffShuffled = cueArrayshuffledAttIn(percBurstWMAttCA>1) - cueArrayshuffledAttOut(percBurstWMAttCA>1);
[H,P,CI,STATS] = ttest(cueArrayDiff,cueArrayDiffPoisson)
[H,P,CI,STATS] = ttest(cueArrayDiff,cueArrayDiffShuffled)
[H,P,CI,STATS] = ttest(cueArrayDiff)

% repeat for SCM
cueArrayDiffSCM = cueArrayAttInSCM(percBurstWMAttCASCM>1) - cueArrayAttOutSCM(percBurstWMAttCASCM>1);
cueArrayDiffPoissonSCM = cueArraypoissonAttInSCM(percBurstWMAttCASCM>1) - cueArraypoissonAttOutSCM(percBurstWMAttCASCM>1);
cueArrayDiffShuffledSCM = cueArrayshuffledAttInSCM(percBurstWMAttCASCM>1) - cueArrayshuffledAttOutSCM(percBurstWMAttCASCM>1);
[H,P,CI,STATS] = ttest(cueArrayDiffSCM,cueArrayDiffPoissonSCM)
[H,P,CI,STATS] = ttest(cueArrayDiffSCM,cueArrayDiffShuffledSCM)
[H,P,CI,STATS] = ttest(cueArrayDiffSCM)

% look at first bin of 200ms post-cue until array onset
cueTargetDiff = cueTargetAttIn(percBurstWMAttCT>1) - cueTargetAttOut(percBurstWMAttCT>1);
cueTargetDiffPoisson = cueTargetpoissonAttIn(percBurstWMAttCT>1) - cueTargetpoissonAttOut(percBurstWMAttCT>1);
cueTargetDiffShuffled = cueTargetshuffledAttIn(percBurstWMAttCT>1) - cueTargetshuffledAttOut(percBurstWMAttCT>1);
[H,P,CI,STATS] = ttest(cueTargetDiff,cueTargetDiffPoisson)
[H,P,CI,STATS] = ttest(cueTargetDiff,cueTargetDiffShuffled)
[H,P,CI,STATS] = ttest(cueTargetDiff)

% repeat for SCM
cueTargetDiffSCM = cueTargetAttInSCM(percBurstWMAttCTSCM>1) - cueTargetAttOutSCM(percBurstWMAttCTSCM>1);
cueTargetDiffPoissonSCM = cueTargetpoissonAttInSCM(percBurstWMAttCTSCM>1) - cueTargetpoissonAttOutSCM(percBurstWMAttCTSCM>1);
cueTargetDiffShuffledSCM = cueTargetshuffledAttInSCM(percBurstWMAttCTSCM>1) - cueTargetshuffledAttOutSCM(percBurstWMAttCTSCM>1);
[H,P,CI,STATS] = ttest(cueTargetDiffSCM,cueTargetDiffPoissonSCM)
[H,P,CI,STATS] = ttest(cueTargetDiffSCM,cueTargetDiffShuffledSCM)
[H,P,CI,STATS] = ttest(cueTargetDiffSCM)



% eof