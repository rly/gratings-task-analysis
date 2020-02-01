clc
clear all 
close all
cd('/Users/labmanager/Documents/MATLAB/Data locally')
nUnits = dir('spikeTimes2use*');

binWidth4hist = 5;
plotOutput = 1;

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

for uniti = 1:length(nUnits)
    if uniti == 6
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
    
    saveFileName = ['HistogramISIrealCueTargetDim' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
    [percentageBurst(uniti), percBurstWMAttCA(uniti), percBurstWMPoissonCA(uniti), cueArrayAttIn(uniti), cueArrayAttOut(uniti),...
    cueArrayshuffledAttIn(uniti), cueArrayshuffledAttOut(uniti),cueArraypoissonAttIn(uniti), cueArraypoissonAttOut(uniti)] = histWDattCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
        arrayOnsetTrialAttIn, arrayOnsetTrialAttOut, binWidth4hist, plotOutput,saveFileName);
    saveFileName = ['HistogramISIrealCueArrayOnset' nUnits(uniti).name(19:end-4) '_' num2str(binWidth4hist) 'ms.png'];
    [~, percBurstWMAttCT(uniti), percBurstWMPoissonCT(uniti), cueTargetAttIn(uniti), cueTargetAttOut(uniti),...
    cueTargetshuffledAttIn(uniti), cueTargetshuffledAttOut(uniti),cueTargetpoissonAttIn(uniti), cueTargetpoissonAttOut(uniti)] = histWDattCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
        targetDimTrialAttIn, targetDimTrialAttOut, binWidth4hist, plotOutput,saveFileName);
    clear UE endTime firstSpikeTimes spikeTimes2use startTime ...
        cueOnsetTrialAttIn cueOnsetTrialAttOut arrayOnsetTrialAttIn ...
        arrayOnsetTrialAttOut targetDimTrialAttIn targetDimTrialAttOut
    close all
end

cd('/Users/labmanager/Documents/MATLAB/BurstSep4all')
save(['AttNw3_' num2str(binWidth4hist) 'ms'],'cueArray*',...
    'cueTarget*','perc*')

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





% eof