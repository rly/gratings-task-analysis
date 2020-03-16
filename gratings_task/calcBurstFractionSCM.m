function [percentageBurst,percBurstWMAttWhole,percBurstWMAttIn,percBurstWMAttOut,...
    percBurstWMAttShuffled,percBurstWMPoisson,histBased,...
    SPNAAttIn,SPNAAttOut,BRIAttIn,BRIAttOut] = calcBurstFractionSCM(epochedSpikeTimesAttIn,...
    epochedSpikeTimesAttOut,poissonBinaryAttIn,poissonBinaryAttOut,datatmpShuffledAttIn,...
    datatmpShuffledAttOut,spikesInTrialsAttIn,spikesInTrialsAttOut,...
    binWidth4hist,maxLag,Poisson,Shuffled)


% FOR BURST FRACTION
spikeTimes2useAttIn = [];
datatmpAttIn = []; datatmpAttInPoisson = []; lastSpikeTime = 0;
for n = 1:size(epochedSpikeTimesAttIn,2)
    spikeTimes2useAttIn = [spikeTimes2useAttIn (epochedSpikeTimesAttIn(n).spikeTimes+lastSpikeTime)];
    if ~isempty(epochedSpikeTimesAttIn(n).spikeTimes)
        lastSpikeTime = spikeTimes2useAttIn(end);
    end
    datatmpAttIn = [datatmpAttIn diff(find(epochedSpikeTimesAttIn(n).binData))];
    if Poisson
        datatmpAttInPoisson = [datatmpAttInPoisson diff(find(poissonBinaryAttIn(n).binData))];
    end
end

spikeTimes2useAttOut = [];
datatmpAttOut = []; datatmpAttOutPoisson = []; lastSpikeTime = 0;
for n = 1:size(epochedSpikeTimesAttOut,2)
    spikeTimes2useAttOut = [spikeTimes2useAttOut (epochedSpikeTimesAttOut(n).spikeTimes+lastSpikeTime)];
    if ~isempty(epochedSpikeTimesAttOut(n).spikeTimes)
        lastSpikeTime = spikeTimes2useAttOut(end);
    end
    datatmpAttOut = [datatmpAttOut diff(find(epochedSpikeTimesAttOut(n).binData))];
    if Poisson
        datatmpAttOutPoisson = [datatmpAttOutPoisson diff(find(poissonBinaryAttOut(n).binData))];
    end
end

spikeTimes2use = [spikeTimes2useAttIn (spikeTimes2useAttOut+spikeTimes2useAttIn(end))];
% calculate percentage of spikes >= 200 Hz firing rate
percentageBurst =  sum(diff(spikeTimes2use)<=5) * 100 / size(spikeTimes2use,2);

totalAttConds = [datatmpAttIn datatmpAttOut];
percBurstWMAttWhole.four = sum(totalAttConds <= 4) *100 / size(totalAttConds,2);
percBurstWMAttWhole.five = sum(totalAttConds <= 5) *100 / size(totalAttConds,2);
percBurstWMAttWhole.ten = sum(totalAttConds <= 10) *100 / size(totalAttConds,2);
percBurstWMAttIn.four = sum(datatmpAttIn <= 4) *100 / size(datatmpAttIn,2);
percBurstWMAttIn.five = sum(datatmpAttIn <= 5) *100 / size(datatmpAttIn,2);
percBurstWMAttIn.ten = sum(datatmpAttIn <= 10) *100 / size(datatmpAttIn,2);
percBurstWMAttOut.four = sum(datatmpAttOut <= 4) *100 / size(datatmpAttOut,2);
percBurstWMAttOut.five = sum(datatmpAttOut <= 5) *100 / size(datatmpAttOut,2);
percBurstWMAttOut.ten = sum(datatmpAttOut <= 10) *100 / size(datatmpAttOut,2);

% HISTOGRAM BASED METHOD
attIn4stats = histogram(datatmpAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
histBased.realAttIn = attIn4stats.Values(1);
attOut4stats = histogram(datatmpAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
histBased.realAttOut = attOut4stats.Values(1);

if Shuffled
    totalAttShuffledConds = [datatmpShuffledAttIn datatmpShuffledAttOut];
    percBurstWMAttShuffled = sum(totalAttShuffledConds <= 5) * 100 / size(totalAttShuffledConds,2);
    attIn4statsSHUF = histogram(datatmpShuffledAttIn,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    histBased.shuffledAttIn = attIn4statsSHUF.Values(1);
    attOut4statsSHUF = histogram(datatmpShuffledAttOut,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    histBased.shuffledAttOut = attOut4statsSHUF.Values(1);
else
    percBurstWMAttShuffled = 0;
end
if Poisson
    datatmpPoisson = [datatmpAttInPoisson datatmpAttOutPoisson];
    percBurstWMPoisson = sum(datatmpPoisson <= 5) * 100 / size(datatmpPoisson,2);
    attIn4statsPoisson = histogram(datatmpAttInPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    histBased.poissonAttIn = attIn4statsPoisson.Values(1);
    attOut4statsPoisson = histogram(datatmpAttOutPoisson,'BinWidth',binWidth4hist,'BinLimits',[0 800],'Normalization','probability');
    histBased.poissonAttOut = attOut4statsPoisson.Values(1);
else
    percBurstWMPoisson = 0;
end

% AUTOCORRELATION BASED METHOD (BRI)
[SPNAAttIn,~,~,~] = computeSPNA(spikesInTrialsAttIn, maxLag);
[SPNAAttOut,~,~,~] = computeSPNA(spikesInTrialsAttOut, maxLag);

BRILagStartInd = maxLag + 1 + 1000 * 0.001;
BRILagEndInd = maxLag + 1 + 1000 * 0.004;
BRIAttIn = mean(SPNAAttIn(BRILagStartInd:BRILagEndInd));
BRIAttOut = mean(SPNAAttOut(BRILagStartInd:BRILagEndInd));

end