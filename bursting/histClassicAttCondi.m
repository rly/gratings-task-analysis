function [firingRatefun, percClassicTwoSpikeBurstfun, percClassicThreeSpikeBurstfun, numInterSpikeTimesfun, numBurstLongfun,...
    numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
    percClassicTwoSpikeBurstAttInfun, percClassicTwoSpikeBurstAttOutfun] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut, tempResolve)

            
% [firingRatefun, percClassicTwoSpikeBurstfun, percClassicThreeSpikeBurstfun, numInterSpikeTimesfun, numBurstLongfun,...
%     numBurstShortfun, frAttInPerTrialfun,frAttOutPerTrialfun,...
%     percClassicTwoSpikeBurstAttInfun, percClassicThreeSpikeBurstAttInfun,...
%     numInterSpikeTimesAttInfun, numBurstLongAttInfun, numBurstShortAttInfun,...
%     percClassicTwoSpikeBurstAttOutfun, percClassicThreeSpikeBurstAttOutfun,...
%     numInterSpikeTimesAttOutfun, numBurstLongAttOutfun, numBurstShortAttOutfun] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
%                 endTimeAttIn, endTimeAttOut, tempResolve)
%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;

% calculate traditional bursting in all spike times
interSpikeTimes = diff(spikeTimes2use);
interSpikeTimes = interSpikeTimes*1000; % in ms
firingRatefun = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

% indices of 100ms quiet periods
quietInd = find(interSpikeTimes>100);
% find times after 100ms quiet that are <= to 4ms
if quietInd(end) == numel(interSpikeTimes); quietInd = quietInd(1:end-1); end
burstInd = find(interSpikeTimes(quietInd+1)<=4);
percClassicTwoSpikeBurstfun = numel(interSpikeTimes(quietInd(burstInd))) * 100 / numel(interSpikeTimes);

if ~isempty(burstInd)
    if quietInd(burstInd(end))+2 >= numel(interSpikeTimes); burstInd = burstInd(1:end-1); end
    burstIndLong = find(interSpikeTimes(quietInd(burstInd)+2)<=4);
    percClassicThreeSpikeBurstfun = numel(burstIndLong) * 100 / numel(interSpikeTimes);
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = numel(burstIndLong);
    numBurstShortfun = numel(burstInd);
else
    percClassicThreeSpikeBurstfun = 0;
    numInterSpikeTimesfun = numel(interSpikeTimes);
    numBurstLongfun = 0;
    numBurstShortfun = 0;
end

% repeat for the task times
data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

kernelSigma = .01;
if tempResolve
%         [cueOnset.spikeTimes,cueOnset.spikeIndices] = createnonemptydatamatpt(spikeTimes2use, ...
%             floor(cueOnsetTrialAttIn*Fs)+1, cueOnset.window);

%     cueOnset.window = [0.8 1.2]; % seconds before, after
%     cueOnset.spdfWindowOffset = [-0.7 1.1];
%     [cueOnset.spikeTimes,cueOnset.spikeIndices] = createnonemptydatamatpt(spikeTimes2use(quietInd(burstInd)+1), ...
%             cueOnsetTrialAttIn, cueOnset.window);
%     cueOnset.t = computeTForSpdf(cueOnset.window(1), cueOnset.spdfWindowOffset, kernelSigma);
%     [cueOnset.spdf,~,cueOnset.spdfErr,cueOnset.singleTrialSpdf] = fixedPsth(cueOnset.spikeTimes, kernelSigma, 2, cueOnset.t);
end

% % binarySpikeTrain = zeros(1,length(data_pts));
% % binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,4),.0000000582)) = 1;
% % sum(binarySpikeTrain) 
% % numel(spikeTimes2use)
% 
% if numel(spikeTimes2use)==sum(binarySpikeTrain)
%     discrep = -1;
% else
%     discrep = 1;
% end

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
        isiPerTrialAttIn(n).isi = diff(find(data(indx)));
        nSpikesAttInPerTrial(n) = sum(data(indx));
        frAttInPerTrialfun(n) = nSpikesAttInPerTrial(n) / ((indx(end) - indx(1))/Fs);
%         attIn(n).spikeTimes=find(data(indx));
        clear indx
    end
end
for n=1:NEAttOut;
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        isiPerTrialAttOut(n).isi = diff(find(data(indx)));
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        frAttOutPerTrialfun(n) = nSpikesAttOutPerTrial(n) / ((indx(end) - indx(1))/Fs);
%         attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end

% for att in
% per trial
burstIndAttInPerTrial = [];
for n = 1:numel(isiPerTrialAttIn)
    qPeriod = isiPerTrialAttIn(n).isi>100;
    if sum(qPeriod)>=1
        if sum(find(qPeriod)+1<=numel(isiPerTrialAttIn(n).isi)) < sum(qPeriod)
            qPeriod = qPeriod(1:end-1);
        end
        if isiPerTrialAttIn(n).isi(find(qPeriod)+1) <= 4
            burstIndAttInTmp = [find(qPeriod) find(qPeriod)+1];
            burstIndAttInPerTrial = [burstIndAttInPerTrial burstIndAttInTmp];
        end
    end
end
if ~isempty(burstIndAttInPerTrial)
    percClassicTwoSpikeBurstAttInfun = numel(burstIndAttInPerTrial) * 100 / sum(nSpikesAttInPerTrial);
else
    percClassicTwoSpikeBurstAttInfun = 0;
end
      

% % indices of 100ms quiet periods
% quietIndAttIn = find(datatmpAttIn>100);
% % find times after 100ms quiet that are <= to 4ms
% if ~isempty(quietIndAttIn)
%     if quietIndAttIn(end) == numel(datatmpAttIn); quietIndAttIn = quietIndAttIn(1:end-1); end
%     burstIndAttIn = find(datatmpAttIn(quietIndAttIn+1)<=4);
%     percClassicTwoSpikeBurstAttInfun = numel(datatmpAttIn(quietIndAttIn(burstIndAttIn))) * 100 / sum(nSpikesAttInPerTrial);
% else
%     burstIndAttIn = [];
%     percClassicTwoSpikeBurstAttInfun = NaN;
% end

% if ~isempty(burstIndAttIn)
%     if quietIndAttIn(burstIndAttIn(end))+2 >= numel(datatmpAttIn); burstIndAttIn = burstIndAttIn(1:end-1); end
%     burstIndLongAttIn = find(datatmpAttIn(quietIndAttIn(burstIndAttIn)+2)<=4);
%     percClassicThreeSpikeBurstAttInfun = numel(burstIndLongAttIn) * 100 / sum(nSpikesAttInPerTrial);
%     numInterSpikeTimesAttInfun = numel(datatmpAttIn);
%     numBurstLongAttInfun = numel(burstIndLongAttIn);
%     numBurstShortAttInfun = numel(burstIndAttIn);
% else
%     percClassicThreeSpikeBurstAttInfun = 0;
%     numInterSpikeTimesAttInfun = numel(datatmpAttIn);
%     numBurstLongAttInfun = 0;
%     numBurstShortAttInfun = 0;
% end

% for att out
burstIndAttOutPerTrial = [];
for n = 1:numel(isiPerTrialAttOut)
    qPeriod = isiPerTrialAttOut(n).isi>100;
    if sum(qPeriod)>=1
        if sum(find(qPeriod)+1<=numel(isiPerTrialAttOut(n).isi)) < sum(qPeriod)
            qPeriod = qPeriod(1:end-1);
        end
        if isiPerTrialAttOut(n).isi(find(qPeriod)+1) <= 4
            burstIndAttOutTmp = [find(qPeriod) find(qPeriod)+1];
            burstIndAttOutPerTrial = [burstIndAttOutPerTrial burstIndAttOutTmp];
        end
    end
end
if ~isempty(burstIndAttOutPerTrial)
    percClassicTwoSpikeBurstAttOutfun = numel(burstIndAttOutPerTrial) * 100 / sum(nSpikesAttOutPerTrial);
else
    percClassicTwoSpikeBurstAttOutfun = 0;
end
% % indices of 100ms quiet periods
% quietIndAttOut = find(datatmpAttOut>100);
% % find times after 100ms quiet that are <= to 4ms
% if ~isempty(quietIndAttOut)
%     if quietIndAttOut(end) == numel(datatmpAttOut); quietIndAttOut = quietIndAttOut(1:end-1); end
%     burstIndAttOut = find(datatmpAttOut(quietIndAttOut+1)<=4);
%     percClassicTwoSpikeBurstAttOutfun = numel(datatmpAttOut(quietIndAttOut(burstIndAttOut))) * 100 / sum(nSpikesAttOutPerTrial);
% else
%     burstIndAttOut = [];
%     percClassicTwoSpikeBurstAttOutfun = NaN;
% end
% 
% if ~isempty(burstIndAttOut)
%     if quietIndAttOut(burstIndAttOut(end))+2 >= numel(datatmpAttOut); burstIndAttOut = burstIndAttOut(1:end-1); end
%     burstIndLongAttOut = find(datatmpAttOut(quietIndAttOut(burstIndAttOut)+2)<=4);
%     percClassicThreeSpikeBurstAttOutfun = numel(burstIndLongAttOut) * 100 / sum(nSpikesAttOutPerTrial);
%     numInterSpikeTimesAttOutfun = numel(datatmpAttOut);
%     numBurstLongAttOutfun = numel(burstIndLongAttOut);
%     numBurstShortAttOutfun = numel(burstIndAttOut);
% else
%     percClassicThreeSpikeBurstAttOutfun = 0;
%     numInterSpikeTimesAttOutfun = numel(datatmpAttOut);
%     numBurstLongAttOutfun = 0;
%     numBurstShortAttOutfun = 0;
% end
