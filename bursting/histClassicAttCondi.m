function [firingRate, percClassicTwoSpikeBurst, percClassicThreeSpikeBurst, numInterSpikeTimes, numBurstLong,...
    numBurstShort, discrep, frAttInPerTrial,frAttOutPerTrial,...
    percClassicTwoSpikeBurstAttIn, percClassicThreeSpikeBurstAttIn,...
    numInterSpikeTimesAttIn, numBurstLongAttIn, numBurstShortAttIn,...
    percClassicTwoSpikeBurstAttOut, percClassicThreeSpikeBurstAttOut,...
    numInterSpikeTimesAttOut, numBurstLongAttOut, numBurstShortAttOut] = histClassicAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut)

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;

% calculate traditional bursting in all spike times
interSpikeTimes = diff(spikeTimes2use);
interSpikeTimes = interSpikeTimes*1000; % in ms
firingRate = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

% indices of 100ms quiet periods
quietInd = find(interSpikeTimes>100);
% find times after 100ms quiet that are <= to 4ms
if quietInd(end) == numel(interSpikeTimes); quietInd = quietInd(1:end-1); end
burstInd = find(interSpikeTimes(quietInd+1)<=4);
interSpikeTimeInd = 1:length(interSpikeTimes);
percClassicTwoSpikeBurst = numel(interSpikeTimeInd(quietInd(burstInd))) * 100 / numel(interSpikeTimes);

if ~isempty(burstInd)
    if quietInd(burstInd(end))+2 >= numel(interSpikeTimes); burstInd = burstInd(1:end-1); end
    burstIndLong = find(interSpikeTimes(quietInd(burstInd)+2)<=4);
    percClassicThreeSpikeBurst = numel(burstIndLong) * 100 / numel(interSpikeTimes);
    numInterSpikeTimes = numel(interSpikeTimes);
    numBurstLong = numel(burstIndLong);
    numBurstShort = numel(burstInd);
else
    percClassicThreeSpikeBurst = NaN;
    numInterSpikeTimes = numel(interSpikeTimes);
    numBurstLong = NaN;
    numBurstShort = NaN;
end

% repeat for the task times
data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

% binarySpikeTrain = zeros(1,length(data_pts));
% binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,4),.0000000582)) = 1;
% sum(binarySpikeTrain) 
% numel(spikeTimes2use)

if numel(spikeTimes2use)==sum(binarySpikeTrain)
    discrep = -1;
else
    discrep = 1;
end

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
        frAttInPerTrial(n) = nSpikesAttInPerTrial(n) / ((indx(end) - indx(1))/Fs);
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
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        frAttOutPerTrial(n) = nSpikesAttOutPerTrial(n) / ((indx(end) - indx(1))/Fs);
%         attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end

% for att in
% indices of 100ms quiet periods
quietIndAttIn = find(datatmpAttIn>100);
% find times after 100ms quiet that are <= to 4ms
if ~isempty(quietIndAttIn)
    if quietIndAttIn(end) == numel(datatmpAttIn); quietIndAttIn = quietIndAttIn(1:end-1); end
    burstIndAttIn = find(datatmpAttIn(quietIndAttIn+1)<=4);
    percClassicTwoSpikeBurstAttIn = numel(datatmpAttIn(quietIndAttIn(burstIndAttIn))) * 100 / numel(datatmpAttIn);
else
    burstIndAttIn = [];
    percClassicTwoSpikeBurstAttIn = NaN;
end

if ~isempty(burstIndAttIn)
    if quietIndAttIn(burstIndAttIn(end))+2 >= numel(datatmpAttIn); burstIndAttIn = burstIndAttIn(1:end-1); end
    burstIndLongAttIn = find(datatmpAttIn(quietIndAttIn(burstIndAttIn)+2)<=4);
    percClassicThreeSpikeBurstAttIn = numel(burstIndLongAttIn) * 100 / numel(datatmpAttIn);
    numInterSpikeTimesAttIn = numel(datatmpAttIn);
    numBurstLongAttIn = numel(burstIndLongAttIn);
    numBurstShortAttIn = numel(burstIndAttIn);
else
    percClassicThreeSpikeBurstAttIn = NaN;
    numInterSpikeTimesAttIn = numel(datatmpAttIn);
    numBurstLongAttIn = NaN;
    numBurstShortAttIn = NaN;
end

% for att in
% indices of 100ms quiet periods
quietIndAttOut = find(datatmpAttOut>100);
% find times after 100ms quiet that are <= to 4ms
if ~isempty(quietIndAttOut)
    if quietIndAttOut(end) == numel(datatmpAttOut); quietIndAttOut = quietIndAttOut(1:end-1); end
    burstIndAttOut = find(datatmpAttOut(quietIndAttOut+1)<=4);
    percClassicTwoSpikeBurstAttOut = numel(datatmpAttOut(quietIndAttOut(burstIndAttOut))) * 100 / numel(datatmpAttOut);
else
    burstIndAttOut = [];
    percClassicTwoSpikeBurstAttOut = NaN;
end

if ~isempty(burstIndAttOut)
    if quietIndAttOut(burstIndAttOut(end))+2 >= numel(datatmpAttOut); burstIndAttOut = burstIndAttOut(1:end-1); end
    burstIndLongAttOut = find(datatmpAttOut(quietIndAttOut(burstIndAttOut)+2)<=4);
    percClassicThreeSpikeBurstAttOut = numel(burstIndLongAttOut) * 100 / numel(datatmpAttOut);
    numInterSpikeTimesAttOut = numel(datatmpAttOut);
    numBurstLongAttOut = numel(burstIndLongAttOut);
    numBurstShortAttOut = numel(burstIndAttOut);
else
    percClassicThreeSpikeBurstAttOut = NaN;
    numInterSpikeTimesAttOut = numel(datatmpAttOut);
    numBurstLongAttOut = NaN;
    numBurstShortAttOut = NaN;
end
