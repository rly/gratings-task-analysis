function [firingRate, percWDTwoSpikeBurst, percWDThreeSpikeBurst, numInterSpikeTimes, numBurstLong,...
    numBurstShort, frAttInPerTrial,frAttOutPerTrial,...
    percWDTwoSpikeBurstAttIn, percWDThreeSpikeBurstAttIn,...
    numInterSpikeTimesAttIn, numBurstLongAttIn, numBurstShortAttIn,...
    percWDTwoSpikeBurstAttOut, percWDThreeSpikeBurstAttOut,...
    numInterSpikeTimesAttOut, numBurstLongAttOut, numBurstShortAttOut] = histWomelsdorfAttCondi(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut)

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;

% calculate traditional bursting in all spike times
interSpikeTimes = diff(spikeTimes2use);
interSpikeTimes = interSpikeTimes*1000; % in ms
firingRate = numel(spikeTimes2use) / (spikeTimes2use(end) - spikeTimes2use(1));

burstInd = find(interSpikeTimes<=5);
percWDTwoSpikeBurst = numel(burstInd) * 100 / numel(interSpikeTimes);

if ~isempty(burstInd)
    if sum(burstInd+1 >= numel(interSpikeTimes)); burstInd = burstInd(1:end-1); end
    burstIndLong = find(interSpikeTimes(burstInd+1)<=5);
    percWDThreeSpikeBurst = numel(burstIndLong) * 100 / numel(interSpikeTimes);
    numInterSpikeTimes = numel(interSpikeTimes);
    numBurstLong = numel(burstIndLong);
    numBurstShort = numel(burstInd);
else
    percWDThreeSpikeBurst = 0;
    numInterSpikeTimes = numel(interSpikeTimes);
    numBurstLong = 0;
    numBurstShort = 0;
end

% repeat for the task times
data_pts = round(spikeTimes2use(1),3):1/1000:round(spikeTimes2use(end),3);
binarySpikeTrain = zeros(1,length(data_pts));
binarySpikeTrain(ismembertol(data_pts,round(spikeTimes2use,3),.00000001)) = 1;

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
burstIndAttIn = find(datatmpAttIn<=5);
percWDTwoSpikeBurstAttIn = numel(burstIndAttIn) * 100 / sum(nSpikesAttInPerTrial);

if ~isempty(burstIndAttIn)
    if sum(burstIndAttIn+1 >= numel(datatmpAttIn)); burstIndAttIn = burstIndAttIn(1:end-1); end
    burstIndLongAttIn = find(datatmpAttIn(burstIndAttIn+1)<=5);
    percWDThreeSpikeBurstAttIn = numel(burstIndLongAttIn) * 100 / sum(nSpikesAttInPerTrial);
    numInterSpikeTimesAttIn = numel(datatmpAttIn);
    numBurstLongAttIn = numel(burstIndLongAttIn);
    numBurstShortAttIn = numel(burstIndAttIn);
else
    percWDThreeSpikeBurstAttIn = 0;
    numInterSpikeTimesAttIn = numel(datatmpAttIn);
    numBurstLongAttIn = 0;
    numBurstShortAttIn = 0;
end

% for att out
burstIndAttOut = find(datatmpAttOut<=5);
percWDTwoSpikeBurstAttOut = numel(burstIndAttOut) * 100 / sum(nSpikesAttOutPerTrial);

if ~isempty(burstIndAttOut)
    if sum(burstIndAttOut+1 >= numel(datatmpAttOut)); burstIndAttOut = burstIndAttOut(1:end-1); end
    burstIndLongAttOut = find(datatmpAttOut(burstIndAttOut+1)<=5);
    percWDThreeSpikeBurstAttOut = numel(burstIndLongAttOut) * 100 / sum(nSpikesAttOutPerTrial);
    numInterSpikeTimesAttOut = numel(datatmpAttOut);
    numBurstLongAttOut = numel(burstIndLongAttOut);
    numBurstShortAttOut = numel(burstIndAttOut);
else
    percWDThreeSpikeBurstAttOut = 0;
    numInterSpikeTimesAttOut = numel(datatmpAttOut);
    numBurstLongAttOut = 0;
    numBurstShortAttOut = 0;
end