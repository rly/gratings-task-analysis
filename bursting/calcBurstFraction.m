function [percentageBurst,percBurstWMAttWhole,percBurstWMAttIn,percBurstWMAttOut] = calcBurstFraction(spikeTimes2use, cueOnsetTrialAttIn, cueOnsetTrialAttOut,...
                endTimeAttIn, endTimeAttOut)

%endTimeAttIn = arrayOnsetTrialAttIn; endTimeAttOut =arrayOnsetTrialAttOut;
%endTimeAttIn = targetDimTrialAttIn; endTimeAttOut =targetDimTrialAttOut;
% calculate percentage of spikes >= 200 Hz firing rate
percentageBurst =  sum(diff(spikeTimes2use)<(1/200)) * 100 / size(spikeTimes2use,2);

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
for n=1:NEAttIn
    %     nwinl=round(0.200*Fs);
    nwinr=round(endTimeAttIn(n)*Fs);
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
    nwinr=round(endTimeAttOut(n)*Fs);
    indx=nEAttOut(n):nwinr-1;
    if length(indx) >1 && indx(end) < size(data,2)
        datatmpAttOut=[datatmpAttOut diff(find(data(indx)))];
        nSpikesAttOutPerTrial(n) = sum(data(indx));
        attOut(n).spikeTimes=find(data(indx));
        clear indx
    end
end

totalAttConds = [datatmpAttIn datatmpAttOut];
percBurstWMAttWhole.four = sum(totalAttConds <= 4)  / size(totalAttConds,2);
percBurstWMAttWhole.five = sum(totalAttConds <= 5)  / size(totalAttConds,2);
percBurstWMAttWhole.ten = sum(totalAttConds <= 10)  / size(totalAttConds,2);
percBurstWMAttIn.four = sum(datatmpAttIn <= 4)  / size(datatmpAttIn,2);
percBurstWMAttIn.five = sum(datatmpAttIn <= 5)  / size(datatmpAttIn,2);
percBurstWMAttIn.ten = sum(datatmpAttIn <= 10)  / size(datatmpAttIn,2);
percBurstWMAttOut.four = sum(datatmpAttOut <= 4)  / size(datatmpAttOut,2);
percBurstWMAttOut.five = sum(datatmpAttOut <= 5)  / size(datatmpAttOut,2);
percBurstWMAttOut.ten = sum(datatmpAttOut <= 10)  / size(datatmpAttOut,2);