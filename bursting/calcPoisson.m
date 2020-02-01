function [poissonBinary,datatmpPoisson] = calcPoisson(nTrials,spikeTimesInEvents,eventEnd,eventStart,Fs,data)
% create Poisson distributed interspike intervals or binary data 

% for debugging:
%     spikeTimesInEvents = nEAttIn;
%     eventEnd = targetDimTrialAttIn;
%     eventStart = cueOnsetTrialAttIn;

poissonBinary = [];
datatmpPoisson = [];
for n=1:nTrials;
%     nwinl=round(0.200*Fs);
    nwinrAttIn=round(eventEnd(n)*Fs);
    indxAttIn=spikeTimesInEvents(n):nwinrAttIn-1;
    if indxAttIn(end) > size(data,2)
        indxAttIn = indxAttIn(1:sum(indxAttIn<=size(data,2)));
    end
    durationS = eventEnd(n) - eventStart(n); % length of simulation
    spikesPerS = sum(data(indxAttIn))/durationS; % avg firing rate
    timeStepS = 0.001;                  % 1 msec
    times = 0:timeStepS:durationS;	% a vector with each time step		
    vt = rand(size(times)); % use rand to set probobility of spike
    spikes = (spikesPerS*timeStepS) > vt;
    poissonBinary = [poissonBinary spikes];
    datatmpPoisson=[datatmpPoisson diff(find(spikes))];
end
