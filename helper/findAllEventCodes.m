function [eventCodesByPort,allEventTimesCombined] = findAllEventCodes(allEventTimes)
% allEventTimes - N x 1 cell matrix for event times of EVT1 to EVTN
% N must be a multiple of 8

allEventTimesCombined = unique(round(cell2mat(allEventTimes'), 4));

% tolerance level for simultaneous signals: if two events are within this
% number of seconds from each other, they represent the same signal.
simultSignalTol = 0.001;

nChannel = numel(allEventTimes);
nPort = nChannel / 8;
assert(mod(nChannel, 8) == 0);

% initialize
eventCodesByPort = cell(nPort, 1);
for j = 1:nPort
    eventCodesByPort{j} = nan(numel(allEventTimesCombined), 1);
end

for i = 1:numel(allEventTimesCombined)
    isFoundEvt = false(nChannel, 1);
    foundEvts = cell(nChannel, 1);
    foundEvtTimes = cell(nChannel, 1);
    
    % find the event on each channel that is at the same time (or within
    % tolerance) of current event time
    for j = 1:nChannel
        foundEvts{j} = find(abs(allEventTimesCombined(i) - allEventTimes{j}) < simultSignalTol);
        foundEvtTimes{j} = allEventTimes{j}(foundEvts{j});
        assert(numel(foundEvts{j}) <= 1);
        isFoundEvt(j) = ~isempty(foundEvts{j});
    end
    assert(any(isFoundEvt));
    
    % ensure that all found event times are within the tolerance away from
    % each other
    sortedFoundEvtTimes = sort(cell2mat(foundEvtTimes));
    assert(sortedFoundEvtTimes(end) - sortedFoundEvtTimes(1) < simultSignalTol);
    
    % convert binary code across channels to decimal code (1 to 255)
    for j = 1:nPort
        % j == 1: ports are 1:8, j == 2: ports are 9:16
        startPort = (j-1)*8+1;
        endPort = j*8;
        % convert from reversed binary
        % if event 8 is 1 and events 1-7 are 0, then code is 128
        eventCodesByPort{j}(i) = bin2dec(num2str(isFoundEvt(endPort:-1:startPort)'));
    end
end