function [allEventTimesCombined,eventCodesPort1,eventCodesPort2] = findAllEventCodes(allEventTimes)
% allEventTimes - N x 1 cell matrix for event times of EVT1 to EVTN

allEventTimesCombined = unique(sort(round(cell2mat(allEventTimes'), 4)));

simultSignalTol = 0.001; % tolerance level for simultaneous signals --
% if two events are within this number of seconds from each other, they 
% represent the same signal. 

nChannels = numel(allEventTimes);
assert(nChannels == 16);

eventCodesPort1 = nan(numel(allEventTimesCombined), 1);
eventCodesPort2 = nan(numel(allEventTimesCombined), 1);

fprintf('Found %d events.\n', numel(allEventTimesCombined));
for i = 1:numel(allEventTimesCombined)
    isFoundEvt = false(nChannels, 1);
    foundEvts = cell(nChannels, 1);
    foundEvtTimes = cell(nChannels, 1);
    for j = 1:nChannels
        foundEvts{j} = find(abs(allEventTimesCombined(i) - allEventTimes{j}) < simultSignalTol);
        foundEvtTimes{j} = allEventTimes{j}(foundEvts{j});
        assert(numel(foundEvts{j}) <= 1)
        isFoundEvt(j) = ~isempty(foundEvts{j});
    end
    assert(any(isFoundEvt));
    
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(cell2mat(foundEvtTimes)) < simultSignalTol));
    
    eventCodesPort1(i) = bin2dec(num2str(isFoundEvt(8:-1:1)'));
    eventCodesPort2(i) = bin2dec(num2str(isFoundEvt(16:-1:9)'));
end