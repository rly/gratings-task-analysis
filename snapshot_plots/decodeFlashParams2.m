function [flashOnsets,decimalCodes] = decodeFlashParams2(flashEventTimes, ...
        allEventTimes)

simultSignalTol = 0.001; % tolerance level for simultaneous signals --
% if two events are within this number of  seconds from each other, they 
% represent the same signal. 

flashOnsets = nan(size(flashEventTimes));
decimalCodes = nan(size(flashEventTimes));

for i = 1:numel(flashEventTimes)   
    isFoundEvt = 0;
    foundEvts = cell(16, 1);
    for j = 9:16
        foundEvts{j} = find(abs(flashEventTimes(i) - allEventTimes{j}) < simultSignalTol);
        % there should be at most 1 event9, event10, etc within the
        % simultaneous signal tolerance of each flashEventTime
        assert(numel(foundEvts{j}) <= 1)
        if ~isempty(foundEvts{j})
            isFoundEvt = 1;
        end
    end
    if ~isFoundEvt
        fprintf('huh...\n');
    end

    allFoundEvtTimes = [allEventTimes{9}(foundEvts{9}); 
            allEventTimes{10}(foundEvts{10}); 
            allEventTimes{11}(foundEvts{11}); 
            allEventTimes{12}(foundEvts{12}); 
            allEventTimes{13}(foundEvts{13});
            allEventTimes{14}(foundEvts{14}); 
            allEventTimes{15}(foundEvts{15}); 
            allEventTimes{16}(foundEvts{16})];
    if (isempty(allFoundEvtTimes))
        continue;
    end;
    % TODO: do pairwise comparisons instead of consecutive comparisons
    % but this should be good enough
    assert(all(diff(allFoundEvtTimes) < simultSignalTol));
    flashOnsets(i) = allFoundEvtTimes(1);
    
    foundEvtsLogical = cellfun(@any, foundEvts);
    
    binaryCode = foundEvtsLogical(9:16)';
    decimalCodes(i) = bi2de(binaryCode);
end

assert(sum(isnan(flashOnsets)) == 0);
assert(sum(isnan(decimalCodes)) == 0);
