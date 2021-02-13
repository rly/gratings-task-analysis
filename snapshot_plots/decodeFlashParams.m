function flashParams = decodeFlashParams(flashEventTimes, ...
        allEventTimes, distsToFix, polarAngles, gratingAngles)

simultSignalTol = 0.001; % tolerance level for simultaneous signals --
% if two events are within this number of  seconds from each other, they 
% represent the same signal. 

theta = nan(size(flashEventTimes));
distance = nan(size(flashEventTimes));
 % (x, y, grating angle, distsToFix index, polarAngles index, gratingAngles index)
flashParams = nan(numel(flashEventTimes), 6);
flashOnsets = nan(size(flashEventTimes));

for i = 1:numel(flashEventTimes)
    % convert binary representation of grating location and orientation 
    % indices to decimal. 
    distBinaryCode = 0;
    polarAngleBinaryCode = 0;
    gratingAngleBinaryCode = 0;
    
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
    
    if numel(foundEvts{9}) == 1
        distBinaryCode = distBinaryCode + 1;
    end
    if numel(foundEvts{10}) == 1
        distBinaryCode = distBinaryCode + 2;
    end
    if numel(foundEvts{11}) == 1
        polarAngleBinaryCode = polarAngleBinaryCode + 1;
    end
    if numel(foundEvts{12}) == 1
        polarAngleBinaryCode = polarAngleBinaryCode + 2;
    end
    if numel(foundEvts{13}) == 1
        polarAngleBinaryCode = polarAngleBinaryCode + 4;
    end
    if numel(foundEvts{14}) == 1
        polarAngleBinaryCode = polarAngleBinaryCode + 8;
    end
    if numel(foundEvts{15}) == 1
        gratingAngleBinaryCode = gratingAngleBinaryCode + 1;
    end
    if numel(foundEvts{16}) == 1
        gratingAngleBinaryCode = gratingAngleBinaryCode + 2;
    end
    
    % there is always a coded polar angle for every flash. 
    % there is not always a coded distance or grating angle (i.e. may be 0).
    if polarAngleBinaryCode ~= 0
        distance(i) = distsToFix(distBinaryCode + 1);
        theta(i) = polarAngles(polarAngleBinaryCode);
        flashParams(i,1:2) = distance(i) * [cos(theta(i)) sin(theta(i))];
        flashParams(i,3) = gratingAngles(gratingAngleBinaryCode + 1);
        flashParams(i,4:6) = [distBinaryCode + 1 polarAngleBinaryCode gratingAngleBinaryCode + 1];
    end
end

assert(~any(isnan(flashParams(:))));
assert(sum(isnan(flashOnsets)) == 0);
assert(all((flashOnsets - flashEventTimes) < simultSignalTol));
