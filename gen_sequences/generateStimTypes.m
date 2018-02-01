function stimTypeSequence = generateStimTypes(eccs, polarAngleOpts, ...
        polarAngleSequencesFileName, polarAngleMultFactor, seqIndex, ...
        gratingAngles, contrasts, driftDirections, mode)

numStimTypes = numel(eccs) * numel(polarAngleOpts) * ...
        numel(gratingAngles) * numel(contrasts) * numel(driftDirections);
stimTypes = nan(numStimTypes, 7);

% get the polar angle sequences -- different one per row
fID = fopen(polarAngleSequencesFileName);
polarAngleSequences = fscanf(fID, '%d', [numStimTypes Inf])';
fclose(fID);
polarAngleSeq = polarAngleSequences(seqIndex,:)'; % pick line seqIndex

[sortedPolarAngleSeq,sortIndex] = sort(polarAngleSeq);
% randomize order of sort index within each group of polar angle
randSortIndex = nan(size(sortIndex));
numInGroup = numStimTypes/numel(polarAngleOpts);
for i = 1:numel(polarAngleOpts)
    starti = (i-1)*numInGroup + 1;
    endi = i*numInGroup;
    origOrder = sortIndex(starti:endi);
    randSortIndex(starti:endi) = origOrder(randperm(numInGroup));
end
revSortIndex(randSortIndex) = 1:numel(polarAngleSeq);

count = 0;
for b = 1:numel(polarAngleOpts) % in order
    for a = 1:numel(eccs)
        for c = 1:numel(gratingAngles)
            for d = 1:numel(contrasts)
                for e = 1:numel(driftDirections)
                    count = count + 1;
                    diameter = round((3.5 / 5 * eccs(a) - 2 * 28)*0.85);
                    if mode == 6 && a < 3 % make diameter same as for eccs(3) in this mode
                    	diameter = round((3.5 / 5 * eccs(3) - 2 * 28)*0.85);
                    end
                    stimTypes(count,:) = ...
                            [count eccs(a) polarAngleOpts(b) diameter ...
                            gratingAngles(c) contrasts(d) driftDirections(e)];
                end
            end
        end
    end
end

assert(~any(any(isnan(stimTypes))));
assert(all(sortedPolarAngleSeq == stimTypes(:,3)));
stimTypeSequence = stimTypes(revSortIndex,:);
stimTypeSequence(:,3) = stimTypeSequence(:,3) * polarAngleMultFactor;
