polarAngles = -9:9;
minDiff = 9;
nSequences = 1;

remainingVals = polarAngles;
sequence = nan(numel(remainingVals), 1);
sequence(1) = remainingVals(1);
remainingVals(1) = [];
possibleOpts = unique(remainingVals(abs(remainingVals - sequence(1)) >= minDiff));

fID = fopen(sprintf('sequence%dto%dAtLeastDiff%dn%d.txt', ...
        polarAngles(1), polarAngles(end), minDiff, numel(sequence)), 'w');

checkNextPos(remainingVals, possibleOpts, sequence, minDiff, fID, nSequences);

fclose(fID);