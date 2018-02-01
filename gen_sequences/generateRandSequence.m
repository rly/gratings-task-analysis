% generate a random sequence of N=4 subsequences where each subsequence is
% a random permutation of the numbers -8, -7, ..., 8 (polarAngles) where
% the difference between consecutive values is >= minDiff=7. The difference
% between the beginning and end of a subsequence also has to be >= 7, and
% the difference between the beginning of a subsequence and the end of the
% previous subsequence also has to be >= 7. The sequences of N=4
% subsequences and perhaps a shorter final sequence are output to the file
% randSequence-8to8AtLeastDiff7n17stitchn4.txt.

%% getting just one sequence of 19*4 is REALLY SLOW
% polarAngles = -9:9;
% minDiff = 8;
% nSequences = 1;
% 
% remainingVals = [polarAngles(randperm(numel(polarAngles))) ...
%         polarAngles(randperm(numel(polarAngles))) ...
%         polarAngles(randperm(numel(polarAngles))) ...
%         polarAngles(randperm(numel(polarAngles)))];
% sequence = nan(numel(remainingVals), 1);
% sequence(1) = remainingVals(1);
% remainingVals(1) = [];
% possibleOpts = unique(remainingVals(abs(remainingVals - sequence(1)) >= minDiff));
% fID = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%d.txt', ...
%         polarAngles(1), polarAngles(end), minDiff, numel(sequence)), 'a'); % append seqs
% checkNextPos(remainingVals, possibleOpts, sequence, minDiff, fID, nSequences);
% 
% fclose(fID);

%% better to run this a few times and append the appendable results
actualNSequences = 1000;
polarAngles = [-11 -8:8 11];
minDiff = 8;
nSequences = 1;
numValuesInSequences = numel(polarAngles);

%%
for i = 1:actualNSequences
    remainingVals = polarAngles(randperm(numel(polarAngles)));
    sequence = nan(numValuesInSequences, 1);
    sequence(1) = remainingVals(1);
    remainingVals(1) = [];
    possibleOpts = unique(remainingVals(abs(remainingVals - sequence(1)) >= minDiff));
    fID = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%d.txt', ...
            polarAngles(1), polarAngles(end), minDiff, numValuesInSequences), 'a'); % append seqs
    checkNextPos(remainingVals, possibleOpts, sequence, minDiff, fID, nSequences);

    fclose(fID);
end

%% read and append
maxSeqsPerLine = 3; % mode 1
fID = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%d.txt', ...
            polarAngles(1), polarAngles(end), minDiff, numValuesInSequences), 'r'); % append seqs
fIDLong = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%dstitchn%d.txt', ...
            polarAngles(1), polarAngles(end), minDiff, numValuesInSequences, maxSeqsPerLine), 'w'); % append seqs
seqs = fscanf(fID, '%d', [numel(polarAngles) Inf])';

fprintf(fIDLong, '%d ', seqs(1,:));
prevSeq = seqs(1,:);
countSeqsPerLine = 1;
% simple append or drop based on whether start of current line is
% compatible with end of previous line
for i = 2:size(seqs, 1)
    if abs(seqs(i,1) - prevSeq(end)) >= minDiff
        fprintf(fIDLong, '%d ', seqs(i,:));
        prevSeq = seqs(i,:);
        countSeqsPerLine = countSeqsPerLine + 1;
        if countSeqsPerLine == maxSeqsPerLine
            fprintf(fIDLong, '\n');
            countSeqsPerLine = 0;
        end
    end
end

fclose(fID);
fclose(fIDLong);


%% read and append
maxSeqsPerLine = 6; % mode 6
fID = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%d.txt', ...
            polarAngles(1), polarAngles(end), minDiff, numValuesInSequences), 'r'); % append seqs
fIDLong = fopen(sprintf('randSequence%dto%dAtLeastDiff%dn%dstitchn%d.txt', ...
            polarAngles(1), polarAngles(end), minDiff, numValuesInSequences, maxSeqsPerLine), 'w'); % append seqs
seqs = fscanf(fID, '%d', [numel(polarAngles) Inf])';

fprintf(fIDLong, '%d ', seqs(1,:));
prevSeq = seqs(1,:);
countSeqsPerLine = 1;
% simple append or drop based on whether start of current line is
% compatible with end of previous line
for i = 2:size(seqs, 1)
    if abs(seqs(i,1) - prevSeq(end)) >= minDiff
        fprintf(fIDLong, '%d ', seqs(i,:));
        prevSeq = seqs(i,:);
        countSeqsPerLine = countSeqsPerLine + 1;
        if countSeqsPerLine == maxSeqsPerLine
            fprintf(fIDLong, '\n');
            countSeqsPerLine = 0;
        end
    end
end

fclose(fID);
fclose(fIDLong);