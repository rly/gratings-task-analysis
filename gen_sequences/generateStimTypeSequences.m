% generate random shufflings of stimuli types where the polar angle
% parameter is ordered according to the polar angle sequence file input so
% that consecutive stimuli are placed at least pi/2 from each other. final
% shufflings are output to the files stimTypes_mode1_nSeq12_run{runID}.txt.

%% set up conditions 
mode = 1;
% polarAngleSequencesFileName = 'randSequence-8to8AtLeastDiff7n17stitchn4.txt';
polarAngleSequencesFileName = 'randSequence-11to11AtLeastDiff8n19stitchn3.txt';

eccs = [180 240 320]; % should match options in file
polarAngleOpts = [-11 -8:8 11]; % should match options in file
polarAngleMultFactor = pi/14;
gratingAngles = 45;
contrasts = 0.5;
driftDirections = 0;

numRuns = 5;
numSequencesPerRun = 12;
for runIndex = 1:numRuns
    startSeqIndex = (runIndex-1)*numSequencesPerRun + 1;
    endSeqIndex = runIndex * numSequencesPerRun;
    for seqIndex = startSeqIndex:endSeqIndex
        stimTypeSequence = generateStimTypes(eccs, polarAngleOpts, ...
                polarAngleSequencesFileName, polarAngleMultFactor, ...
                seqIndex, gratingAngles, contrasts, driftDirections, mode);
        numStimTypes = size(stimTypeSequence, 1);

        fileName = sprintf('stimTypes_mode%d_seq%d.txt', mode, seqIndex);
        fprintf('Writing to file: %s\n', fileName);
        fID = fopen(fileName, 'w');
        for j = 1:numStimTypes
            fprintf(fID, '%d\t%f\t%f\t%f\t%f\t%f\t%f\n', stimTypeSequence(j,:));
        end
        fclose(fID);
    end

    %% append all the sequences together
    % the polar angle differences between sequences should be OK (see
    % generateRandSequence.m)
    allSeqsFileName = sprintf('stimTypes_mode%d_nSeq%d_run%d.txt', mode, numSequencesPerRun, runIndex);
    fIDAllSeqs = fopen(allSeqsFileName, 'w');
    fprintf('Combining in file: %s\n', allSeqsFileName);
    for seqIndex = startSeqIndex:endSeqIndex
        fileName = sprintf('stimTypes_mode%d_seq%d.txt', mode, seqIndex);
        fID = fopen(fileName, 'r');
        fwrite(fIDAllSeqs, fread(fID));
        fclose(fID);
    end
    fclose(fIDAllSeqs);
end


%% set up conditions 
mode = 6;
% polarAngleSequencesFileName = 'randSequence-8to8AtLeastDiff7n17stitchn4.txt';
polarAngleSequencesFileName = 'randSequence-11to11AtLeastDiff8n19stitchn6.txt';

eccs = [90 110 140 180 240 320]; % should match options in file
polarAngleOpts = [-11 -8:8 11]; % should match options in file
polarAngleMultFactor = pi/14;
gratingAngles = 45;
contrasts = 0.5;
driftDirections = 0;

numRuns = 5;
numSequencesPerRun = 6;
for runIndex = 1:numRuns
    startSeqIndex = (runIndex-1)*numSequencesPerRun + 1;
    endSeqIndex = runIndex * numSequencesPerRun;
    for seqIndex = startSeqIndex:endSeqIndex
        stimTypeSequence = generateStimTypes(eccs, polarAngleOpts, ...
                polarAngleSequencesFileName, polarAngleMultFactor, ...
                seqIndex, gratingAngles, contrasts, driftDirections, mode);
        numStimTypes = size(stimTypeSequence, 1);

        fileName = sprintf('stimTypes_mode%d_seq%d.txt', mode, seqIndex);
        fprintf('Writing to file: %s\n', fileName);
        fID = fopen(fileName, 'w');
        for j = 1:numStimTypes
            fprintf(fID, '%d\t%f\t%f\t%f\t%f\t%f\t%f\n', stimTypeSequence(j,:));
        end
        fclose(fID);
    end

    %% append all the sequences together
    % the polar angle differences between sequences should be OK (see
    % generateRandSequence.m)
    allSeqsFileName = sprintf('stimTypes_mode%d_nSeq%d_run%d.txt', mode, numSequencesPerRun, runIndex);
    fIDAllSeqs = fopen(allSeqsFileName, 'w');
    fprintf('Combining in file: %s\n', allSeqsFileName);
    for seqIndex = startSeqIndex:endSeqIndex
        fileName = sprintf('stimTypes_mode%d_seq%d.txt', mode, seqIndex);
        fID = fopen(fileName, 'r');
        fwrite(fIDAllSeqs, fread(fID));
        fclose(fID);
    end
    fclose(fIDAllSeqs);
end
