function findTimeSeriesMotifs(data, numReferences, subsequenceLength, numDeadSamplesAroundSubsequence)
    % adapted from cpp code for MK algorithm by Abdullah Al Mueen
    % doesn't seem to capture the lowest distance motifs... may be
    % something wrong in my adaptation    
    
    smallestDistanceSoFar = Inf;
    N = numel(data);
    motifPair = nan(1, 2);
    
    fprintf('\n');
    fprintf('Length of time series: %d\n', N);
    fprintf('Length of subsequences: %d\n', subsequenceLength);
    fprintf('Number of references: %d\n\n', numReferences);
    tic;
    
    % break data into subsequences
    numSubsequences = N - subsequenceLength + 1;
    subsequences = nan(numSubsequences, subsequenceLength);
    for i = 1:numSubsequences
        subsequences(i,:) = data(i:i+subsequenceLength-1);
    end
    
    % normalize data
    normSubsequences = (subsequences - mean(subsequences, 2)) ./ std(subsequences, 0, 2);
    
    % generate the reference time series randomly
    refs = nan(numReferences, subsequenceLength);
    distancesToRefs = nan(numReferences, numSubsequences);
    
    randPermutation = randperm(numSubsequences);
    randPicks = randPermutation(1:numReferences); % ensure unique references
    for r = 1:numReferences
        refs(r,:) = normSubsequences(randPicks(r),:);
        
        for i = 1:numSubsequences
            if i == randPicks(r)
                distancesToRefs(r,i) = NaN;
                continue;
            end
            distancesToRefs(r,i) = computeEuclideanDistanceEarlyStop(refs(r,:), normSubsequences(i,:));
            
            % if the subsequence is too close, skip this pair
            if abs(i - randPicks(r)) <= numDeadSamplesAroundSubsequence
                continue;
            end
            
            % update smallestDistanceSoFar
            if distancesToRefs(r,i) < smallestDistanceSoFar
                smallestDistanceSoFar = distancesToRefs(r,i);
                motifPair = [i randPicks(r)];
                fprintf('New smallest distance is %0.2f and (%d, %d) is the new motif pair\n', smallestDistanceSoFar, motifPair);
            end
        end
    end
    fprintf('References have been picked and distances have been computed.\n');

    % sort references by SD, largest to smallest
    sdDistancesToRefs = nanstd(distancesToRefs, 0, 2);
    [~,sdDistancesToRefSortOrder] = sort(sdDistancesToRefs, 1, 'descend'); % could just use max? but useful for early abort 
    
    % sort subsequences for reference with largest SD by distance to
    % reference, smallest to largest
    [~,distancesToRefSortOrder] = sort(distancesToRefs(sdDistancesToRefSortOrder(1),:));
    
    fprintf('Orderings have been computed and search has begun!\n\n');

    % compute the distances between a pair of time series only when there
    % is a potential motif

    offset = 0;
    abandon = false;
    while (~abandon && offset < numSubsequences)
        abandon = true;
        offset = offset + 1;

        % for each subsequence, pick the subsequence with smallest distance
        % to ref and one with a larger distance (offset)
        for i = 1:numSubsequences-offset
            seq1 = distancesToRefSortOrder(i);
            seq2 = distancesToRefSortOrder(i + offset);
            if seq1 > seq2 % aesthetics -- swap to make sure seq < seq2
                temp = seq1;
                seq1 = seq2;
                seq2 = temp;
                clear temp;
            end
            
            % if the subsequences are too close, skip this pair
            if abs(seq1 - seq2) <= numDeadSamplesAroundSubsequence
                continue;
            end
            
            % compute the difference between the distance ofright to ref
            % and the distance of left to ref. if all differences are
            % smaller than the smallest distance so far, update the
            % smallest distance so far, and re-run the search with a larger
            % offset. if any of the differences are larger than the
            % smallest distance so far and this is the case for all
            % subsequences, then abandon the search
            
            % TODO code early abort -- matters if many references
            diffDistancesToRef = abs(distancesToRefs(:,seq1) - distancesToRefs(:,seq2));
            diffDistancesToRef(isnan(diffDistancesToRef)) = [];

            if all(diffDistancesToRef < smallestDistanceSoFar)
                abandon = false;
                d = computeEuclideanDistanceEarlyStop(normSubsequences(seq1,:), normSubsequences(seq2,:), smallestDistanceSoFar);
                if d < smallestDistanceSoFar
                    smallestDistanceSoFar = d;
                    motifPair = [seq1 seq2];
                    fprintf('New smallest distance is %0.2f and (%d, %d) is the new motif pair\n', smallestDistanceSoFar, motifPair);
                end
            end
        end
    end
    fprintf('\nFinal Motif is the pair (%d, %d) and the Motif distance is %0.2f\n', motifPair, smallestDistanceSoFar);
%     subsequences(motifPair,:)
%     computeEuclideanDistanceEarlyStop(normSubsequences(motifPair(1),:), normSubsequences(motifPair(2),:))
    fprintf('Execution Time was: %0.2f seconds\n', toc);
end

% Calculates the Euclidean distance between two time series x and y. If the distance is
% larger than the best so far, it stops computing and returns the approximate
% distance. To get exact distance the stopDistance argument should be omitted.
% **not sure if early stop is faster than computing distance using vectors -- 
% depends on length of x**
function distance = computeEuclideanDistanceEarlyStop(x, y, stopDistance)
    if nargin < 3
        stopDistance = Inf;
    end
    N = numel(x);
    assert(N == numel(y));
    
    sumSquares = 0;
    stopDistanceSquared = stopDistance^2;
    for i = 1:N
        sumSquares = sumSquares + (x(i) - y(i)).^2;
        if sumSquares >= stopDistanceSquared
            break;
        end
    end
    distance = sqrt(sumSquares);
end

