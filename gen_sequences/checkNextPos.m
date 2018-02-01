function count = checkNextPos(remainingVals, possibleOpts, order, minDiff, fID, nSequences)

assert(sum(isnan(order)) == numel(remainingVals));
possibleOpts = possibleOpts(randperm(numel(possibleOpts)));

currentPos = find(isnan(order), 1, 'first');
count = 0;
for i = 1:numel(possibleOpts)
    order(currentPos) = possibleOpts(i);
    if currentPos == 2
        fprintf('2nd position try %d/%d...', i, numel(possibleOpts));
    end
    remainingValsNew = remainingVals;
    remainingValsNew(find(remainingValsNew == order(currentPos), 1, 'first')) = [];
    if currentPos < numel(order) - 1
        % unique() is slow
%         possibleOptsNext = unique(remainingValsNew(abs(remainingValsNew - order(currentPos)) >= minDiff));
        possibleOptsNext = remainingValsNew(abs(remainingValsNew - order(currentPos)) >= minDiff);
    elseif currentPos == numel(order) - 1
        % one left
        possibleOptsNext = remainingValsNew;
        assert(numel(possibleOptsNext) == 1);
        if ~(abs(remainingValsNew - order(currentPos)) >= minDiff && abs(remainingValsNew - order(1)) >= minDiff)
            return;
        end
    else
        % none left. print and return
        for j = 1:numel(order)
            fprintf(fID, '%d ', order(j));
        end
        fprintf(fID, '\n');
        count = 1;
        return;
    end
    
    added = checkNextPos(remainingValsNew, possibleOptsNext, order, minDiff, fID, nSequences);
    if currentPos == 2
        fprintf('%d added...\n', added);
    end
    count = count + added;
    if count >= nSequences
        return;
    end
end