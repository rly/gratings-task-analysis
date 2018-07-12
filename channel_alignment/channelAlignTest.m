%% setup 1D
clear;

actual = [1 3 6 8 12 14 13 9 7 4 1 3 5 6 8]';
nVals = 8;
xind1 = 8:15; % 8 values
xind2 = 3:10; % actual shift is -5
ref1 = 3;
ref2 = -2;
sigma1 = 2;
sigma2 = 2;
x1 = actual(xind1) + ref1 + randn(nVals, 1) * sigma1;
x2 = actual(xind2) + ref2 + randn(nVals, 1) * sigma2;

fprintf('x1:\n');
disp(x1)
fprintf('x2:\n');
disp(x2)

%% brute-force search
fprintf('\n');

minShift = -10;
maxShift = 10;
shiftToTest = minShift:maxShift;
nShiftToTest = numel(shiftToTest);
refDiffToTest = -10:0.5:10;
sseBest = Inf * ones(nShiftToTest, 1);
nMatchesBest = nan(nShiftToTest, 1);
for j = 1:nShiftToTest
    shift = shiftToTest(j);
    for refDiff = refDiffToTest
        sse = 0;
        nMatches = 0;
        for i = 1:nVals
            if i - shift >= 1 && i - shift <= nVals
                sse = sse + (x1(i) - x2(i - shift) - refDiff)^2;
                nMatches = nMatches + 1;
            end
        end
        if sse < sseBest(j)
            sseBest(j) = sse;
            nMatchesBest(j) = nMatches;
        end
    end
    fprintf('shift = %d, nMatches = %d, sse = %0.2f\n', shift, nMatchesBest(j), sseBest(j));
end

figure;
plot(shiftToTest, sseBest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup 2D
clear;

nActualCh = 50;
nTime = 100;
actual = 100*rand(nActualCh, nTime);
nVals = 20;
xind1 = 10:29; % 8 values
xind2 = 2:21; % actual shift is -8
ref1 = repmat(3 + randn(1, nTime), nVals, 1); % same ref shift across channels
ref2 = repmat(-2 + randn(1, nTime), nVals, 1);
actualRefDiff = ref2 - ref1;
actualRefDiff = actualRefDiff(1,:);
sigma1 = 50;
sigma2 = 50;
x1 = actual(xind1,:) + ref1 + randn(nVals, nTime) * sigma1;
x2 = actual(xind2,:) + ref2 + randn(nVals, nTime) * sigma2;

fprintf('x1:\n');
disp(x1)
fprintf('x2:\n');
disp(x2)

%% brute-force search
fprintf('\n');

minShift = -nVals+1;
maxShift = nVals-1;
shiftToTest = minShift:maxShift;
nShiftToTest = numel(shiftToTest);
nRandRefDiffsToTest = 1;%10000;
sigmaRefDiffToTest = 0;%20;
sseBest = Inf * ones(nShiftToTest, 1);
nMatchesBest = nan(nShiftToTest, 1);
refDiffBest = nan(nShiftToTest, nTime);
for j = 1:nShiftToTest
    shift = shiftToTest(j);
    for k = 1:nRandRefDiffsToTest
        refDiff = sigmaRefDiffToTest*randn(1, nTime); % 0 mean
        sse = 0;
        nMatches = 0;
        for i = 1:nVals
            if i - shift >= 1 && i - shift <= nVals
                sse = sse + sum((x1(i,:) - x2(i - shift,:) - refDiff).^2);
                nMatches = nMatches + 1;
            end
        end
        if sse < sseBest(j)
            sseBest(j) = sse;
            nMatchesBest(j) = nMatches;
            refDiffBest(j,:) = refDiff;
        end
    end
    fprintf('j = %d, shift = %d, nMatches = %d, sse = %0.2f\n', j, shift, nMatchesBest(j), sseBest(j));
end

figure;
plot(shiftToTest, sseBest);