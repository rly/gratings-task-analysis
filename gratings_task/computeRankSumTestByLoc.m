function statsByLoc = computeRankSumTestByLoc(group1, group2ByLoc, statsByLoc)
% check assumptions
% test difference in medians
% does not test difference in variance (see levene's test)

if nargin < 3
    statsByLoc = struct();
end

nLoc = numel(group2ByLoc);
for i = 1:nLoc
    if ~isempty(group2ByLoc{i})
        [statsByLoc(i).rankSum.p,~,stats] = ranksum(group1, group2ByLoc{i});
        statsByLoc(i).rankSum.zval = stats.zval;
        statsByLoc(i).rankSum.ranksum = stats.ranksum;
    else
        statsByLoc(i).rankSum.p = NaN;
        statsByLoc(i).rankSum.zval = NaN;
        statsByLoc(i).rankSum.ranksum = NaN;
    end
end