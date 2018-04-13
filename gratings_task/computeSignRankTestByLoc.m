function statsByLoc = computeSignRankTestByLoc(group1ByLoc, group2ByLoc, statsByLoc)
% check assumptions
% test difference in medians
% does not test difference in variance (see levene's test)

if nargin < 3
    statsByLoc = struct();
end

nLoc = numel(group2ByLoc);
assert(nLoc == numel(group1ByLoc));
for i = 1:nLoc
    if ~isempty(group1ByLoc{i}) && ~isempty(group2ByLoc{i})
        assert(numel(group1ByLoc{i}) == numel(group2ByLoc{i}));
        [statsByLoc(i).signRank.p,~,stats] = signrank(group1ByLoc{i} - group2ByLoc{i});
        statsByLoc(i).signRank.zval = stats.zval;
        statsByLoc(i).signRank.signedrank = stats.signedrank;
    else
        statsByLoc(i).signRank.p = NaN;
        statsByLoc(i).signRank.zval = NaN;
        statsByLoc(i).signRank.signedrank = NaN;
    end
end