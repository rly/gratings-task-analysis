function m = getMinPValueOverLoc(statsByLoc)

m = min(arrayfun(@(x) x.p, statsByLoc));