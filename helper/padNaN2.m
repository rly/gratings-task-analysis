function out = padNaN2(data)
% Creates a padded data matrix from input structural array of spike times
% pads with NaN
% Usage: data=padNaN(data)
% Input:
% data : structural array of spike times
% Output:
% data : data matrix (NaN padded)
NC = length(data);
Nsp = cellfun(@length, {data.times});
out(1:max(Nsp),1:NC)=NaN;
for c = 1:NC;
    out(1:Nsp(c),c) = data(c).times;
end;
