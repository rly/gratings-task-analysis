function struct2varfn(s)
% take a structure and convert its fields into workspace variables
% source: http://stackoverflow.com/a/3470731
fn = fieldnames(s);
sc = struct2cell(s);
for i = 1:numel(fn)
    assignin('caller', fn{i}, sc{i});
end
