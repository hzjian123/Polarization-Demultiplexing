function s1 = mergestruct(s1, s2)

if isa(s1, 'struct') && isa(s2, 'Struct')
    s2 = s2.toStruct;
end

if isa(s1, 'Struct') && isa(s2, 'struct')
    s2 = Struct(s2);
end

isstruct = @(x) isa(x, 'struct') || isa(x, 'Struct');

if ~isstruct(s1)
    error('arg#1 must be struct')
end

if ~isstruct(s2)
    error('arg#2 must be struct')
end

f1 = fieldnames(s1);
f2 = fieldnames(s2);

for i = 1:numel(f2)
    f = f2{i};
    if ismember(f, f1)
        if isstruct(s1.(f)) && isstruct(s2.(f))
            s1.(f) = mergestruct(s1.(f), s2.(f));
        end
    else
        s1.(f) = s2.(f);
    end
end
