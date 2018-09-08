
function targ_struct = safeStructAssign(targ_struct, ref_struct, exclude_fields, cutdown)
% 
% function targ_struct = safeStructAssign(targ_struct, ref_struct, exclude_fields, cutdown)
% 
% 

if nargin<3 || isempty(exclude_fields)
    exclude_fields = {};
elseif ~iscell(exclude_fields)
    exclude_fields = {exclude_fields};
end

if nargin<4
    cutdown = true;
end

    %ensure targ_struct and ref_struct are structures and that targ_struct has the same
    %size as ref_struct
if ~isstruct(targ_struct) || ~isstruct(ref_struct)
    return;
else
    if length(targ_struct)<length(ref_struct)
        for i = length(targ_struct)+1 : length(ref_struct)
            targ_struct(i) = targ_struct(end);
        end
    elseif length(targ_struct)>length(ref_struct) && cutdown
        targ_struct = targ_struct(1:length(ref_struct));        
    end
end

    %override existing values of targ_struct with those in ref_struct
names = fieldnames(ref_struct);
for i = 1 : length(names)
    if isfield(targ_struct,names{i}) && ~ismember(names{i},exclude_fields)
        for j = 1 : length(ref_struct)
            targ_struct(j).(names{i}) = ref_struct(j).(names{i});
        end
    end
end




