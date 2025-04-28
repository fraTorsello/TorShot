% utils/getfield_safe.m
function val = getfield_safe(s, field, default)
% Restituisce il valore del campo 'field' dalla struct 's'.
% Se il campo non esiste o è vuoto, restituisce il valore 'default'.
% INPUT:
%   s: La struct da cui leggere.
%   field: Il nome del campo (stringa).
%   default: Il valore da restituire se il campo manca o è vuoto.
% OUTPUT:
%   val: Il valore del campo o il valore di default.

    if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end