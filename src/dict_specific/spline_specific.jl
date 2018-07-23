
using CompactTranslatesDict: degree
function first_index(b::BSplineTranslatesBasis, x::Real)
    ii, on_edge = interval_index(b, x)
    d = degree(b)
    if d == 0
        return ii, 1
    end
    if on_edge
        return mod(ii-2, length(b))+1, d
    else
        return mod(ii-1, length(b))+1, d+1
    end
end

_element_spans_one(b::BSplineTranslatesBasis) = degree(b) == 0
