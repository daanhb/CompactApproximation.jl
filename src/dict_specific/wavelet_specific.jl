
function first_index(b::WaveletBasis, x::Real)
    ii, on_edge = interval_index(b, x)
    s = support(side(b), kind(b), wavelet(b))
    s1 = Int(s[1])
    L = Int(s[2]) - s1
    if L == 1
        return mod(ii-s1-1,length(b))+1, 1
    end
    if on_edge
        return mod(ii-s1-2,length(b))+1, L - 1
    else
        return mod(ii-s1-1,length(b))+1, L
    end
end

_element_spans_one(b::WaveletBasis) = support_length(side(b), kind(b), wavelet(b)) == 1
