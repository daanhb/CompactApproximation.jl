
"""
Number of basis elements overlapping with a point.
"""
no_overlapping_elements(dict::Dictionary) = ceil(Int,support_length_of_compact_function(dict)*length(dict))

"""
Makes sure that (i-1)/N <= x < i/N holds.
Return (i, true) if x ≈ i/N
else (i, false)
"""
function interval_index(B::Dictionary,x::Real)
    L = length(B)
    s = x*L
    r =  round(s)
    floor(Int,s)+1, s≈r
end

function restriction_operator(dict::Dictionary, mask::BitArray)
     indices = (VERSION < v"0.7-") ? find(mask) : LinearIndices(mask)[findall(mask)]
     IndexRestrictionOperator(dict, dict[indices], indices)
end


# include("wavelet_specific.jl")
include("spline_specific.jl")
