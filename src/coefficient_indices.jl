
function _coefficient_index_limits_of_overlapping_elements(B::Dictionary1d, x::Real)
    # The init_index is the starting index of all spline elements that overlap with x
    init_index, no_elements = first_index(B,x)
    (init_index-no_elements+1, init_index)
end

"""
Limits of the indices of the coefficients of B that overlap with x.
This is a tupple of (number of elements in tuple depends is equal to dimension) of two element tuples.
"""
coefficient_index_limits_of_overlapping_elements(B::Dictionary, x::Real) =
    tuple(_coefficient_index_limits_of_overlapping_elements(B, x))

coefficient_index_limits_of_overlapping_elements(B::TensorProductDict, x::SVector{N}) where {N} =
    [_coefficient_index_limits_of_overlapping_elements(Bi,xi) for (Bi, xi) in zip(elements(B), x)]

"""
Cartesian index limits of the coefficients of B that overlap with x.
"""
function coefficient_cartesian_index_limits_of_overlapping_elementst(B::Dictionary, x)
    t = coefficient_index_limits_of_overlapping_elements(B, x)
    CartesianIndex([i[1]for i in t]...), CartesianIndex([i[2]for i in t]...)
end

"""
Range of coefficient indices of B that overlap with the point x.
"""
coefficient_index_range_of_overlapping_elements(B::Dictionary, x) =
    ModCartesianRange(size(B), coefficient_cartesian_index_limits_of_overlapping_elementst(B, x)...)

coefficient_index_mask_of_overlapping_elements(d::Dictionary, g::AbstractGrid) =
    coefficient_index_mask_of_overlapping_elements!(BitArray(size(d)), d, g)

function coefficient_index_mask_of_overlapping_elements!(mask::BitArray, B::Dictionary, g::AbstractGrid)
    fill!(mask, 0)
    for x in g, i in coefficient_index_range_of_overlapping_elements(B, x)
        mask[i] = 1
    end
    mask
end

coefficient_indices_of_overlapping_elements(dict::Dictionary, boundary::AbstractGrid) =
    find(coefficient_index_mask_of_overlapping_elements(dict, boundary))
