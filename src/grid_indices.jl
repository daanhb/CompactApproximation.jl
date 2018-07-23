function _grid_index_limits_in_element_support(B::Dictionary, g::AbstractEquispacedGrid, i)
    dx = stepsize(g)
    x0 = g[1]
    s = support(B,i)
    if isa(s,AbstractInterval)
        start = ceil(Int,(infimum(s)-x0)/dx)
        stop = floor(Int,(supremum(s)-x0)/dx)
        !_element_spans_one(B) && ((infimum(s)-x0)/dx ≈ start) && (start += 1)
        ((supremum(s)-x0)/dx ≈ stop) && (stop -= 1)
        return (start+1, stop+1)
    else
        interval = elements(s)[1]
        # start = 0
        stop = floor(Int,(supremum(interval)-x0)/dx)
        ((supremum(interval)-x0)/dx ≈ stop) && (stop -= 1)
        # push!(indices,(start+1:stop+1)...)
        interval = elements(s)[2]
        start = ceil(Int,(infimum(interval)-x0)/dx)
        # stop = length(g)-1
        ((infimum(interval)-x0)/dx ≈ start) && (start += 1)
        return (start+1-length(g), stop+1)
    end
end

"""
Limits of the indices of `g` of points in the support of `B[i]`.
This is a tupple of (number of elements in tuple depends is equal to dimension) of two element tuples.
"""
grid_index_limits_in_element_support(B::Dictionary1d, g::AbstractGrid1d, i::Int) =
    tuple(_grid_index_limits_in_element_support(B, g, i))

grid_index_limits_in_element_support(B::TensorProductDict, g::ProductGrid, cartindex::CartesianIndex{N}) where {N} =
    [_grid_index_limits_in_element_support(s,element(g,i),cartindex[i]) for (i,s) in enumerate(elements(B))]

"""
Grid cartesian index limits of `g` of points in the support of `B[index]`.
"""
function grid_cartesian_index_limits_in_element_support(B::Dictionary, g::AbstractGrid, index)
    t = grid_index_limits_in_element_support(B, g, index)
    CartesianIndex([i[1]for i in t]...), CartesianIndex([i[2]for i in t]...)
end

"""
Grid indices of `g` of points in the support of `B[index]`.
"""
grid_index_range_in_element_support(B::Dictionary, g::AbstractGrid, index) =
    ModCartesianRange(size(g), grid_cartesian_index_limits_in_element_support(B, g, index)...)

grid_index_mask_in_element_support(B::Dictionary, g::AbstractGrid, indices) =
    grid_index_mask_in_element_support!(BitArray(size(g)), B, g, indices)

function grid_index_mask_in_element_support!(mask::BitArray, B::Dictionary, g::AbstractGrid, indices)
    fill!(mask, 0)
    for index in indices, i in grid_index_range_in_element_support(B, g, index)
        mask[i] = 1
    end
    mask
end

function grid_index_mask_in_element_support!(mask::BitArray, B::Dictionary, g::AbstractGrid, indices::BitArray)
    fill!(mask, 0)
    for i in eachindex(B)
        if indices[i]
            for j in grid_index_range_in_element_support(B, g, i)
                mask[j] = 1
            end
        end
    end
    mask
end
