SYSTEM_SIZE=8000


"""
Assign a sequence number to each element in bin.

Elements with all numbers even are to be solved first and are assigned 1.
Elements with all numbers odd are to be solved last.
"""
function assign_sequence_nro(bins)
    L = length(bins)
    M = length(bins[1])
    seq = Array{Int}(L)
    a = Array{Bool}(M)
    b = Array{Int}(M)
    for i in 1:L
        b .= bins[i]
        a .= iseven.(b)
        seq[i] = (1<<M)-_int(a)
    end
    seq
end

"Convert array of ones and zeros to integer. [1,1,0] => 110b = 6. "
function _int(a)
    i = 0
    for ai in a
        i = 2*i+ai
    end
    i
end

"""
Clasifies the coefficient indices activated in `coefficient_mask` in `depth` dimensions.
"""
function classified_indices(coefficient_mask::AbstractArray{Bool}, primal::TensorProductDict, gamma::ProductGrid, depth::Int; no_blocks::Union{Int,Nothing}=nothing)
    sprimal = size(primal)
    # The cartesian indices of the activated coefficients
    cart_indices = [CartesianIndex(ind2sub(sprimal, i) ) for i in find(coefficient_mask)]
    cart_indices_matrix = zeros(Int,length(cart_indices), depth)
    for i in 1:depth
        # Transfrom one dimensian of the cartesian indices to an array
        cart_indices_matrix[:,i] .= getindex.(cart_indices,i)
        if no_blocks==nothing
            # Determin the space between two coefficient regions in the ith dimension
            # Take into account the spacing of the collocation points and the width
            # or the primal and the dual basis.
            # (it scales, but other methods may be better)
            g1d = element(gamma,i)
            primal1d = element(primal, i)
            dual1d = wavelet_dual(primal1d)
            OS = cld(length(g1d),length(primal1d))
            no_samples_in_1d = cld(SYSTEM_SIZE,OS*no_overlapping_elements(dual1d))
            no_coeffs_in_1d_other = ceil(Int, fld(no_samples_in_1d, OS)^(1/(max(1,depth-1))))
            no_coeffs_in_1d_mid = ceil(Int, no_overlapping_elements(primal1d)^(1/(max(1,depth-1))))

            # Divide all coefficients in a single dimension in an even n.o. partitions (m).
            m = cld(length(primal1d),max(no_coeffs_in_1d_mid, no_coeffs_in_1d_other))
            m = isodd(m) ? m+1 : m
            # +1 to be on the save side.
            m = (length(primal1d)+1)/m
        else
            isodd(no_blocks) && no_blocks!=1 && (warn("An odd number of blocks may lead to errors"))
            primal1d = element(primal, i)
            m = (length(primal1d)+1)/no_blocks
        end
        # Reuse the array to minimize allocation
        cart_indices_matrix[:,i] .= Int.(cld.(cart_indices_matrix[:,i], m))
    end
    # Create a vector instead of a matrix.
    # Then reducing dimensions is not necessary in methods using this output.
    cart_indices, [tuple(cart_indices_matrix[i,:]...) for i in 1:size(cart_indices,1)]
end
