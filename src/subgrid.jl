
struct NBIndexList{N}
    index::NTuple{N,Int}
    size::NTuple{N,Int}
end
NBIndexList(index::Base.CartesianIndex{N}, size) where {N}  = NBIndexList(index.I, size)
NBIndexList(index::Int, size)  = NBIndexList((index,), size)

if VERSION < v"0.7-"
    @generated function Base.start(l::NBIndexList{N}) where {N}
    	startargs = fill(-1, N)
        stopargs = fill(1, N)
    	:(CartesianRange(CartesianIndex{$N}($(startargs...)), CartesianIndex{$N}($(stopargs...))), CartesianIndex{$N}($(startargs...)))
    end
    @generated function Base.next(l::NBIndexList{N}, state) where {N}
        t = Expr(:tuple, [:(if 1<=idx[$i]+l.index[$i]<=l.size[$i];idx[$i]+l.index[$i];elseif idx[$i]+l.index[$i]==0;l.size[$i];else; 1 ;end) for i in 1:N]...)
        return quote
            iter = state[1]
            iter_state = state[2]
            idx, iter_next_state = next(iter, iter_state)
            (sum(abs.(idx.I)) == 0) && ((idx, iter_next_state) = next(iter, iter_next_state))
            CartesianIndex($t),(iter, iter_next_state)
        end
    end

    @generated function Base.next(l::NBIndexList{1}, state)
        t = :(if 1<=idx[1]+l.index[1]<=l.size[1];idx[1]+l.index[1];elseif idx[1]+l.index[1]==0;l.size[1];else; 1 ;end)
        return quote
            iter = state[1]
            iter_state = state[2]
            idx, iter_next_state = next(iter, iter_state)
            (sum(abs.(idx.I)) == 0) && ((idx, iter_next_state) = next(iter, iter_next_state))
            $t,(iter, iter_next_state)
        end
    end
    function Base.done(g::NBIndexList{N}, state) where {N}
        iter = state[1]
        iter_state = state[2]
        done(iter, iter_state)
    end
else
    function Base.iterate(l::NBIndexList{N}) where {N}
        iter = CartesianIndices(ntuple(k->-1:1, Val(N)))
        CartesianIndex(ntuple(i->(if 1<=-1+l.index[i]<=l.size[i];-1+l.index[i];elseif -1+l.index[i]==0;l.size[i];else; 1 ;end),Val(N))), (iterate(iter)[2], iter)
    end

    function Base.iterate(l::NBIndexList{N}, tuple) where {N}
        item, iter = tuple
        next_tuple = iterate(iter, item)
        if next_tuple != nothing
            idx, next_state = next_tuple
            CartesianIndex(ntuple(i->if 1<=idx[i]+l.index[i]<=l.size[i];idx[i]+l.index[i];elseif idx[i]+l.index[i]==0;l.size[i];else; 1 ;end,Val(N))), (next_state, iter)
        end
    end

    function Base.iterate(l::NBIndexList{1})
        iter = CartesianIndices(ntuple(k->-1:1, Val(1)))
        (if 1<=-1+l.index[1]<=l.size[1];-1+l.index[1];elseif -1+l.index[1]==0;l.size[1];else; 1 ;end), (iterate(iter)[2], iter)
    end

    function Base.iterate(l::NBIndexList{1}, tuple)
        item, iter = tuple
        next_tuple = iterate(iter, item)
        if next_tuple != nothing
            idx, next_state = next_tuple
            (if 1<=idx[1]+l.index[1]<=l.size[1];idx[1]+l.index[1];elseif idx[1]+l.index[1]==0;l.size[1];else; 1 ;end), (next_state, iter)
        end
    end
end


"""
A Masked grid that contains the elements of grid that are on the boundary of the domain
"""
function boundary_grid(grid::AbstractGrid, domain::Domains.Domain)
    mask = boundary_mask(grid, domain);
    MaskedGrid(grid,mask);
end

boundary_grid(grid::MaskedGrid, domain::Domains.Domain) = boundary_grid(supergrid(grid), domain)


function boundary_mask(grid::AbstractGrid, domain::Domains.Domain)
    S = size(grid)
    m = (VERSION < v"0.7-") ? BitArray(S...) : BitArray(undef,S)
    m[:] .= 0
    t = true
    for i in eachindex(grid)
        if Domains.indomain(grid[i], domain)
            t = true
            for bi in NBIndexList(i, S)
                if !Domains.indomain(grid[bi], domain)
                    t = false
                    break
                end
            end
            m[i] = !t
        end
    end
    m
end


# It is assumed that all points of `from` are in `relativeto` and that the supergrids of both grids are equal.
function relative_indices(from::MaskedGrid, relativeto::Union{IndexSubGrid,MaskedGrid})
    # @assert (supergrid(from)) == (supergrid(relativeto))
    @assert length(supergrid(from)) == length(supergrid(relativeto))
    support_index = (VERSION < v"0.7-") ? Array{Int}(length(from)) : Array{Int}(undef, length(from))
    index = 1
    for (i_i,i) in enumerate(subindices(relativeto))
        if is_subindex(i, from)
            support_index[index] = i_i
            index += 1
        end
    end
    support_index
end
