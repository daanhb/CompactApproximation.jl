
if VERSION < v"0.7-"
    struct ModCartesianIndices{N}
        size::NTuple{N,Int}
        range::CartesianIndices{CartesianIndex{N}}
    end
    ModCartesianIndices(size::NTuple{N,Int}, index1::CartesianIndex{N}, index2::CartesianIndex{N}) where {N} =
        ModCartesianIndices(size, CartesianIndices(index1, index2))

    Base.length(m::ModCartesianIndices) = length(m.range)
    Base.start(m::ModCartesianIndices{N}) where {N} = start(m.range)
    @generated function Base.next(m::ModCartesianIndices{N}, state) where N
        t = Expr(:tuple, [:(mod(index[$i]-1,m.size[$i])+1) for i in 1:N]...)
        return quote
            index, state = next(m.range, state)
            CartesianIndex($t), state
        end
    end

    @generated function Base.next(m::ModCartesianIndices{1}, state)
        t = :(mod(index[1]-1,m.size[1])+1)
        return quote
            index, state = next(m.range, state)
            $t, state
        end
    end

    Base.done(m::ModCartesianIndices{N}, state) where N= done(m.range,  state)

else
    struct ModCartesianIndices{N}
        size::NTuple{N,Int}
        range::CartesianIndices{N}
    end
    # ModCartesianIndices(size::NTuple{N,Int}, index1::CartesianIndex{N}, index2::CartesianIndex{N}) where {N} =
    #     ModCartesianIndices(size, CartesianIndices(index1, index2))
    ModCartesianIndices(size::NTuple{N,Int}, index_ranges::AbstractUnitRange{Int64}...) where {N} =
        ModCartesianIndices(size, CartesianIndices(index_ranges))


    Base.length(m::ModCartesianIndices) = length(m.range)

    @generated function Base.iterate(m::ModCartesianIndices{1})
        t = :(mod(index[1]-1,m.size[1])+1)
        return quote
            index, state = iterate(m.range)
            $t, state
        end
    end

    @generated function Base.iterate(m::ModCartesianIndices{1}, state)
        t = :(mod(index[1]-1,m.size[1])+1)
        return quote
            tuple = iterate(m.range, state)
            if tuple!=nothing
                index, state = tuple
                $t, state
            end
        end
    end

    @generated function Base.iterate(m::ModCartesianIndices{N}) where N
        t = Expr(:tuple, [:(mod(index[$i]-1,m.size[$i])+1) for i in 1:N]...)
        return quote
            index, state = iterate(m.range)
            CartesianIndex($t), state
        end
    end

    @generated function Base.iterate(m::ModCartesianIndices{N}, state) where N
        t = Expr(:tuple, [:(mod(index[$i]-1,m.size[$i])+1) for i in 1:N]...)
        return quote
            tuple = iterate(m.range, state)
            if tuple != nothing
                index, state = tuple
                CartesianIndex($t), state
            end
        end
    end

    Base.eltype(::Type{ModCartesianIndices{N}}) where {N} = CartesianIndex{N}

    Base.eltype(::Type{ModCartesianIndices{1}}) = Int64
end
