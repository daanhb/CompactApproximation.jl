struct ModCartesianRange{N}
    size::NTuple{N,Int}
    range::CartesianRange{CartesianIndex{N}}
end
ModCartesianRange(size::NTuple{N,Int}, index1::CartesianIndex{N}, index2::CartesianIndex{N}) where {N} =
    ModCartesianRange(size, CartesianRange(index1, index2))

Base.length(m::ModCartesianRange) = length(m.range)
Base.start(m::ModCartesianRange{N}) where {N} = start(m.range)
@generated function Base.next(m::ModCartesianRange{N}, state) where N
    t = Expr(:tuple, [:(mod(index[$i]-1,m.size[$i])+1) for i in 1:N]...)
    return quote
        index, state = next(m.range, state)
        CartesianIndex($t), state
    end
end

@generated function Base.next(m::ModCartesianRange{1}, state)
    t = :(mod(index[1]-1,m.size[1])+1)
    return quote
        index, state = next(m.range, state)
        $t, state
    end
end

Base.done(m::ModCartesianRange{N}, state) where N= done(m.range,  state)
