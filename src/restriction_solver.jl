"""
A RestrictionSolver is a solver used specifically for B splines,
taking advantage of the sparse structure.
"""
struct RestrictionSolver{ELT} <: FE_Solver{ELT}
    Aop::DictionaryOperator

    A::Matrix{ELT}
    # Mapping operator of dict to its boundary overlapping prolates
    BE::DictionaryOperator
    # Selection operator for the collocation points in the support of the boundary overlapping prolates
    GR::DictionaryOperator

    cutoff

    scratch_b::Array{ELT}
    scratch_b_linear::Vector{ELT}
    scratch_y1_native::Array{ELT}
    y1::Vector{ELT}
    function RestrictionSolver{ELT}(A::DictionaryOperator, BE::DictionaryOperator, GR::DictionaryOperator; cutoff=default_cutoff(A), options...) where {ELT}
        new(A, truncated_svd(matrix(GR*A*BE),cutoff), BE, GR, cutoff, zeros(ELT,size(dest(GR))...), zeros(ELT, length(dest(GR))), zeros(ELT, size(src(BE))...), zeros(ELT, length(src(BE))))
    end
end

src(t::RestrictionSolver) = dest(t.Aop)
dest(t::RestrictionSolver) = src(t.Aop)
inv(t::RestrictionSolver) = t.Aop

RestrictionSolver(A::DictionaryOperator, BE::DictionaryOperator, GR::DictionaryOperator; options...) =
    RestrictionSolver{eltype(A)}(A, BE, GR; options...)

function apply!(S::RestrictionSolver, x, b::Vector)
    apply!(S.GR, S.scratch_b, b)

    mul!(S.y1, S.A, S.scratch_b)
    # y1, rr = LAPACK.gelsy!(S.scratch_A,S.scratch_b, S.cutoff)
    apply!(S.BE, x, S.y1)
end

function truncated_svd(A, cutoff)
    USV = LAPACK.gesdd!('S',A)
    S = USV[2]
    maxind = findlast(S.>(maximum(S)*cutoff))
    Sinv = 1 ./ S[1:maxind]
    USV[3][1:maxind,:]'*(Sinv.*USV[1][:,1:maxind]')
end
