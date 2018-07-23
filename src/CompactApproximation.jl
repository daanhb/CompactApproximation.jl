module CompactApproximation

using Domains, BasisFunctions, CompactTranslatesDict, FrameFun, StaticArrays

using FrameFun: FE_Solver

if VERSION < v"0.7-"
    CartesianIndices = CartesianRange
    Nothing = Void
else
    nothing
end
include("dict_specific/dict_specific.jl")


include("ModCartesianIndices.jl")
include("grid_indices.jl")
include("coefficient_indices.jl")

include("decomposition.jl")
include("azssolver_restriction_operators.jl")

include("azssolver.jl")
include("solve.jl")


# azselection_restriction_operators(primal::Union{BasisFunctions.WaveletBasis,BasisFunctions.WaveletTensorDict}, gamma::AbstractGrid, omega::AbstractGrid, domain::Domains.Domain) =
#     azselection_restriction_operators(primal, BasisFunctions.wavelet_dual(primal), gamma, omega, domain)


export AZSSolver, decomposition_solve

end
