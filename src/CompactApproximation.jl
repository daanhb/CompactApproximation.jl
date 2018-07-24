__precompile__()
module CompactApproximation

using Domains, BasisFunctions, CompactTranslatesDict, FrameFun, StaticArrays

using FrameFun: FE_Solver, default_cutoff, mask, flatten, domain
using BasisFunctions: GenericPlatform, Platform, domain

import BasisFunctions: grid_restriction_operator, restriction_operator, src, dest, apply!

if VERSION < v"0.7-"
    CartesianIndices = CartesianRange
    Nothing = Void
    mul! = A_mul_B!
else
    using LinearAlgebra
    using LinearAlgebra: LAPACK
end
include("dict_specific/dict_specific.jl")
include("subgrid.jl")

include("ModCartesianIndices.jl")
include("grid_indices.jl")
include("coefficient_indices.jl")

include("decomposition.jl")
include("azssolver_restriction_operators.jl")

include("restriction_solver.jl")
include("azssolver.jl")
include("solve.jl")


# azselection_restriction_operators(primal::Union{WaveletBasis,WaveletTensorDict}, gamma::AbstractGrid, omega::AbstractGrid, domain::Domains.Domain) =
#     azselection_restriction_operators(primal, wavelet_dual(primal), gamma, omega, domain)


export AZSSolver, decomposition_solve, RestrictionSolver, azs_solve, boundary_grid

end
