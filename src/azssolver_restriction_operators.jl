

function azselection_restriction_operators(fplatform::GenericPlatform, i; options...)
    platform = fplatform.super_platform
    s = sampler(fplatform, i)
    omega = grid(s)
    gamma = supergrid(omega)
    d = domain(primal(fplatform, i))
    azselection_restriction_operators(primal(platform, i), gamma, omega, d)
end

azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, domain::Domains.Domain) =
    azselection_restriction_operators(primal, primal, gamma, omega, domain)


"""
Grid restriction and dictionary restriction operators to restrict the system when approximating with compact frame elements.

`primal` and `dual` should be the not restricted basis, `gamma` the (oversampled) grid of `primal`, `omega` is `gamma` restricted to `domain`
"""
function azselection_restriction_operators(primal::Dictionary, dual::Dictionary, gamma::AbstractGrid, omega::Union{MaskedGrid,IndexSubGrid}, domain::Domains.Domain)
    bound = boundary_grid(gamma, domain)
    coefficient_mask = coefficient_index_mask_of_overlapping_elements(dual, bound)
    _azselection_restriction_operators(primal, gamma, omega, coefficient_mask)
end

function _azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, coefficient_mask)
    grid_mask = grid_index_mask_in_element_support(primal, gamma, coefficient_mask)
    grid_mask .= grid_mask .& mask(omega)
    __azselection_restriction_operators(primal, gamma, omega, grid_mask)
end

function __azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, grid_mask)
    DMZ = gamma[grid_mask]
    system_coefficient_mask = coefficient_index_mask_of_overlapping_elements(primal, DMZ)
    gr = restriction_operator(gridbasis(omega, coeftype(primal)), gridbasis(DMZ, coeftype(primal)))
    dr = restriction_operator(primal, system_coefficient_mask)
    dr, gr
end

"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
boundary_support_grid(basis, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid}) =
    boundary_support_grid(basis, basis, boundary_grid, omega_grid)

function boundary_support_grid(basis::Dictionary, dual::Dictionary, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid})
    boundary_element_m = coefficient_index_mask_of_overlapping_elements(dual, boundary_grid)
    gamma = supergrid(omega_grid)
    m = grid_index_mask_in_element_support(basis, gamma, boundary_element_m)
    m .= m .& mask(omega_grid)
    MaskedGrid(gamma,m)
end

function grid_restriction_operator(src::Dictionary, dest::Dictionary, src_grid::Union{IndexSubGrid,MaskedGrid}, dest_grid::MaskedGrid)
    @assert supergrid(src_grid) == supergrid(dest_grid)
    IndexRestrictionOperator(src, dest, relative_indices(dest_grid,src_grid))
end


# spline_util_restriction_operators(platform::GenericPlatform, i) =
#     spline_util_restriction_operators(primal(platform, i), sampler(platform, i))
#
# spline_util_restriction_operators(dict::Dictionary, sampler::GridSamplingOperator) =
#     spline_util_restriction_operators(dict, grid(sampler))
#
# spline_util_restriction_operators(frame::ExtensionFrame, grid::Union{IndexSubGrid,MaskedGrid}) =
#     spline_util_restriction_operators(superdict(frame), grid, supergrid(grid), Domains.domain(frame))
#
# spline_util_restriction_operators(dict::TensorProductDict{N,NTuple{N1,DT},S,T}, grid::AbstractGrid) where {N,N1,DT<:ExtensionFrame,S,T} =
#     spline_util_restriction_operators(flatten(dict), grid)
#
# spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, grid::AbstractGrid, domain::Domain) =
#     spline_util_restriction_operators(dict, omega, boundary_grid(grid, domain))
#
# """
# Frame restriction operator and grid restriction operator.
# The former restricts `dict` to the elements that overlap with the boundary and
# the latter restricts `omega` to the points in the span of the dict elements that
# overlap with the boundary.
# """
# spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, boundary::AbstractGrid) =
#     _spline_util_restriction_operators(dict, omega, boundary_support_grid(dict, boundary, omega))
#
# function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, DMZ::AbstractGrid)
#     boundary_indices = coefficient_indices_of_overlapping_elements(dict, DMZ)
#     frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
#     grid_restriction = restriction_operator(gridbasis(grid), gridbasis(DMZ))
#     frame_restriction, grid_restriction
# end
#
# """
# Frame restriction operator and grid restriction operator.
#
# The former restricts `dict` to the elements that overlap with the boundary and
# the latter restricts `grid` to the points in the span of the dict elements that
# overlap with a region defined as the span of the dict elements that overlap with the boundary.
# """
# function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, boundary::AbstractGrid, DMZ::AbstractGrid)
#     boundary_indices = coefficient_indices_of_overlapping_elements(dict, boundary)
#     frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
#     grid_restriction = restriction_operator(gridbasis(grid), gridbasis(DMZ))
#     frame_restriction, grid_restriction
# end
