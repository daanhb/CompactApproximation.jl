function BasisFunctions.DWTSamplingOperator(span::FrameFun.ExtensionFrame, oversampling::Int=1, recursion::Int=0)
    S = BasisFunctions.GridSamplingOperator(gridbasis(BasisFunctions.dwt_oversampled_grid(dictionary(span), oversampling, recursion), coeftype(span)))
    new_oversampling = Int(length(supergrid(grid(S)))/length(span))>>recursion
    E = extension_operator(gridbasis(S), gridbasis(supergrid(grid(S)), coeftype(span)))
    W = BasisFunctions.WeightOperator(FrameFun.basisspan(span), new_oversampling, recursion)
    BasisFunctions.DWTSamplingOperator(S, W*E)
end

function BasisFunctions.sampler(platform::BasisFunctions.GenericPlatform, sampler::BasisFunctions.DWTSamplingOperator, domain::Domain)
    S = BasisFunctions.sampler(platform, sampler.sampler, domain)
    E = extension_operator(gridbasis(S), gridbasis(supergrid(grid(S)), coeftype(primal(platform, 1))))
    BasisFunctions.DWTSamplingOperator(S,sampler.weight*E)
end
