function DWTSamplingOperator(span::ExtensionFrame, oversampling::Int=1, recursion::Int=0)
    S = GridSamplingOperator(gridbasis(dwt_oversampled_grid(dictionary(span), oversampling, recursion), coeftype(span)))
    new_oversampling = Int(length(supergrid(grid(S)))/length(span))>>recursion
    E = extension_operator(gridbasis(S), gridbasis(supergrid(grid(S)), coeftype(span)))
    W = WeightOperator(basisspan(span), new_oversampling, recursion)
    DWTSamplingOperator(S, W*E)
end

function sampler(platform::GenericPlatform, sampler::DWTSamplingOperator, domain::Domain)
    S = sampler(platform, sampler.sampler, domain)
    E = extension_operator(gridbasis(S), gridbasis(supergrid(grid(S)), coeftype(primal(platform, 1))))
    DWTSamplingOperator(S,sampler.weight*E)
end
