"""
Helper functions for constructing Monte Carlo degrees of freedom (dof)
arrays in a single place so that observables share the same logic.
"""

"""
    build_dof(diagpara; include_probe=false)

Return a vector of integer vectors compatible with `MCIntegration.Configuration`
based on the diagram parameters in `diagpara`.

Each entry contains the number of independent imaginary-time variables
(`tau = p.totalTauNum - 1`) and spatial variables (`spatial = p.innerLoopNum - 1`).
Set `include_probe=true` when an observable introduces an additional
probe time (e.g. doublon estimators), which prepends a `1` to every dof entry.
"""
function build_dof(diagpara; include_probe::Bool=false)
    dof = Vector{Vector{Int}}(undef, length(diagpara))
    for (i, p) in enumerate(diagpara)
        tau = p.totalTauNum - 1
        spatial = p.innerLoopNum - 1
        if include_probe
            dof[i] = [tau, spatial, 1]
        else
            dof[i] = [tau, spatial]
        end
    end
    return dof
end
