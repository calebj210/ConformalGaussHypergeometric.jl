# License is MIT: http://julialang.org/license

include("transformations.jl")

@doc raw"""
    conformal2f1(a, b, c, z; rtol = 1e-14, ord = 2, esterr = false)

Compute the Gauss hypergeometric function for arbitrary parameters ``a``, ``b``, ``c``, and argument ``z`` defined by

```math
{_2}F_1(a,b;c;z) = \sum_{n = 0}^\infty \frac{(a)_n (b)_n}{(c)_n} \frac{z^n}{n!}, \quad (x)_n = \frac{\Gamma(x + n)}{\Gamma(n)}, \qquad |z| < 1
```
and by analytic continuation in the whole complex plane. Setting ``esterr = true`` also returns an estimate for the relative error in the evaluation.
end

# Examples
```jldoctest
julia> conformal2f1(1, 1/2, 1/3, 3/2)
-3.0545670655014225 + 2.0623874582040482im

julia> conformal2f1(1.1, 1.2, -1.3, .5 + 3im)
-1.0631414961412355 + 2.684842244462105im

julia> conformal2f1(1.1, 1.2, -1.3, .5 + 3im; esterr = true)
(-1.0631414961412355 + 2.684842244462105im, 2.06328233701442e-15)
```

External links:
[DLMF 15.1](https://dlmf.nist.gov/15.1)

# Implementation by 
the conformal mapping method [(paper link)](linkoncesubmitted)
"""
function conformal2f1(a, b, c, z; rtol = 1e-14, ord = 2, esterr = false)
    if isreal(z)
        z = real(z) - 0im
    end

    if abs2(z) < 1 && abs2(1 - z) < 1
        if real(z) <= 0.5
            trans = [:z, :zoverzminusone]
        else
            trans = [:oneminusz, :oneminusoneoverz]
        end
    else
        if real(z) <= 0.5
            trans = [:oneoverz, :oneoveroneminusz, :zoverzminusone, :z, :oneminusz]
        else
            trans = [:oneoverz, :oneoveroneminusz, :oneminusoneoverz, :z, :oneminusz]
        end
    end

    val = compare(a, b, c, z, trans; rtol, ord, esterr)
    return val
end

# Compare the evaluations of 2F1 using 
function compare(a, b, c, z, trans; ord = 4, esterr = false, kwargs...)
    dif = Inf
    idx = [0,0]

    N = length(trans)
    vals = zeros(ComplexF64, ord*N)
    for n ∈ 1:ord*N - 1
        ordn = ((n - 1) ÷ N) + 1
        trann = trans[((n - 1) % N) + 1]
        if iszero(vals[n])
            vals[n] = transformations[trann](a, b, c, z, ordn)
        end
        for k ∈ (n+1):ord*N
            ordk = ((k - 1) ÷ N) + 1
            trank = trans[((k - 1) % N) + 1]
            if iszero(vals[k])
                vals[k] = transformations[trank](a, b, c, z, ordk)
            end

            if isapprox(vals[n], vals[k]; kwargs...)
                if esterr
                    return ((vals[n] + vals[k]) / 2, abs(vals[n] - vals[k]) / max(abs(vals[n]), abs(vals[k])))
                else
                    return (vals[n] + vals[k]) / 2
                end
            elseif abs(vals[n] - vals[k]) < dif
                dif = abs(vals[n] - vals[k])
                idx = [n,k]
            end
        end
    end

    # @warn "Tolerance not met, answer within a relative tolerance of $(abs(dif / vals[first(idx)]))."

    if idx == [0,0]
        if esterr
            return (NaN + NaN * im, NaN)
        else
            return NaN + NaN * im
        end
    elseif esterr
        return (sum(vals[idx]) / 2, dif / max(abs.(vals[idx])...))
    else
        return sum(vals[idx]) / 2
    end
end

