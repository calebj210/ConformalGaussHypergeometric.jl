# ConformalGaussHypergeometric.jl

[![Build Status](https://github.com/calebj210/ConformalGaussHypergeometric.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/calebj210/ConformalGaussHypergeometric.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package implements the Gauss hypergeometric function `2F1(a, b; c; z)` using the conformal mapping method.

## Usage

### Evaluate 2F1

```julia
julia> conformal2f1(1.1, 1.22, 1.333, .5+.5im)
1.0013565978548216 + 1.0071914083362823im
```

### Estimate relative error 2F1 evaluations

```julia
julia> conformal2f1(1.1, 1.22, 1.333,.5+.5im, esterr=true)
(1.0013565978548216 + 1.0071914083362823im, 2.2109869152586823e-16)
```

## References
[1] Jacobs, C., Piret, C. & Fornberg, B. Efficient and accurate evaluation of the Gauss hypergeometric function. Numer Algor (2026). <https://doi.org/10.1007/s11075-026-02436-2>
