# MixFit.jl

[![License: CC0-1.0](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](http://joshuapritsker.com/MixFit.jl/stable)

MixFit.jl is a Julia package for fitting finite mixture models, and for nonparametric density estimation via finite mixture models ([Wang & Chee, 2012](https://doi.org/10.1177/1471082X1001200104)).  Models are fitted using random-swap expectation maximization ([Zhao, Hautamäki, Kärkkäinen, & Fränti, 2012](https://doi.org/10.1016/j.patrec.2012.06.017)), which avoids the problem of local optimums by randomly replacing components and using the final maximum.  By default, MixFit.jl supports normal, gamma, lognormal, and gumbel mixtures, but it is easy to use your own density function instead.  If you use this, it *would be nice* if you cited me, but I'm putting it under public domain so you can do whatever you want with it.

## Example

```
julia> using MixFit

julia> x = [randn(500); (randn(500).+3).*2]
    1000-element Array{Float64,1}:
     -0.3872542393168341
      0.44218146773613404
     -1.140006685489404  
     -0.11093365262262657
      0.917287805330094  
      1.3997276699755827
      ⋮                  
      3.860428712486353  
      4.924743765884938  
      4.767237225983121  
      3.6111782228201013
      5.498455397469906  

julia> densfit(x, criterion = AIC3)
    ----------------------------------

    Log-likelihood: -2408.562390297034
    AIC: 4827.124780594068
    AIC3: 4832.124780594068
    BIC: 4851.6635569889795

    Component 1:
        α: 0.5093992
        μ: 0.07551261
        σ: 1.0396444

    Component 2:
        α: 0.49060085
        μ: 6.004384
        σ: 2.0322912

    ----------------------------------
    MixModel(Float32[0.5093992, 0.49060085], Float32[0.07551261, 6.004384], Float32[1.0396444, 2.0322912], MixFit.dnorm)

julia>
```
