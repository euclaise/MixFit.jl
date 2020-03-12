const log2π = log(2*π)
const sqrt2π = √(2*π)
const sqrtπ = √π
# Euler–Mascheroni constant
const γ = 0.577215665

# Ramanujan's approximation to the Γ function:
Γ(x::Real) = sqrtπ*((x-1)/exp(1))^(x-1) * (8*(x-1)^3 + 4*(x-1)^2 + (x-1) + 1/30)^(1/6)
