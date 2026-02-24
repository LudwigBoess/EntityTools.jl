using Roots

const A = 1.187
const B = 1.251
const C = 0.714
const D = 0.936

γ_e(E_int) = (A + B*E_int) / ( C + D * E_int )

# Sridhar et al. (2021), Eq. C5
E(E_int, γ_bar, Γ) = (γ_bar - Γ) * Γ / ( 1 + γ_e(E_int) * (Γ^2 - 1) ) - E_int

"""
    get_Temp(N::T, T00::T, V1::T, V2::T, V3::T) where T <: Real

Compute the temperature based on Sridhar et al. (2021, https://arxiv.org/pdf/2107.00263), Appendix C.
It uses the number density `N`, total energy `T00`, and velocity components `V1`, `V2`, and `V3`.
The function uses a root-finding method to solve for the internal energy and then converts it to temperature using the relation T = (2/3) * E_int.
"""
function get_Temp(N::T, T00::T, V1::T, V2::T, V3::T) where T <: Real

    v_bulk = V1^2 + V2^2 + V3^2
    γ_bar = T00 / N
    Γ = 1 / sqrt(1 - v_bulk)

    _E(E_int) = E(E_int, γ_bar, Γ)
    Temp = 0.0
    try
        Temp = find_zero(_E, (-1000000.0, 1000000.0), verbose=false)
    catch
        Temp = NaN
    end

    # internal energy to temperature conversion
    Temp *= 2/3
    return Temp
end