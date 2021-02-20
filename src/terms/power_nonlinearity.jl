"""
Power nonlinearity, with energy C ∫ρ^α where ρ is the density
"""
struct PowerNonlinearity
    C::Real
    α::Real
end
(P::PowerNonlinearity)(basis) = TermPowerNonlinearity(basis, P.C, P.α)

struct TermPowerNonlinearity <: Term
    basis::PlaneWaveBasis
    C::Real
    α::Real
end

function ene_ops(term::TermPowerNonlinearity, ψ, occ; ρ, kwargs...)
    basis = term.basis

    E = term.C * sum(ρ .^ term.α) * term.basis.dvol
    potential = @. term.C * term.α * ρ^(term.α-1)

    # In the case of collinear spin, the potential is spin-dependent
    ops = [RealSpaceMultiplication(basis, kpoint, potential[:, :, :, kpoint.spin])
           for kpoint in basis.kpoints]
    (E=E, ops=ops)
end


_pnl_kernel(C, α, ρ) = @. C * α * (α-1) * ρ^(α-2)

function compute_kernel(term::TermPowerNonlinearity;
                        ρ, kwargs...)
    K = Diagonal(vec(_pnl_kernel(term.C, term.α, ρ)))
end

function apply_kernel(term::TermPowerNonlinearity, dρ;
                      ρ, kwargs...)
    kernel = _pnl_kernel(term.C, term.α, ρ) .* dρ
end
