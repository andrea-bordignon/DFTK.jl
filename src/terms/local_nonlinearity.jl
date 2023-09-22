"""
Local nonlinearity, with energy ∫f(ρ) where ρ is the density. The firts and second derivative of f must be provided.
"""
struct LocalNonlinearity
    f
    fρ
    fρρ
end

struct TermLocalNonlinearity{TF} <: TermNonlinear
    f::TF
    fρ::TF
    fρρ::TF
end

(L::LocalNonlinearity)(::AbstractBasis) = TermLocalNonlinearity(L.f,L.fρ,L.fρρ) # this syntax defines a function that can be called by an object of type ::LocalNonlinearity on an object of type ::AbstractBasis. It seems redundant but it is necessary for the compatibility with the rest of the codebase probabily 

function ene_ops(term::TermLocalNonlinearity, basis::PlaneWaveBasis{T}, ψ, occupation;
                 ρ, kwargs...) where {T}
    E = sum(term.fρ -> convert_dual(T, term.fρ), term.f.(ρ)) * basis.dvol
    potential = convert_dual.(T, term.fp.(ρ))

    # In the case of collinear spin, the potential is spin-dependent
    ops = [RealSpaceMultiplication(basis, kpt, potential[:, :, :, kpt.spin])
           for kpt in basis.kpoints]
    (; E, ops)
end


function compute_kernel(term::TermLocalNonlinearity, ::AbstractBasis{T}; ρ, kwargs...) where {T}
    Diagonal(vec(convert_dual.(T, term.fpp.(ρ))))
end

function apply_kernel(term::TermLocalNonlinearity, ::AbstractBasis{T}, δρ; ρ, kwargs...) where {T}
    convert_dual.(T, term.fpp.(ρ) .* δρ)
end
