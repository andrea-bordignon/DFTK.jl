using SparseArrays

"""
Compute the index mapping between the global grids of two bases.
Returns an iterator of 8 pairs `(block_in, block_out)`. Iterated over these pairs
`x_out_fourier[block_out, :] = x_in_fourier[block_in, :]` does the transfer from
the Fourier coefficients `x_in_fourier` (defined on `basis_in`) to
`x_out_fourier` (defined on `basis_out`, equally provided as Fourier coefficients).
"""
function transfer_mapping(basis_in::PlaneWaveBasis, basis_out::PlaneWaveBasis)
    @assert basis_in.model.lattice == basis_out.model.lattice

    # TODO This logic feels rather convoluted ... think if there are ways to simplify
    idcs = map(basis_in.fft_size, basis_out.fft_size) do fft_in, fft_out
        if fft_in <= fft_out
            a = cld(fft_in, 2)
            b = fld(fft_in, 2)
            return (; in=(1:a,         (a+1):fft_in ),
                     out=(1:a, (fft_out+1-b):fft_out))
        else
            a = cld(fft_out, 2)
            b = fld(fft_out, 2)
            return (; in=(1:a, (fft_in+1-b):fft_in ),
                     out=(1:a,        (a+1):fft_out))
        end
    end

    idcs_in  = CartesianIndices.(Iterators.product(idcs[1].in,  idcs[2].in,  idcs[3].in))
    idcs_out = CartesianIndices.(Iterators.product(idcs[1].out, idcs[2].out, idcs[3].out))
    zip(idcs_in, idcs_out)
end

"""
Compute the index mapping between two bases. Returns two arrays
`idcs_in` and `idcs_out` such that `ψkout[idcs_out] = ψkin[idcs_in]` does
the transfer from `ψkin` (defined on `basis_in` and `kpt_in`) to `ψkout`
(defined on `basis_out` and `kpt_out`).
"""
function transfer_mapping(basis_in::PlaneWaveBasis{T},  kpt_in::Kpoint,
                          basis_out::PlaneWaveBasis{T}, kpt_out::Kpoint) where {T}
    idcs_in  = 1:length(G_vectors(basis_in, kpt_in))  # All entries from idcs_in
    kpt_in == kpt_out && return idcs_in, idcs_in

    # Get indices of the G vectors of the old basis inside the new basis.
    idcs_out = index_G_vectors.(Ref(basis_out), G_vectors(basis_in, kpt_in))

    # In the case where G_vectors(basis_in.kpoints[ik]) are bigger than vectors
    # in the fft_size box of basis_out, we need to filter out the "nothings" to
    # make sure that the index linearization works. It is not an issue to
    # filter these vectors as this can only happen if Ecut_in > Ecut_out.
    if any(isnothing, idcs_out)
        idcs_in  = idcs_in[idcs_out .!= nothing]
        idcs_out = idcs_out[idcs_out .!= nothing]
    end
    idcs_out = getindex.(Ref(LinearIndices(basis_out.fft_size)), idcs_out)

    # Map to the indices of the corresponding G-vectors in
    # G_vectors(basis_out, kpt_out) this array might contains some nothings if
    # basis_out has less G_vectors than basis_in at this k-point
    idcs_out = indexin(idcs_out, kpt_out.mapping)
    if any(isnothing, idcs_out)
        idcs_in  = idcs_in[idcs_out .!= nothing]
        idcs_out = idcs_out[idcs_out .!= nothing]
    end

    idcs_in, idcs_out
end

"""
Return a sparse matrix that maps quantities given on `basis_in` and `kpt_in`
to quantities on `basis_out` and `kpt_out`.
"""
function compute_transfer_matrix(basis_in::PlaneWaveBasis{T}, kpt_in::Kpoint,
                                 basis_out::PlaneWaveBasis{T}, kpt_out::Kpoint) where {T}
    idcs_in, idcs_out = transfer_mapping(basis_in, kpt_in, basis_out, kpt_out)
    sparse(idcs_out, idcs_in, true)
end

"""
Return a list of sparse matrices (one per ``k``-point) that map quantities given in the
`basis_in` basis to quantities given in the `basis_out` basis.
"""
function compute_transfer_matrix(basis_in::PlaneWaveBasis{T}, basis_out::PlaneWaveBasis{T}) where {T}
    @assert basis_in.model.lattice == basis_out.model.lattice
    @assert length(basis_in.kpoints) == length(basis_out.kpoints)
    @assert all(basis_in.kpoints[ik].coordinate == basis_out.kpoints[ik].coordinate
                for ik in 1:length(basis_in.kpoints))
    [compute_transfer_matrix(basis_in, kpt_in, basis_out, kpt_out)
     for (kpt_in, kpt_out) in zip(basis_in.kpoints, basis_out.kpoints)]
end

"""
Transfer an array `ψk` defined on basis_in ``k``-point kpt_in to basis_out ``k``-point kpt_out.
"""
function transfer_blochwave_kpt(ψk_in, basis_in::PlaneWaveBasis{T}, kpt_in::Kpoint,
                                basis_out::PlaneWaveBasis{T}, kpt_out::Kpoint) where {T}
    kpt_in == kpt_out && return copy(ψk_in)
    @assert length(G_vectors(basis_in, kpt_in)) == size(ψk_in, 1)
    idcsk_in, idcsk_out = transfer_mapping(basis_in, kpt_in, basis_out, kpt_out)

    n_bands = size(ψk_in, 2)
    ψk_out  = similar(ψk_in, length(G_vectors(basis_out, kpt_out)), n_bands)
    ψk_out .= 0
    ψk_out[idcsk_out, :] .= ψk_in[idcsk_in, :]

    ψk_out
end

"""
Transfer an array `ψk_in` expanded on `kpt_in`, and produce ``ψ(r) e^{i ΔG·r}`` expanded on
`kpt_out`. It is mostly useful for phonons.
Beware: `ψk_out` can lose information if the shift `ΔG` is large or if the `G_vectors`
differ between `k`-points.
"""
function transfer_blochwave_kpt(ψk_in, basis::PlaneWaveBasis, kpt_in, kpt_out, ΔG)
    ψk_out = zeros(eltype(ψk_in), length(G_vectors(basis, kpt_out)), size(ψk_in, 2))
    for (iG, G) in enumerate(G_vectors(basis, kpt_in))
        # e^i(kpt_in + G)r = e^i(kpt_out + G')r, where
        # kpt_out + G' = kpt_in + G = kpt_out + ΔG + G, and
        # G' = G + ΔG
        idx_Gp_in_kpoint = index_G_vectors(basis, kpt_out, G - ΔG)
        if !isnothing(idx_Gp_in_kpoint)
            ψk_out[idx_Gp_in_kpoint, :] = ψk_in[iG, :]
        end
    end
    ψk_out
end

"""
Transfer Bloch wave between two basis sets. Limited feature set.
"""
function transfer_blochwave(ψ_in, basis_in::PlaneWaveBasis{T},
                            basis_out::PlaneWaveBasis{T}) where {T}
    @assert basis_in.model.lattice == basis_out.model.lattice
    @assert length(basis_in.kpoints) == length(basis_out.kpoints)
    @assert all(basis_in.kpoints[ik].coordinate == basis_out.kpoints[ik].coordinate
                for ik in 1:length(basis_in.kpoints))

    # If, for some kpt ik, basis_in has less vectors than basis_out, then idcs_out[ik] is
    # the array of the indices of the G_vectors from basis_in in basis_out.
    # It is then of size G_vectors(basis_in.kpoints[ik]) and the transfer can be done with
    # ψ_out[ik] .= 0
    # ψ_out[ik][idcs_out[ik], :] .= ψ_in[ik]

    # Otherwise, if, for some kpt ik, basis_in has more vectors than basis_out, then
    # idcs_out[ik] just keep the indices of the G_vectors from basis_in that are in basis_out.
    # It is then of size G_vectors(basis_out.kpoints[ik]) and the transfer can be done with
    # ψ_out[ik] .= ψ_in[ik][idcs_in[ik], :]

    map(enumerate(basis_out.kpoints)) do (ik, kpt_out)
        transfer_blochwave_kpt(ψ_in[ik], basis_in, basis_in.kpoints[ik], basis_out, kpt_out)
    end
end

@doc raw"""
Transfer density (in real space) between two basis sets.

This function is fast by transferring only the Fourier coefficients from the small basis
to the big basis.

Note that this implies that for even-sized small FFT grids doing the
transfer small -> big -> small is not an identity (as the small basis has an unmatched
Fourier component and the identity ``c_G = c_{-G}^\ast`` does not fully hold).

Note further that for the direction big -> small employing this function does not give
the same answer as using first `transfer_blochwave` and then `compute_density`.
"""
function transfer_density(ρ_in, basis_in::PlaneWaveBasis{T},
                          basis_out::PlaneWaveBasis{T}) where {T}
    @assert basis_in.model.lattice == basis_out.model.lattice
    @assert length(size(ρ_in)) ∈ (3, 4)

    ρ_freq_in  = fft(basis_in, ρ_in)
    ρ_freq_out = zeros_like(ρ_freq_in, basis_out.fft_size..., size(ρ_in, 4))

    for (block_in, block_out) in transfer_mapping(basis_in, basis_out)
        ρ_freq_out[block_out, :] .= ρ_freq_in[block_in, :]
    end

    irfft(basis_out, ρ_freq_out)
end

"""
Find the equivalent index of the coordinate `kcoord` ∈ ℝ³ in a list `kcoords` ∈ [-½, ½)³.
`ΔG` is the vector of ℤ³ such that `kcoords[index] = kcoord + ΔG`.
"""
function find_equivalent_kpt(basis::PlaneWaveBasis{T}, kcoord, spin; tol=sqrt(eps(T))) where {T}
    kcoord_red = map(kcoord) do k
                     k = mod(k, 1)              # coordinate in [0, 1)³
                     k ≥ 0.5 - tol ? k - 1 : k  # coordinate in [-½, ½)³
                 end

    ΔG = kcoord_red - kcoord
    # ΔG should be an integer.
    @assert all(is_approx_integer.(ΔG))
    ΔG = round.(Int, ΔG)

    indices_σ = krange_spin(basis, spin)
    kcoords_σ = getfield.(basis.kpoints[indices_σ], :coordinate)
    # Unique by construction.
    index::Int = findfirst(isapprox(kcoord_red; atol=tol), kcoords_σ) + (indices_σ[1] - 1)

    return (; index, ΔG)
end

"""
Return the Fourier coefficients for `ψk · e^{i q·r}` in the basis of `kpt_out`, where `ψk`
is defined on a basis `kpt_in`.
"""
function multiply_by_expiqr(basis, kpt_in, q, ψk)
    shifted_kcoord = kpt_in.coordinate .+ q  # coordinate of ``k``-point in ℝ
    index, ΔG = find_equivalent_kpt(basis, shifted_kcoord, kpt_in.spin)
    kpt_out = basis.kpoints[index]
    return transfer_blochwave_kpt(ψk, basis, kpt_in, kpt_out, ΔG)
end

"""
Return the indices of the `kpoints` shifted by `q`. That is for each `kpoint` of the `basis`:
`kpoints[ik].coordinate + q = kpoints[indices[ik]].coordinate`.
"""
function k_to_kpq_mapping(basis::PlaneWaveBasis, qcoord)
    kpoints = basis.kpoints
    indices = [find_equivalent_kpt(basis, kpt.coordinate + qcoord, kpt.spin).index
               for kpt in kpoints]
    @assert sort(indices) == 1:length(basis.kpoints)
    indices
end
