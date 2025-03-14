# TODO: Temporary, explanations too scarce. To be changed with proper phonon computations.
using Test
using DFTK
using DFTK: energy_forces_pairwise
using LinearAlgebra
using ForwardDiff
using StaticArrays
using Random

@testset "Phonons" begin

# Convert back and forth between Vec3 and columnwise matrix
fold(x)   = hcat(x...)
unfold(x) = Vec3.(eachcol(x))

function prepare_system(; n_scell=1)
    positions = [Vec3([0.0, 0.0, 0.0])]
    for i in 1:n_scell-1
        push!(positions, Vec3(i * ones(3) / n_scell))
    end

    a = 5. * length(positions)
    lattice = a * [[1 0 0.]; [0 0 0.]; [0 0 0.]]

    Linv = DFTK.compute_inverse_lattice(lattice)
    directions = [[Linv * [i==j, 0.0, 0.0] for i in 1:n_scell] for j in 1:n_scell]

    params = Dict((:X, :X) => (; ε=1, σ=a / length(positions) /2^(1/6)))
    V(x, p) = 4*p.ε * ((p.σ/x)^12 - (p.σ/x)^6)

    (; positions, lattice, directions, params, V)
end

function prepare_3d_system(; n_scell=1)
    positions = [[0.0, 0.0, 0.0]]
    for i in 1:n_scell-1
        push!(positions, i * ones(3) / n_scell)
    end

    a = 5. * length(positions)
    lattice = a * rand(3, 3)

    (; positions, lattice)
end

# Compute phonons for a one-dimensional pairwise potential for a set of `q = 0` using
# supercell method
function test_supercell_q0(; n_scell=1, max_radius=1e3)
    case = prepare_system(; n_scell)
    Linv  = DFTK.compute_inverse_lattice(case.lattice)
    n_atoms = length(case.positions)

    directions = [reshape(vcat([[i==j, 0.0, 0.0] for i in 1:n_atoms]...), 3, :)
                  for j in 1:n_atoms]

    Φ = zeros(length(directions), n_atoms)
    for (i, direction) in enumerate(directions)
        Φ[i, :] = - ForwardDiff.derivative(0.0) do ε
            new_positions = unfold(fold(case.positions) .+ ε .* Linv * direction)
            forces = energy_forces_pairwise(case.lattice, fill(:X, n_atoms), new_positions,
                                            case.V, case.params; max_radius).forces
            [(Linv * f)[1] for f in forces]
        end
    end
    sqrt.(abs.(eigvals(Φ)))
end

# Compute phonons for a one-dimensional pairwise potential for a set of `q`-points
function test_ph_disp(; n_scell=1, max_radius=1e3, n_points=2)
    case = prepare_system(; n_scell)
    Linv  = DFTK.compute_inverse_lattice(case.lattice)
    n_atoms = length(case.positions)

    function pairwise_ph(q, d)
        energy_forces_pairwise(case.lattice, fill(:X, n_atoms), case.positions, case.V,
                               case.params; q=[q, 0, 0], ph_disp=d, max_radius).forces
    end

    ph_bands = []
    for q in range(-1/2, 1/2, length=n_points+1)
        as = ComplexF64[]
        for d in case.directions
            res = -ForwardDiff.derivative(0.0) do ε
                forces = pairwise_ph(q, ε*d)
                [Linv' * f for f in forces]
            end
            [push!(as, r[1]) for r in res]
        end
        M = reshape(as, n_atoms, :)
        @assert norm(imag.(eigvals(M))) < 1e-8
        push!(ph_bands, sqrt.(abs.(real(eigvals(M)))))
    end
    ph_bands
end

"""
Real-space equivalent of `transfer_blochwave_kpt`.
"""
function transfer_blochwave_kpt_real(ψk_in, basis::PlaneWaveBasis, kpt_in, kpt_out, ΔG)
    ψk_out = zeros(eltype(ψk_in), length(kpt_out.G_vectors), size(ψk_in, 2))
    exp_ΔGr = DFTK.cis2pi.(-dot.(Ref(ΔG), r_vectors(basis)))
    for n in 1:size(ψk_in, 2)
        ψk_out[:, n] = fft(basis, kpt_out, exp_ΔGr .* ifft(basis, kpt_in, ψk_in[:, n]))
    end
    ψk_out
end


@testset "Phonon consistency" begin
    tol = 1e-6
    n_points = 10
    max_radius = 1e3

    ph_bands = [test_ph_disp(; n_scell, max_radius, n_points) for n_scell in 1:3]

    # Recover the same extremum for the system whatever case we test
    for n_scell in 2:3
        @test minimum(fold(ph_bands[1])) ≈ minimum(fold(ph_bands[n_scell])) atol=tol
        @test maximum(fold(ph_bands[1])) ≈ maximum(fold(ph_bands[n_scell])) atol=tol
    end

    # Test consistency between supercell method at `q = 0` and direct `q`-points computations
    for n_scell in 1:3
        r_q0 = test_supercell_q0(; n_scell, max_radius)
        @assert length(r_q0) == n_scell
        ph_band_q0 = ph_bands[n_scell][n_points÷2+1]
        @test norm(r_q0 - ph_band_q0) < tol
    end
end

@testset "Shifting functions" begin
    Random.seed!()
    tol = 1e-12

    case = prepare_3d_system()

    X = ElementGaussian(1.0, 0.5, :X)
    atoms = [X for _ in case.positions]
    n_atoms = length(case.positions)

    model = Model(case.lattice, atoms, case.positions; n_electrons=n_atoms,
                  symmetries=false, spin_polarization=:collinear)
    kgrid = rand(2:10, 3)
    k1, k2, k3 = kgrid
    basis = PlaneWaveBasis(model; Ecut=100, kgrid)

    # We consider a smooth periodic function with Fourier coefficients given if the basis
    # e^(iG·x)
    ψ = rand(ComplexF64, size(r_vectors(basis)))

    # Random `q` shift
    q0 = rand(basis.kpoints).coordinate
    ishift = [rand(-k1*2:k1*2), rand(-k2*2:k2*2), rand(-k3*2:k3*2)]
    q = Vec3(q0 .* ishift)
    @testset "Transfer function" begin
        for kpt in unique(rand(basis.kpoints, 4))
            ψk = fft(basis, kpt, ψ)

            ψk_out_four = DFTK.multiply_by_expiqr(basis, kpt, q, ψk)
            ψk_out_real = let
                shifted_kcoord = kpt.coordinate .+ q
                index, ΔG = DFTK.find_equivalent_kpt(basis, shifted_kcoord, kpt.spin)
                kpt_out = basis.kpoints[index]
                transfer_blochwave_kpt_real(ψk, basis, kpt, kpt_out, ΔG)
            end
            @testset "Testing kpoint $(kpt.coordinate) on kgrid $kgrid" begin
                @test norm(ψk_out_four - ψk_out_real) < tol
            end
        end
    end

    @testset "Ordering function" begin
        kpoints_plus_q = DFTK.k_to_kpq_mapping(basis, q)
        ordering(kdata) = kdata[kpoints_plus_q]
        kcoords = getfield.(basis.kpoints, :coordinate)
        for (ik, kcoord) in enumerate(kcoords)
            @test mod.(kcoord + q .- tol, 1) ≈ mod.(ordering(kcoords)[ik] .- tol, 1)
        end
    end
end

end
