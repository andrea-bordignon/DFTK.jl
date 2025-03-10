using DFTK
import DFTK: mpi_mean!
using Test
using Random
using MPI
include("testcases.jl")

# TODO: refactor tests to avoid code repetition
@testset "Forces on silicon" begin
    function energy_forces(positions)
        model = model_DFT(silicon.lattice, silicon.atoms, positions, [:lda_x, :lda_c_pw])
        basis = PlaneWaveBasis(model; Ecut=7, kgrid=[2, 2, 2], kshift=[0, 0, 0],
                               symmetries_respect_rgrid=true,
                               fft_size=(18, 18, 18))  # FFT chosen to match QE
        is_converged = DFTK.ScfConvergenceDensity(1e-11)
        scfres = self_consistent_field(basis; is_converged)
        scfres.energies.total, compute_forces(scfres), compute_forces_cart(scfres)
    end

    # symmetrical positions, forces should be 0
    _, F0, _ = energy_forces([(ones(3)) / 8, -ones(3) / 8])
    @test norm(F0) < 1e-4

    pos1 = [([1.01, 1.02, 1.03]) / 8, -ones(3) / 8]  # displace a bit from equilibrium
    disp = rand(3)
    mpi_mean!(disp, MPI.COMM_WORLD)  # must be identical on all processes
    ε = 1e-5
    pos2 = [pos1[1] + ε * disp, pos1[2]]
    pos3 = [pos1[1] - ε * disp, pos1[2]]

    # second-order finite differences for accurate comparison
    # TODO switch the other tests to this too
    E1, F1, Fc1 = energy_forces(pos1)
    E2,  _,  _  = energy_forces(pos2)
    E3,  _,  _  = energy_forces(pos3)

    diff_findiff = -(E2 - E3) / (2ε)
    diff_forces = dot(F1[1], disp)
    @test abs(diff_findiff - diff_forces) < 1e-7

    # Test against QE v7.1 using LibXC v5.3.2 lda_x+lda_c_pw
    # (see testcases_QuantumESPRESSO/silicon_LDA_forces)
    reference = [[-5.809880257762980e-3, -4.601261807107919e-3, -3.374609176201532e-3],
                 [ 5.809880257762980e-3,  4.601261807107919e-3,  3.374609176201532e-3]]
    @test maximum(v -> maximum(abs, v), reference - Fc1) < 1e-5
end

@testset "Forces on silicon with non-linear core correction" begin
    function energy_forces(positions)
        Si = ElementPsp(silicon.atnum, :Si, load_psp(silicon.psp_upf))
        atoms = fill(Si, length(silicon.atoms))
        model = model_DFT(silicon.lattice, atoms, positions, [:lda_x, :lda_c_pw])
        basis = PlaneWaveBasis(model; Ecut=7, kgrid=[2, 2, 2], kshift=[0, 0, 0],
                               symmetries_respect_rgrid=true,
                               fft_size=(18, 18, 18))  # FFT chosen to match QE
        is_converged = DFTK.ScfConvergenceDensity(1e-11)
        scfres = self_consistent_field(basis; is_converged)
        scfres.energies.total, compute_forces(scfres), compute_forces_cart(scfres)
    end

    # symmetrical positions, forces should be 0
    _, F0, _ = energy_forces([(ones(3)) / 8, -ones(3) / 8])
    @test norm(F0) < 1e-4

    pos1 = [([1.01, 1.02, 1.03]) / 8, -ones(3) / 8]  # displace a bit from equilibrium
    disp = rand(3)
    mpi_mean!(disp, MPI.COMM_WORLD)  # must be identical on all processes
    ε = 1e-5
    pos2 = [pos1[1] + ε * disp, pos1[2]]
    pos3 = [pos1[1] - ε * disp, pos1[2]]

    # second-order finite differences for accurate comparison
    # TODO switch the other tests to this too
    E1, F1, Fc1 = energy_forces(pos1)
    E2,  _,  _  = energy_forces(pos2)
    E3,  _,  _  = energy_forces(pos3)

    diff_findiff = -(E2 - E3) / (2ε)
    diff_forces = dot(F1[1], disp)
    @test abs(diff_findiff - diff_forces) < 1e-7

    # Test against Abinit v9.6.2 using LibXC v4.3.2 lda_x+lda_c_pw
    # (see testcases_ABINIT/silicon_NLCC_forces)
    reference = [[-0.00574838157984, -0.00455216015517, -0.00333786048065],
                 [ 0.00574838157984,  0.00455216015517,  0.00333786048065]]
    @test maximum(v -> maximum(abs, v), reference - Fc1) < 1e-5
end

@testset "Forces on silicon with spin and temperature" begin
    function silicon_energy_forces(positions; smearing=Smearing.FermiDirac())
        model  = model_DFT(silicon.lattice, silicon.atoms, positions, :lda_xc_teter93;
                           temperature=0.03, smearing, spin_polarization=:collinear)
        basis  = PlaneWaveBasis(model; Ecut=4, kgrid=[4, 1, 2], kshift=[1/2, 0, 0])
        scfres = self_consistent_field(basis; is_converged=DFTK.ScfConvergenceDensity(5e-10))
        scfres.energies.total, compute_forces(scfres)
    end

    pos1 = [([1.01, 1.02, 1.03]) / 8, -ones(3) / 8]  # displace a bit from equilibrium
    disp = rand(3)
    mpi_mean!(disp, MPI.COMM_WORLD)  # must be identical on all processes
    ε = 1e-6
    pos2 = [pos1[1] + ε * disp, pos1[2]]

    for (tol, smearing) in [(0.003, Smearing.FermiDirac()), (5e-5, Smearing.Gaussian())]
        E1, F1 = silicon_energy_forces(pos1; smearing)
        E2, _  = silicon_energy_forces(pos2; smearing)

        diff_findiff = -(E2 - E1) / ε
        diff_forces  = dot(F1[1], disp)
        @test abs(diff_findiff - diff_forces) < tol
    end
end

@testset "Forces on oxygen with spin and temperature" begin
    function oxygen_energy_forces(positions)
        magnetic_moments = [1.0, 1.0]
        model = model_PBE(diagm([7.0, 7.0, 7.0]), o2molecule.atoms, positions;
                          temperature=0.02, smearing=Smearing.Gaussian(), magnetic_moments)
        basis = PlaneWaveBasis(model; Ecut=4, kgrid=[1, 1, 1])

        scfres = self_consistent_field(basis;
                                       is_converged=DFTK.ScfConvergenceDensity(1e-7),
                                       ρ=guess_density(basis, magnetic_moments),
                                       damping=0.7, mixing=SimpleMixing())
        scfres.energies.total, compute_forces(scfres)
    end

    pos1 = [[0, 0, 0.1155], [0.01, -2e-3, -0.2]]
    disp = rand(3)
    mpi_mean!(disp, MPI.COMM_WORLD)  # must be identical on all processes
    ε = 1e-6
    pos2 = [pos1[1] + ε * disp, pos1[2]]

    E1, F1 = oxygen_energy_forces(pos1)
    E2, _  = oxygen_energy_forces(pos2)

    diff_findiff = -(E2 - E1) / ε
    diff_forces  = dot(F1[1], disp)

    @test abs(diff_findiff - diff_forces) < 1e-4
end
