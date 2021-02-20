# # Input and output formats
#md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/guide/@__NAME__.ipynb)
#md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/guide/@__NAME__.ipynb)

# This section provides an overview of the input and output formats
# supported by DFTK, usually via integration with a third-party library.
#
# ## Reading input formats supported by ASE
# ASE is short for the
# [atomistic simulation environment](https://wiki.fysik.dtu.dk/ase/index.html),
# a Python package to simplify the process of setting up, running and
# analysing results from atomistic simulations across different programs.
# If it is installed it is automatically used by DFTK in order to read a wide range
# of input files, such as xyz, CIF, input files of various other codes
# (e.g. Quantum Espresso, VASP, ABINIT, CASTEP, ...).
# The full list of formats
# can be found in the [ASE IO documentation](https://wiki.fysik.dtu.dk/ase/ase/io/io.html).
#
# As an example we start the calculation of a simple antiferromagnetic iron crystal
# using a Quantum-Espresso input file, [Fe_afm.pwi](Fe_afm.pwi).
# From this file the lattice, atomic positions and the initial magnetisation are read.
# For more details about calculations on magnetic systems
# using collinear spin, see [Collinear spin and magnetic systems](@ref).
#
# First we parse the Quantum Espresso input file to DFTK datastructures:

using DFTK

qe_input = "Fe_afm.pwi"
atoms    = load_atoms(qe_input)
lattice  = load_lattice(qe_input)
magnetic_moments = load_magnetic_moments(qe_input)

# At this point a file of any format supported by ASE could be passed instead,
# e.g. an `xyz` file or an ABINIT input file. Behind the scenes ASE takes care
# to select the right parser and extract the required structural information.

# Next we attach the pseudopotential information, since this information is currently
# not exposed inside the ASE datastructures.
# See [Creating slabs with ASE](@ref) for more details.

atoms = [ElementPsp(el.symbol, psp=load_psp(el.symbol, functional="pbe")) => position
         for (el, position) in atoms];

# Finally we run the calculation.

model = model_LDA(lattice, atoms, magnetic_moments=magnetic_moments, temperature=0.01)
basis = PlaneWaveBasis(model, 10; kgrid=(2, 2, 2))
ρ0 = guess_density(basis, magnetic_moments)
scfres = self_consistent_field(basis, tol=1e-4, ρ=ρ0, mixing=KerkerMixing());

# !!! note "DFTK and ASE"
#     DFTK also supports using ASE to setup a calculation
#     and provides two-way conversion routines between key DFTK
#     and ASE datastructures. See [Creating slabs with ASE](@ref) for details.
#
# ## Writing VTK files for visualization
# For visualizing the density or the Kohn-Sham orbitals DFTK supports storing
# the result of an SCF calculations in the form of VTK files.
# These can afterwards be visualized using tools such
# as [paraview](https://www.paraview.org/).
# Using this feature requires
# the [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl/) Julia package.

using WriteVTK
save_scfres("iron_afm.vts", scfres; save_ψ=true);

# This will save the iron calculation above into the file `iron_afm.vts`,
# using `save_ψ=true` to also include the KS orbitals.

# ## Writing and reading JLD2 files
# The full state of a DFTK self-consistent field calculation can be
# stored on disk in form of an [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) file.
# This file can be read from other Julia scripts
# as well as other external codes supporting the HDF5 file format
# (since the JLD2 format is based on HDF5).

using JLD2
save_scfres("iron_afm.jld2", scfres);

# Since such JLD2 can also be read by DFTK to start or continue a calculation,
# these can also be used for checkpointing or for transferring results
# to a different computer.
# See [Saving SCF results on disk and SCF checkpoints](@ref) for details.

# (Cleanup files generated by this notebook)
rm("iron_afm.vts")
rm("iron_afm.jld2")
