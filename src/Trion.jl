module Trion 

using SpecialFunctions
using QuadGK
using LinearAlgebra
using Optim

include("Indexing.jl")
include("Potential.jl")
include("Hamiltonian.jl")
include("Spectrum.jl")
include("Nonlinearity.jl")

end