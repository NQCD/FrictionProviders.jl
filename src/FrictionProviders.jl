module FrictionProviders

using DataInterpolations: CubicSpline
using DelimitedFiles: readdlm
using UnitfulAtomic: austrip, auconvert
using Unitful: @u_str, ustrip
using NQCBase: PeriodicCell, apply_cell_boundaries!, au_to_ang, au_to_eV
using NQCModels: NQCModels, FrictionModels, Model, Subsystem
using JuLIP: set_positions!
using LinearAlgebra

"""
    friction_matrix_indices(model, indices)

Returns the indices of the friction matrix corresponding to the given Atom indices. 
"""
function friction_matrix_indices(indices, dofs)
	dof_range=collect(1:dofs)
	return vcat(broadcast(x->x.+dof_range, dofs .* (indices .- 1))...)  
end

include("ldfa_friction.jl")
include("cube.jl")
include("cube_ldfa.jl")
include("scikit_ldfa.jl")
include("ace_ldfa.jl")
include("odf_friction.jl")
include("schnet_odf.jl")
include("ace_odf.jl")

export LDFAFriction
export AceLDFA
export CubeLDFA
export SciKitLDFA
export ODFriction
export SchNetODF
export ACEdsODF


end
