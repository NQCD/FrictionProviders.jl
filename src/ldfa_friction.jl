struct LDFAFriction{D,T,S} <: NQCModels.FrictionModels.ElectronicFrictionProvider
    "Density model"
    density::D
    "Temporary array for storing the electron density."
    rho::Vector{T}
    "Temporary array for the Wigner-Seitz radii."
    radii::Vector{T}
    "Splines fitted to the numerical LDFA data."
    splines::S
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
    "Degrees of freedom for each atom. (Should be 3)"
    ndofs::Int
end

"""
    LDFAFriction(density, atoms; friction_atoms=collect(Int, range(atoms)))

FrictionProvider for the Local Density Friction Approximation. 
    
    # Arguments
    
    ## density
    This should be an interface to the prediction of the electronic density as a function of atomic positions. 
    It should support a `density` method. 
    
    ## atoms
    `NQCBase.Atoms` for the structure in question. This determines the correct fitting curve between Wigner-Seitz radii and electronic friction as shown in Gerrits2020. 
    Usually, this is the same `Atoms` object as you use in your dynamics simulation. 

    ## friction_atoms
    
    The atom indices which electronic friction should be applied to. 

"""
function LDFAFriction(density, atoms; friction_atoms=collect(Int, range(atoms)))
    ldfa_data, _ = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true)
    r = ldfa_data[:, 1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:, atoms.numbers[i].+1]
        indices = η .!= ""
        ri = convert(Vector{Float64}, r[indices])
        η = convert(Vector{Float64}, η[indices])
        push!(ri, 10.0) # Ensure it goes through 0.0 for large r.
        push!(η, 0.0)
        push!(splines, CubicSpline(η, ri; extrapolate=true))
    end

    rho = zeros(length(atoms))
    radii = zero(rho)

    LDFAFriction(density, rho, radii, splines, friction_atoms, 3)
end

"""
    get_friction_matrix(model::LDFAFriction, R::AbstractMatrix)

get_friction_matrix uses the specified density model to predict the friction matrix for `friction_atoms` and return it for just those atoms. 
Units of friction are mass-weighted, and the atomic unit of friction is: [E_h / ħ / m_e]
Convert common friction units such as ps^-1 or meV ps Å^-2 using `UnitfulAtomic.austrip`. 

This behaviour is different to NQCModels.friction!, which returns friction for the whole system, not just `friction_atoms`. 
"""
function get_friction_matrix(model::LDFAFriction, R::AbstractMatrix)
    density!(model.density, model.rho, R, model.friction_atoms)
    clamp!(model.rho, 0, Inf)
    @. model.radii = 1 / cbrt(4 / 3 * π * model.rho)
    if any(model.radii .< 1.5)
        @warn "LDFA model density exceeds the fitting regime from Gerrits2020 (rₛ<1.5). This usually happens when interatomic distances are unusually short and should be investigated. rₛ=1.5 is clamped here so dynamics can continue without errors." maxlog=1
        clamp!(model.radii, 1.5, 10) # Ensure Wigner-Seitz radii are within spline bounds. 
    end
    η(r) = r < 10 ? model.splines[1](r) : 0.0
    eft_diagonal = repeat(η.(model.radii[model.friction_atoms]), inner=NQCModels.ndofs(model)) # Evaluate friction splines for each atom, then repeat the values per atom for all degrees of freedom. 
    return eft_diagonal |> diagm |> Diagonal # Return the EFT as a diagonal Matrix with the LinearAlgebra type for efficiency. 
end

export get_friction_matrix

function FrictionModels.friction!(model::LDFAFriction, F::AbstractMatrix, R::AbstractMatrix)
    friction_atom_indices = friction_matrix_indices(model.friction_atoms, NQCModels.ndofs(model))
    F[friction_atom_indices, friction_atom_indices] .= get_friction_matrix(model, R)
end

# Overload friction function to map friction_atoms to Subsystem indices
function FrictionModels.friction!(system::Subsystem{<:LDFAFriction}, F::AbstractMatrix, R::AbstractMatrix)
    F .= get_friction_matrix(system.model, R)
end
