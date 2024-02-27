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

function LDFAFriction(density, atoms; friction_atoms=collect(range(atoms))) 
    ldfa_data, _ = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true)
    r = ldfa_data[:,1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:,atoms.numbers[i].+1]
        indices = η .!= ""
        ri = convert(Vector{Float64}, r[indices])
        η = convert(Vector{Float64}, η[indices])
        push!(ri, 10.0) # Ensure it goes through 0.0 for large r.
        push!(η, 0.0)
        push!(splines, CubicSpline(η, ri))
    end

    rho = zeros(length(atoms))
    radii = zero(rho)

    LDFAFriction(density, rho ,radii, splines, friction_atoms, 3)
end

function get_friction_matrix(model::LDFAFriction, R::AbstractMatrix)
    density!(model.density, model.rho, R, model.friction_atoms)
    clamp!(model.rho, 0, Inf)
    @. model.radii = 1 / cbrt(4/3 * π * model.rho)
    η(r)=r < 10 ? model.splines[1](r) : 0.0
    return Diagonal(diagm(repeat(η.(model.radii[model.friction_atoms]), inner=NQCModels.ndofs(model))))
end

export get_friction_matrix

function FrictionModels.friction!(model::LDFAFriction, F::AbstractMatrix, R::AbstractMatrix)
    friction_atom_indices=friction_matrix_indices(model.friction_atoms, NQCModels.ndofs(model))
    F[friction_atom_indices, friction_atom_indices] .= get_friction_matrix(model, R)
end

# Overload friction function to map friction_atoms to Subsystem indices
function FrictionModels.friction!(system::Subsystem{<:LDFAFriction}, F::AbstractMatrix, R::AbstractMatrix)
    F .= get_friction_matrix(system.model, R)
end