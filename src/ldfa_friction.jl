struct LDFAFriction{D,T,S}
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
end

"""
    LDFAFriction(density, atoms; friction_atoms=collect(range(atoms)))

Constructor for an LDFA friction model. This type of model uses the relationship between Wigner-Seitz radii and electron densities to predict electronic friction coefficients based on the electron density. 

## Arguments

1. `density` - An electron density model that can be evaluated using the `density!` function.

2. `atoms` - Atom types in the structure the model is used for. This determines the relationship between electron density and electronic friction coefficient. 

3. `friction_atoms` - Atom indices for which the electronic friction coefficent should be returned. All other atoms will return η=0. 

"""
function LDFAFriction(density, atoms; friction_atoms=collect(range(atoms)))
    ldfa_data, _ = readdlm(joinpath(@__DIR__, "ldfa.txt"), ',', header=true) # electron density to electronic friction coefficent relationships taken from Gerrits et al, 2020: https://doi.org/10.1103/PhysRevB.102.155130
    r = ldfa_data[:,1]
    splines = []
    for i in range(atoms)
        η = ldfa_data[:,atoms.numbers[i].+1]
        indices = η .!= ""
        ri = convert(Vector{Float64}, r[indices])
        η = convert(Vector{Float64}, η[indices])
        push!(ri, 10.0) # Ensure it goes through 0.0 for large r.
        push!(η, 0.0)
        push!(splines, CubicSpline(η, ri; extrapolate=false))
    end

    rho = zeros(length(atoms))
    radii = zero(rho)

    LDFAFriction(density, rho ,radii, splines, friction_atoms)
end

NQCModels.ndofs(model::LDFAFriction) = 3

function FrictionModels.friction!(model::LDFAFriction, F::AbstractMatrix, R::AbstractMatrix)
    density!(model.density, model.rho, R, model.friction_atoms)
    clamp!(model.rho, 0, Inf)
    @. model.radii = 1 / cbrt(4/3 * π * model.rho)
    if any(model.radii .< 1.5)
        @warn "LDFA model density exceeds the fitting regime from Gerrits2020 (rₛ<1.5). This usually happens when interatomic distances are unusually short and should be investigated. rₛ=1.5 is clamped here so dynamics can continue without errors." maxlog=1
        clamp!(model.radii, 1.5, 10) # Ensure Wigner-Seitz radii are within spline bounds. 
    end
    DoFs = size(R, 1)
    for i in model.friction_atoms
        η = model.splines[i](model.radii[i])
        for j in axes(R, 1)
            F[(i-1)*DoFs+j, (i-1)*DoFs+j] = η
        end
    end
    return F
end
