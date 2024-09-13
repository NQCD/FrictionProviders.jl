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
        push!(splines, CubicSpline(η, ri; extrapolate=true))
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
    DoFs = size(R, 1)
    for i in model.friction_atoms
        η = model.radii[i] < 10 ? model.splines[i](model.radii[i]) : 0.0
        for j in axes(R, 1)
            F[(i-1)*DoFs+j, (i-1)*DoFs+j] = η
        end
    end
    return F
end
