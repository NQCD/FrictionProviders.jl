struct ODFriction{C,F,D}
    "Friction model"
    model::C
    "Temporary vector for storing friction output of model"
    friction::F
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
    "Cutoff distance"
    cutoff::D
end

function ODFriction(model; friction_atoms=collect(range(atoms)), cutoff=Float64(5.0))
    natoms = length(friction_atoms)
    friction = zeros(natoms*3,natoms*3)

    ODFriction(model, friction, friction_atoms, cutoff)
end

NQCModels.ndofs(model::ODFriction) = 3

function FrictionModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)
    friction!(model.model, R, model.friction, model.friction_atoms, model.cutoff)
    
    DoFs = size(R, 1)
    F[(model.friction_atoms[1]-1)*DoFs+1:(model.friction_atoms[2])*DoFs, (model.friction_atoms[1]-1)*DoFs+1:(model.friction_atoms[2])*DoFs] = model.friction
    
    return F
end