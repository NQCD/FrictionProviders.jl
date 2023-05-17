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

function ODFriction(model; friction_atoms=collect(range(atoms)), cutoff=5.0)
    nfrict = length(friction_atoms)
    friction = zeros(nfrict*3,nfrict*3)

    ODFriction(model, friction, friction_atoms, cutoff)
end

NQCModels.ndofs(model::ODFriction) = 3

function FrictionModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)
    friction!(model.model, R, model.friction, model.friction_atoms, model.cutoff)

    DoFs = size(R, 1)
    for i in axes(R, 1)
        F[(model.friction_atoms[1]-1)*DoFs+i, (model.friction_atoms[2]-1)*DoFs+i] = model.friction
    end
    return F
end