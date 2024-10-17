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

    for i in 1:length(model.friction_atoms)
        for j in 1:length(model.friction_atoms)
            atom_i = model.friction_atoms[i]
            atom_j = model.friction_atoms[j]
            
            F[(atom_i-1)*DoFs+1:atom_i*DoFs, (atom_j-1)*DoFs+1:atom_j*DoFs] = model.friction[(i-1)*DoFs+1:i*DoFs, (j-1)*DoFs+1:j*DoFs]
        end
    end

    
    return F
end