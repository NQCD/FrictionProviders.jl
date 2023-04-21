struct ODFriction{M,T,V}
    "Friction model"
    friction::M
    "Temporary vector for storing friction output of model"
    f::T
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
    cutoff::V
end

function ODFriction(friction, atoms; friction_atoms=collect(range(atoms)), cutoff=[Inf])
    
    nfrict = length(friction_atoms)
    f = zeros(nfrict*3,nfrict*3)

    ODFriction(friction, f, friction_atoms, cutoff)
end

NQCModels.ndofs(model::ODFriction) = 3


function NQCModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)

    # for i in model.friction_atoms

    
    skfriction!(model.friction,R, model.f, model.friction_atoms)

    counter = 1
    for i in model.friction_atoms
        for j=1:3
            F[(i-1)*3+j, (i-1)*3+j] = model.f[(counter-1)*3+j, (counter-1)*3+j]

            if R[3,i] > model.cutoff[i]
                F[(i-1)*3+j, (i-1)*3+j] = 0.0
            end
        end
        counter = counter + 1
    end
    
    return F
end
