struct ODFriction{M,T}
    "Friction model"
    friction::M
    "Temporary vector for storing friction output of model"
    f::T
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
end

function ODFriction(friction, atoms; friction_atoms=collect(range(atoms)))
    
    nfrict = length(friction_atoms)
    f = zeros(nfrict*3,nfrict*3)

    ODFriction(friction, f, friction_atoms)
end

NQCModels.ndofs(model::ODFriction) = 3



function FrictionModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)
    skfriction!(model.friction,R, model.f, model.friction_atoms)

    # for i in model.friction_atoms
    #     for j in axes(R, 1)
    #         F[(i-1)*DoFs+j, (i-1)*DoFs+j] = f
    #     end
    # end

    counter = 1
    for i in model.friction_atoms
        for j=1:3
            F[(i-1)*3+j, (i-1)*3+j] = model.f[(counter-1)*3+j, (counter-1)*3+j]
        end
        counter = counter + 1
    end

    return F
end
