struct ODFriction{T}
    "Temporary vector for storing friction output of model"
    f::T
    "Indices of atoms that should have friction applied."
    friction_atoms::Vector{Int}
end

function ODFriction(atoms; friction_atoms=collect(range(atoms)))
    
    ODFriction(friction_atoms)
end

NQCModels.ndofs(model::ODFriction) = 3



function FrictionModels.friction!(model::ODFriction, F::AbstractMatrix, R::AbstractMatrix)
    skfriction!(R,f, model.friction_atoms)

    println(size(f))
    # for i in model.friction_atoms
    #     for j in axes(R, 1)
    #         F[(i-1)*DoFs+j, (i-1)*DoFs+j] = f
    #     end
    # end
    return F
end
