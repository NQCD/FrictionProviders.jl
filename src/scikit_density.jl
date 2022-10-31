"""
This uses a SciKit ML models to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

struct SciKitDensity{D,L,A}
    "Descriptors used for scikit models"
    descriptors::D
    "Loaded ML model"
    loaded_model::L
    "Atoms"
    atoms::A
end

function set_coordinates!(model::SciKitDensity, R)
    model.atoms.set_positions(ustrip.(auconvert.(u"Å", R')))
end

function density!(model::SciKitDensity, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    for i in friction_atoms
        set_coordinates!(model, R)
        r_desc = model.descriptors.create(model.atoms, positions=[i-1], n_jobs=1) #n_threads)
        rho[i] = model.loaded_model.predict(r_desc)[end]
    end
end
