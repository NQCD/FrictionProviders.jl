"""
This uses a SciKit ML models to attach friction coefficients to existing models, where the Sk model
provides the friction coefficients directly. (unlike models that provide the density)
"""

struct SciKitFriction{D,L,A,S,U,E}
    "Descriptors used for scikit models"
    descriptors::D
    "Sci-Kit ML model"
    ml_model::L
    "Atoms"
    atoms::A
    "Scaler for descriptor"
    scaler::S
    "Units"
    friction_unit::U
    "Indices of atoms for descriptor"
    descriptor_atoms::E
end

function SciKitFriction(descriptors, ml_model, atoms, scaler; friction_unit=u"ps^-1", descriptor_atoms=[0])
    SciKitFriction(descriptors, ml_model, atoms, scaler, friction_unit, descriptor_atoms)
end


function set_coordinates!(model::SciKitFriction, R)
    model.atoms.set_positions(ustrip.(auconvert.(u"Å", R')))
end

function skfriction!(model::SciKitFriction,f::AbstractMatrix, R::AbstractMatrix, friction_atoms::AbstractVector, )

    set_coordinates!(model, R)
    r_desc = model.descriptors.create(atoms, positions=model.descriptor_atoms.-1, n_jobs=-1) #n_threads)
    f = austrip(model.ml_model.predict(r_desc) * model.friction_unit)

    f
end

"Need to convert atoms to SOAP descriptor"

"Need to load scaling from python? and apply same"

"Load model from python?"

"Predict friction from model"

"Need to figure out if its giving just diagonal or full tensor"

