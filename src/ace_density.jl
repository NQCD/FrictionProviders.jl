"""
This uses a ACE ML models to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

struct ACEDensity{L,U}
    "Sci-Kit ML model"
    ml_model::L
    "Units"
    density_unit::U
end

function ACEDensity(ml_model; density_unit=u"Å^-3")
    ACEDensity(ml_model, density_unit)
end


function density!(model::ACEDensity, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    for i in friction_atoms
        xR = R[:,[1:end-2,i]]
        rho[i] = austrip(NQCModels.potential(model.ml_model, xR) * model.density_unit)
    end
end