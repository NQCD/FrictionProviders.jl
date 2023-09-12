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
    for i,f_a in enumerate(friction_atoms)
        xR = R[:,1:end-length(friction_atoms)+1]
        xR[:,end]=R[:,f_a]
        rho[i] = austrip(NQCModels.potential(model.ml_model, xR) * model.density_unit)
    end
end
