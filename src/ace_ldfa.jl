"""
This uses a ACE ML models to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

struct AceLDFA{L,U}
    "Sci-Kit ML model"
    ml_model::L
    "Units"
    density_unit::U
end

function AceLDFA(ml_model; density_unit=u"Å^-3")
    AceLDFA(ml_model, density_unit)
end


function density!(model::AceLDFA, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    for i in friction_atoms
        xR = R[:,1:end-length(friction_atoms)+1] # we assume that friction atoms are at the end
        xR[:,end]=R[:,i]
        rho[i] = austrip(au_to_eV(NQCModels.potential(model.ml_model, xR)) * model.density_unit)
    end
end
