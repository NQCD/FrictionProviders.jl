"""
This uses a SciKit ML models to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

struct SciKitLDFA{D,L,A,U,S}
    "Descriptors used for scikit models"
    descriptors::D
    "Sci-Kit ML model"
    ml_model::L
    "Atoms"
    atoms_ase::A
    "Units"
    density_unit::U
    "scaler"
    scaler::S
end

function SciKitLDFA(descriptors, ml_model, atoms_ase; density_unit=u"Å^-3", scaler=nothing)
    SciKitLDFA(descriptors, ml_model, atoms_ase, density_unit, scaler)
end


function density!(model::SciKitLDFA, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    for i in friction_atoms
        model.atoms_ase.set_positions(au_to_ang.(R'))
        density_atoms = model.atoms_ase.copy()
        friction_atoms_srtd = sort(friction_atoms, rev=true)
        for j=1:length(friction_atoms_srtd)
            density_atoms.pop(i=friction_atoms_srtd[j]-1)
        end
        density_atoms.append(model.atoms_ase[i])
        r_desc = model.descriptors.create(density_atoms, centers=[length(density_atoms)-1], n_jobs=1) #n_threads)
        if model.scaler != nothing
            r_desc = model.scaler.transform(r_desc)
        end
        rho[i] = austrip(model.ml_model.predict(r_desc)[end] * model.density_unit)
    end
end

