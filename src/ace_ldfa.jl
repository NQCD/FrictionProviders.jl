"""
This uses a ACE ML models to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

struct AceLDFA{L,U} <: ElectronDensityProvider
    friction_IP::L
    "Units"
    density_unit::U
end

"""
    AceLDFA(friction_IP; density_unit=u"Å^-3")

ACE.jl electron density friction provider. 
    
    # Arguments
    
    ## friction_IP
    This is the JuLIP potential created from the fitted ACE model, which is used to predict the electron density as a function of atomic positions. 
        
    ## density_unit
    
    The unit in which the ACE model predicts the electron density. By default, this is in Å^-3 and should be converted to atomic units of Bohr^-3. 

"""
function AceLDFA(friction_IP; density_unit=u"Å^-3")
    AceLDFA(friction_IP, density_unit)
end


function density!(model::AceLDFA, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    # Needs a structure with only one "friction_atom" atom type for evaluation. 
    one_frictionatom_indices = length(friction_atoms)≥2 ? symdiff(1:size(R,2), friction_atoms[2:end]) : 1:size(R,2)
    for i in friction_atoms
        xR = R[:,one_frictionatom_indices] 
        xR[:,end]=R[:,i]
        rho[i] = austrip(au_to_eV(NQCModels.potential(model.friction_IP, xR)) * model.density_unit)
    end
end
