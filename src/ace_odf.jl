"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct ACEdsODF{M,G,A,U}
    "ACE/friction_tensor calculator"
    model::M
    "Gamma ACEds function"
    gamma::G
    "JuLIP atoms object"
    atoms_julip::A
    "Units"
    friction_unit::U
end

function ACEdsODF(model, gamma, atoms_julip; friction_unit=u"ps^-1")
    ACEdsODF(model, gamma, atoms_julip, friction_unit)
end

function friction!(model::ACEdsODF, R::AbstractMatrix, friction::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)
    set_positions!(model.atoms_julip, au_to_ang.(R))
    DoFs = size(R, 1)
    friction .= reinterpret(Matrix,Matrix(model.gamma(model.model, model.atoms_julip)[friction_atoms, friction_atoms]))
    mass_weights = zeros(length(friction_atoms)*DoFs,length(friction_atoms)*DoFs)
    for fx in 1:size(mass_weights,1)
        for fy in 1:size(mass_weights,2)
            mass_weights[fx,fy] = sqrt(model.atoms_julip.M[friction_atoms[Int(ceil(fx/DoFs,digits=0))]])*sqrt(model.atoms_julip.M[friction_atoms[Int(ceil(fy/DoFs,digits=0))]])
        end
    end
    friction .= austrip.(friction .* model.friction_unit)
    friction .*= mass_weights
    friction .= austrip.(friction .* u"u")
end
