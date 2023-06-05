"""
This script includes using SchNet ML model based on electronic_friction package, to generate electronic friction tensor
"""

struct SchNetODF{C,A,U}
    "SchNet/friction_tensor calculator"
    calculator::C
    "ASE atoms object"
    atoms_ase::A
    "Units"
    friction_unit::U
end

function SchNetODF(calculator, atoms_ase; friction_unit=u"ps^-1")
    SchNetODF(calculator, atoms_ase, friction_unit)
end


function friction!(model::SchNetODF, R::AbstractMatrix, friction::AbstractMatrix, friction_atoms::AbstractVector, cutoff::Float64)
    model.atoms_ase.set_positions(au_to_ang.(R'))
    model.calculator.calculate(model.atoms_ase)
    friction .= model.calculator.get_friction_tensor()
    mass_weights = zeros(length(friction_atoms)*3,length(friction_atoms)*3)
    for fx in 1:size(mass_weights,1)
        for fy in 1:size(mass_weights,2)
            mass_weights = sqrt(model.atoms_ase[friction_atoms[Int(ceil(fx/3,digits=0))]].mass)*sqrt(model.atoms_ase[friction_atoms[Int(ceil(fy/3,digits=0))]].mass)
        end
    end
    friction .= austrip.(friction .* model.friction_unit)
    friction .*= mass_weights
    friction .= austrip.(friction *. u"u")
end
