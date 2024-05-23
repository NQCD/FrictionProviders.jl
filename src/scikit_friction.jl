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

function SciKitFriction(descriptors, ml_model, atoms, scaler; friction_unit=u"ps^-1", descriptor_atoms=[1])
    SciKitFriction(descriptors, ml_model, atoms, scaler, friction_unit, descriptor_atoms)
end


function set_coordinates!(model::SciKitFriction, R)
    model.atoms.set_positions(ustrip.(auconvert.(u"Å", R')))
end


# function spin_singularity(z;z0=13.9,Gamma=0.01)
#     lorentzian =   (4/π)  * ((Gamma) /  ((z - z0)^2 +(Gamma)^2 ))
#     lorentzian
# end
function F_PearsonVII(x, p_H, p_x0, p_omega, p_sigma)
    f =  p_H / ((1 + ((2 * (x - p_x0) * sqrt(2^ (1 / p_omega) - 1)) /p_sigma)^ 2)^p_omega)
end

function monoExp(x, m, t, b, c)
    f = m * exp(-t * x + c) + b
end

function fitted_function(z)

    popt1 = [1.00049208,13.74298538,  2.81859272,  0.01748954]
    popt2 = [ 1.05460886, 13.7397012,   0.64277216,  0.0451869 ]
    popt3 = [1.31754445e+01, 6.08067956e+00, 0., 7.89657321e+01]

    zz = 0

    if (z>=13.71) && (z<=13.74)

        zz = F_PearsonVII(z,popt1...)

    elseif (z>13.74) && (z<=14.0)

        zz =  F_PearsonVII(z,popt2...)

    else #(z>14.0) 
        # && (z<=20)

        zz = monoExp(z,popt3...)

    # else
    #     zz=0
    end

    return zz * 156.764979 
end


function find_min_xy_dist(atoms,idx1,idx2)

    new_atoms = atoms
    new_atoms.positions[:,3] .= 0.

    distances = new_atoms.get_distances(idx1.-1,idx2.-1, mic=true)

    minimum(distances)
end

function skfriction!(model::SciKitFriction,R::AbstractMatrix, f::AbstractMatrix, friction_atoms::AbstractVector)

    #apply_cell_boundaries!(model.atoms.cell, R)
    set_coordinates!(model, R)


    r_desc = model.descriptors.create(model.atoms, positions=model.descriptor_atoms.-1, n_jobs=-1) #n_threads)
    r_desc = model.scaler.transform(r_desc)
    m_out = model.ml_model.predict(r_desc) 



    # ODF-S1 model3:  This model shifts the z coordinate of the spin peak.
    # It gave very poor results during a trajectory, so we no longer use it.


    # Find minimum xy distance to Pt atom
    #idx1 = [1]
    #idx2 = [22,23,24,25]
    #dist_to_top = find_min_xy_dist(model.atoms, idx1, idx2)
    #offset = (1.6163596822847763-dist_to_top)*0.125

    # Evaluate spin peak fit shifted due to proximity to top site
    #if (model.atoms.positions[1,3])>=(13.71)
    #    fit = fitted_function((model.atoms.positions[1,3]+offset))
    #    if fit>1.
    #        m_out[3] = fit
    #    elseif (model.atoms.positions[1,3])>=(13.93)
    #        m_out[3] = fit 
    #    end
    #end 


    # ODF-S1 model2
    # This model has the spin peak at a fixed z coordinate.
    if (model.atoms.positions[1,3])>=13.71
        m_out[3] = fitted_function((model.atoms.positions[1,3]))
    end

    m_out = austrip.(m_out.* model.friction_unit)


    if size(m_out)[1] == 1
        for i=1:size(m_out)[2]
            f[i,i] = austrip(m_out[1,i] * 1u"u") #Hardcode for H atom mass
        end
    else
        f = m_out
    end
  
    f
end

