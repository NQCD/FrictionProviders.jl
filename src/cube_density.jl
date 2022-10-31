"""
This uses a cube file to attach friction coefficients to existing models by fitting the data
provided by Gerrits et al. in PHYSICAL REVIEW B 102, 155130 (2020).
"""

using UnitfulAtomic: austrip
using NQCBase: PeriodicCell, apply_cell_boundaries!

struct CubeDensity{T}
    "Cube file reader."
    cube::Cube{T}
    "Temporary array for storing the electron density."
    cell::PeriodicCell{T}
end

function CubeDensity(filename, cell; cell_matching_rtol=1e-3)
    cube = Cube(filename)
    if !isapprox(cell.vectors, cube.cell.vectors, rtol=cell_matching_rtol)
        error("the cube file cell vectors do not match the simulation cell vectors.\n",
            "  Simulation vectors: ", cell.vectors, "\n",
            "  Cube vectors: ", cube.cell.vectors
        )
    end

    CubeDensity(cube, cell)
end

function density!(model::CubeDensity, rho::AbstractVector, R::AbstractMatrix, friction_atoms::AbstractVector)
    for i in friction_atoms
        r = R[:,i]
        apply_cell_boundaries!(model.cube.cell, r)
        rho[i] = austrip(model.cube(r + model.cube.origin) * u"Å^-3")
    end
end
