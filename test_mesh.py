from dolfin import *
import PyOmega_h as omega_h

def test_mesh():
    mesh = UnitSquareMesh(32, 32)
    mesh_2 = UnitSquareMesh(12, 12)

    assert mesh.hash() != mesh_2.hash()

    mesh_osh = omega_h.new_empty_mesh()

    omega_h.mesh_from_dolfin_unit_square(mesh_osh, mesh)

    omega_h.mesh_to_dolfin(mesh_2, mesh_osh);

    assert mesh.hash() == mesh_2.hash()
