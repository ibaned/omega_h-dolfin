from dolfin import *
import PyOmega_h as omega_h
import numpy as np

# TODO assert is missing
def test_functionspace():
    comm_osh = omega_h.world()
    mesh_osh = omega_h.build_box(comm_osh, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 32, 32, 0)
    mesh_osh.set_parting(omega_h.Parting.GHOSTED, 1)
    omega_h.add_implied_metric_tag(mesh_osh)
    mesh_osh.set_parting(omega_h.Parting.ELEM_BASED, 0)

    mesh = Mesh()

    opts = omega_h.AdaptOpts(mesh_osh)
    omega_h.adapt(mesh_osh, opts)

    omega_h.mesh_to_dolfin(mesh, mesh_osh)

    V = FunctionSpace(mesh, "Lagrange", 1)


def test_mesh():
    mesh = UnitSquareMesh(32, 32)
    mesh_2 = UnitSquareMesh(12, 12)

    assert mesh.hash() != mesh_2.hash()

    mesh_osh = omega_h.new_empty_mesh()

    # TODO Does not work
    #omega_h.mesh_to_dolfin(mesh, mesh_osh)
    omega_h.mesh_from_dolfin_unit_square(mesh_osh, mesh)

    omega_h.mesh_to_dolfin(mesh_2, mesh_osh);

    assert mesh.hash() == mesh_2.hash()
    np.testing.assert_array_equal(mesh.cells(), mesh_2.cells())


# TODO assert is missing
def test_function():
    mesh = UnitSquareMesh(32, 32)
    mesh_osh = omega_h.new_empty_mesh()

    omega_h.mesh_from_dolfin_unit_square(mesh_osh, mesh)

    V = FunctionSpace(mesh, "Lagrange", 1)
    u = Function(V)

    omega_h.function_from_dolfin(mesh_osh, u._cpp_object, "u")

#def test_adapt_new_functionspace():
#    mesh = UnitSquareMesh(32, 32)
#    mesh_osh = omega_h.new_empty_mesh()
#
#    omega_h.mesh_from_dolfin(mesh_osh, mesh)
#
#    omega_h.add_implied_metric_tag(mesh_osh);
#
#    opts = omega_h.AdaptOpts(mesh_osh)
#    omega_h.adapt(mesh_osh, opts)
#
#    omega_h.mesh_to_dolfin(mesh, mesh_osh);
#
#    V = FunctionSpace(mesh, "Lagrange", 1)

