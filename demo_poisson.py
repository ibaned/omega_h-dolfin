from dolfin import *
import PyOmega_h as omega_h

class MeshContainer(object):
    def __init__(self):
        comm_osh = omega_h.world()
        self.mesh_osh = omega_h.build_box(comm_osh, omega_h.Family.SIMPLEX, 1.0, 1.0, 0.0, 32, 32, 0)
        self.mesh_osh.set_parting(omega_h.Parting.GHOSTED, 1)
        omega_h.add_implied_metric_tag(self.mesh_osh)
        self.mesh_osh.set_parting(omega_h.Parting.ELEM_BASED, 0)

        self.mesh = Mesh()
        omega_h.mesh_to_dolfin(self.mesh, self.mesh_osh)


def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS


meshes = MeshContainer()
mesh = meshes.mesh
mesh_osh = meshes.mesh_osh

u0 = Constant(0.0)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)

file = File("poisson.pvd")

i = 0
n = 3

def create_weakforms(V):
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx + g*v*ds

    return a, L

def create_metric_input():
    source = omega_h.MetricSource(omega_h.VARIATION, 2e-3, "u")
    metric_input = omega_h.MetricInput()
    metric_input.add_source(source)
    metric_input.should_limit_lengths = True
    metric_input.max_length = 1.0 / 2.0
    metric_input.should_limit_gradation = True

    return metric_input

while (True):
    omega_h.mesh_to_dolfin(mesh, mesh_osh);

    V = FunctionSpace(mesh, "Lagrange", 1)
    bc = DirichletBC(V, u0, boundary)
    a, L = create_weakforms(V)

    u = Function(V)
    solve(a == L, u, bc)

    omega_h.function_from_dolfin(mesh_osh, u._cpp_object, "u")

    file << (u,i)

    i = i + 1
    if (i == n):
        break

    mesh_osh.set_parting(omega_h.GHOSTED, 1)

    metric_input = create_metric_input()
    omega_h.generate_target_metric_tag(mesh_osh, metric_input)
    opts = omega_h.AdaptOpts(mesh_osh)
    opts.verbosity = omega_h.EXTRA_STATS

    while (omega_h.approach_metric(mesh_osh, opts)):
        omega_h.adapt(mesh_osh, opts)
