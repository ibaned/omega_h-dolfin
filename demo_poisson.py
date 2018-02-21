from dolfin import *
import PyOmega_h as omega_h

mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Lagrange", 1)

def boundary(x):
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS

u0 = Constant(0.0)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)

mesh_osh = omega_h.new_empty_mesh()
omega_h.mesh_from_dolfin(mesh_osh, mesh)
omega_h.classify_by_angles(mesh_osh, 3.14159 / 4.0)
mesh_osh.set_parting(omega_h.GHOSTED, 1)
omega_h.add_implied_metric_tag(mesh_osh);
mesh_osh.set_parting(omega_h.ELEM_BASED, 0);

i = 0
n = 3
while (True):
  omega_h.mesh_to_dolfin(mesh, mesh_osh);
  V = FunctionSpace(mesh, "Lagrange", 1)

  bc = DirichletBC(V, u0, boundary)

  u = TrialFunction(V)
  v = TestFunction(V)
  a = inner(grad(u), grad(v))*dx
  L = f*v*dx + g*v*ds

  u = Function(V)
  solve(a == L, u, bc)

  omega_h.function_from_dolfin(mesh_osh, u._cpp_object, "u")

  file = File("poisson.pvd")
  file << u

  i = i + 1
  if (i == n):
    break

  mesh_osh.set_parting(omega_h.GHOSTED, 1)

  source = omega_h.MetricSource(omega_h.VARIATION, 2e-3, "u")

  metric_input = omega_h.MetricInput()
  metric_input.add_source(source)
  metric_input.should_limit_lengths = True
  metric_input.max_length = 1.0 / 2.0
  metric_input.should_limit_gradation = True

  omega_h.generate_target_metric_tag(mesh_osh, metric_input)

  opts = omega_h.AdaptOpts(mesh_osh)
  opts.verbosity = omega_h.EXTRA_STATS

  while (omega_h.approach_metric(mesh_osh, opts)):
    omega_h.adapt(mesh_osh, opts)
