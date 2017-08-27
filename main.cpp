#include <Omega_h_dolfin.hpp>
#include <Omega_h_library.hpp>
#include <Omega_h_mesh.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_class.hpp>
#include <Omega_h_adapt.hpp>
#include "Poisson.h"

using namespace dolfin;

class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = sin(5*x[0]);
  }
};

class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

int main(int argc, char** argv)
{
  auto lib_osh = Omega_h::Library(&argc, &argv);

  std::shared_ptr<Mesh> mesh;

  mesh = std::make_shared<UnitSquareMesh>(32, 32);
  auto V = std::make_shared<Poisson::FunctionSpace>(mesh);

  auto u0 = std::make_shared<Constant>(0.0);
  auto boundary = std::make_shared<DirichletBoundary>();
  DirichletBC bc(V, u0, boundary);

  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  auto f = std::make_shared<Source>();
  auto g = std::make_shared<dUdN>();
  L.f = f;
  L.g = g;

  Function u(V);
  solve(a == L, u, bc);

  File file("poisson.pvd");
  file << u;

  Omega_h::Mesh mesh_osh(&lib_osh);
  Omega_h::from_dolfin(&mesh_osh, *mesh);
  Omega_h::from_dolfin(&mesh_osh, u, "u");
  Omega_h::classify_by_angles(&mesh_osh, Omega_h::PI / 4.0);
  Omega_h::MetricInput metric_input;
  auto source = Omega_h::MetricSource(OMEGA_H_VARIATION, 1e-3, "u");
  metric_input.sources.push_back(source);
  metric_input.should_limit_lengths = true;
  metric_input.min_length = 1.0 / 32.0;
  metric_input.max_length = 1.0 / 2.0;
  metric_input.should_limit_gradation = true;
  Omega_h::generate_target_metric_tag(&mesh_osh, metric_input);
  Omega_h::add_implied_metric_tag(&mesh_osh);
  Omega_h::vtk::write_vtu("before.vtu", &mesh_osh);
  Omega_h::AdaptOpts opts(&mesh_osh);
  opts.xfer_opts.type_map["u"] = OMEGA_H_LINEAR_INTERP;
  while (Omega_h::approach_metric(&mesh_osh, opts)) {
    Omega_h::adapt(&mesh_osh, opts);
  }
  Omega_h::vtk::write_vtu("after.vtu", &mesh_osh);

  // Plot solution
//plot(u);
//interactive();

  return 0;
}
