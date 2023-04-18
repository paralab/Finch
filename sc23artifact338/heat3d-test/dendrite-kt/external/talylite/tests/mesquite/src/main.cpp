/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //

#include <talyfem/talyfem.h>
#include <talyfem/grid/mesquite.h>

int main(int argc, char** args) {
  PetscInitialize(&argc, &args, nullptr, nullptr);

  InputData input_data;
  if (!input_data.ReadFromFile())
    throw TALYException() << "Error reading config file";
  if (!input_data.CheckInputData())
    throw TALYException() << "Invalid config file";

  std::shared_ptr<GRID> grid;
  GridField<NODEData> field;
  CreateGrid(grid, &input_data);
  field.redimNodeData(grid.get());

  // set boundaries
  /*const double eps = 1e-6;
  for (LocalNodeID i = 0; i < grid->n_nodes(); i++) {
    const ZEROPTV& p = grid->GetNode(i)->location();
    if (fabs(p.x() - 5.0) < eps || fabs(p.x() + 5.0) < eps
        || fabs(p.y() - 5.0) < eps || fabs(p.y() + 5.0) < eps) {
      grid->GetNode(i)->addIndicatorNum(1);
    }
  }*/

  namespace msq = Mesquite2;

  TalyMesqMesh mesh(grid.get());
  //msq::PlanarDomain domain(msq::PlanarDomain::XY);
  msq::PlanarDomain domain(msq::Vector3D(0, 0, -1), msq::Vector3D(0, 0, -5));

  msq::MeshDomainAssoc mesh_and_domain(&mesh, &domain);

  msq::MsqError err;
  /*msq::LaplaceWrapper mesh_quality_algorithm;
  mesh_quality_algorithm.run_instructions(&mesh_and_domain, err);*/

  /*msq::ShapeImprovementWrapper mesh_quality_algo;
  mesh_quality_algo.run_instructions(&mesh_and_domain, err);*/

  msq::IdealWeightInverseMeanRatio inverse_mean_ratio(err);
  msq::LPtoPTemplate obj_func(&inverse_mean_ratio, 2.0, err);
  msq::TrustRegion t_region(&obj_func);
  t_region.use_global_patch();

  msq::TerminationCriterion tc_inner;
  tc_inner.add_absolute_gradient_L2_norm(1e-4);
  t_region.set_inner_termination_criterion(&tc_inner);

  msq::QualityAssessor m_ratio_qa(&inverse_mean_ratio);
  msq::InstructionQueue queue;
  queue.add_quality_assessor(&m_ratio_qa, err);
  queue.set_master_quality_improver(&t_region, err);
  queue.add_quality_assessor(&m_ratio_qa, err);

  queue.run_instructions(&mesh_and_domain, err);

  if (err) {
    throw TALYException() << "Mesquite error: " << err.error_message();
  }

  save_gf(&field, &input_data, "output.plt", "out");

  PetscFinalize();
  return 0;
}

