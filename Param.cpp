/* ***************************************************************************** 
Class for input parameters
***************************************************************************** */

#include "Param.h"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;

double PI_val=3.1415926535897932384626433832795;

void Param::loadCellData()
{
    boost::property_tree::ptree pt;
    boost::property_tree::ini_parser::read_ini("Input Data/cell_data.txt", pt);

  // [coupling] --- files for data transfer (not used now)
  fileConcCells = pt.get<string>("coupling.fileConcCells","nil");
  fileCells2FEM = pt.get<string>("coupling.fileCells2FEM","cells.txt");
  fileCellsDensity2FEM = pt.get<string>("coupling.fileCellsDensity2FEM","trias.txt");
  fileFEM2Cells = pt.get<string>("coupling.fileFEM2Cells","concentration_O2.txt");

  // [cells]
  n_initial_cells  = pt.get<int>("cells.n_initial_cells",1);
  n_steps = pt.get<int>("cells.n_steps",1);
  time_step = pt.get<double>("cells.time_step",1.);
  time_death = pt.get<double>("cells.time_death",2880.);
  growth_rate = to_array<double>(pt.get<std::string>("cells.growth_rate","0.1,0.1"));
  alpha_birthrate = to_array<double>(pt.get<std::string>("cells.alpha_birthrate","1e-3,1e-3"));
  alpha_YoungM = to_array<double>(pt.get<std::string>("cells.alpha_YoungM","1e-3,1e-3"));
  alpha_PoissonNo = to_array<double>(pt.get<std::string>("cells.alpha_PoissonNo","0.5,0.5"));
  Gcm = pt.get<double>("cells.Gcm",1e-2);
  alpha_gcm = to_array<double>(pt.get<std::string>("cells.alpha_gcm","1.0,1.0"));
  adhesion_value = to_array<double>(pt.get<std::string>("cells.adhesion_value"));
  variance_motion = pt.get<double>("cells.variance_motion",4e-4);
  variance_phenotype = pt.get<double>("cells.variance_phenotype",0.);
  variance_adhesion = pt.get<double>("cells.variance_adhesion",0.);
  string pheno = "0";
  for (int i=0; i<(n_initial_cells-1); i++) pheno.append(",0");
  ic_phenotype = to_array<int>(pt.get<std::string>("cells.ic_phenotype",pheno));
  string fl = "1";
  for (int i=0; i<(n_initial_cells-1); i++) fl.append(",1");
  ic_follower_leader= to_array<int>(pt.get<std::string>("cells.ic_follower_leader",fl));
  string pol = "0.0";
  for (int i=0; i<(n_initial_cells-1); i++) pol.append(",0.0");
  polarity_x = to_array<double>(pt.get<std::string>("cells.polarity_x",pol));
  polarity_y = to_array<double>(pt.get<std::string>("cells.polarity_y",pol));
  polarity_z = to_array<double>(pt.get<std::string>("cells.polarity_z",pol));
  follower_force = pt.get<double>("cells.follower_force",0.);
  follower_denominator = pt.get<double>("cells.follower_denominator",0.);
  be_multiplier = pt.get<double>("cells.be_multiplier",12.0);
  be_displacement = pt.get<double>("cells.be_displacement",1.5);

  // [oxygen]
  oxygen_response = pt.get<double>("oxygen.oxygen_response",0.);
  threshold_death = pt.get<double>("oxygen.threshold_death",0.7);
  threshold_hypo = pt.get<float>("oxygen.threshold_hypo",7.);

  // [fem] --- pde (oxygen)
  femSolverType = pt.get<int>("fem.femSolverType",0);
  meshdir = pt.get<string>("fem.meshdir","./");
  meshname = pt.get<string>("fem.meshname","rectangle_7x4x1.5_4Knodes.mesh");
  FreeFemCall = pt.get<string>("fem.FreeFemCall","FreeFem++");
  FreeFemFile = pt.get<string>("fem.FreeFemFile","diffusion3d.edp");

  // [geo] --- boxes and geometry
  dimension = pt.get<int>("geo.dimension",3);
  boxesx = pt.get<int>("geo.boxesx",61);
  boxesy = pt.get<int>("geo.boxesy",61);
  boxesz = pt.get<int>("geo.boxesz",61);
  lattice_length_x = pt.get<double>("geo.lattice_length_x",600.1);
  lattice_length_y = pt.get<double>("geo.lattice_length_y",600.1);
  lattice_length_z = pt.get<double>("geo.lattice_length_z",600.1);
  string pos = "300.0";
  for (int i=0; i<(n_initial_cells-1); i++) pos.append(",300.0");
  ic_cell_x = to_array<double>(pt.get<std::string>("geo.ic_cell_x",pos));
  ic_cell_y = to_array<double>(pt.get<std::string>("geo.ic_cell_y",pos));
  ic_cell_z = to_array<double>(pt.get<std::string>("geo.ic_cell_z",pos));
  max_cell = pt.get<double>("geo.max_cell", 100000.);

  // [fibres] --- fibre parameters
  n_sub_domain = 2;
  x_start.resize(n_sub_domain);
  x_end.resize(n_sub_domain);
  x_start[0] = 0.0;
  x_start[1] = lattice_length_x/2.0;
  x_end[0] = lattice_length_x/2.0;
  x_end[1] = lattice_length_x;
  n_initial_fibres.resize(n_sub_domain);
  fibre_orientation_distribution.resize(n_sub_domain);
  fibre_orientation_mean_phi.resize(n_sub_domain);
  fibre_orientation_variance_phi.resize(n_sub_domain);
  fibre_orientation_mean_theta.resize(2);
  fibre_orientation_variance_theta.resize(2);
  fibre_length_mean.resize(2);
  fibre_length_variance.resize(2);
  for (unsigned int i=0; i<n_sub_domain; i++) {
    n_initial_fibres[i] =  pt.get<int>("fibres.n_initial_fibres",0);
    fibre_orientation_distribution[i] =
      pt.get<int>("fibres.fibre_orientation_distribution",0);
    fibre_orientation_mean_phi[i] =
      pt.get<double>("fibres.fibre_orientation_mean_phi",PI_val);
    fibre_orientation_variance_phi[i] =
      pt.get<double>("fibres.fibre_orientation_variance_phi",1.);
    fibre_orientation_mean_theta[i] =
      pt.get<double>("fibres.fibre_orientation_mean_theta",PI_val);
    fibre_orientation_variance_theta[i] =
      pt.get<double>("fibres.fibre_orientation_variance_theta",1.);
    if (dimension==2) {
      fibre_orientation_mean_theta[i] = PI_val/2.;
      fibre_orientation_variance_theta[i] = 0.;
    }
    fibre_length_mean[i] =  pt.get<double>("fibres.fibre_length_mean",50.);
    fibre_length_variance[i] =  pt.get<double>("fibres.fibre_length_variance",2.);
    }
    fibre_radius = pt.get<double>("fibres.fibre_radius",5.);
    fibre_make_gap = pt.get<int>("fibres.fibre_make_gap",0);
    vel_adhesion = pt.get<double>("fibres.vel_adhesion",0.);
    vel_contact = pt.get<double>("fibres.vel_contact",0.);

  // [vessels] --- fibre parameters
  vessel_PoissonNo = pt.get<double>("vessels.vessel_PoissonNo",0.5);
  vessel_YoungM = pt.get<double>("vessels.vessel_YoungM",1e-3);
  vessel_adhesion = pt.get<double>("vessels.vessel_adhesion",3.72e-4);

  // [postprocessing]
  verbose = pt.get<int>("postprocessing.verbose",0);
  writeVtkFiles = pt.get<int>("postprocessing.writeVtkFiles",0);
  outputDirectory = pt.get<string>("postprocessing.outputDirectory","./vtk/");
  testcase = pt.get<string>("postprocessing.testcase","test");
  fileCellsVisualization =
    pt.get<string>("postprocessing.fileCellsVisualization","celulas.txt");
  cellTracking = pt.get<int>("postprocessing.cellTracking",0);
  fileCellsTracking = pt.get<string>("postprocessing.fileCellsTracking","track.txt");
  casename = pt.get<string>("postprocessing.casename","case");
  casedirectory = pt.get<string>("postprocessing.casedirectory","./cell_output/");

  // print parameter database
  print();

}

void Param::print()
{
  cout << endl;
  cout << " # ======================= " << endl;
  cout << " # CELL DATA " << filename << endl;
  cout << " # ======================= " << endl;
  cout << endl;

  cout << "[coupling]" << endl;
  cout << "fileConcCells = " << fileConcCells << endl;
  cout << "fileCells2FEM = " << fileCells2FEM << endl;
  cout << "fileCellsDensity2FEM = " << fileCellsDensity2FEM << endl;
  cout << "fileFEM2Cells = " << fileFEM2Cells << endl;
  cout << endl;

  cout << "[fem]" << endl;
  cout << "meshdir = " << meshdir << endl;
  cout << "meshname = " << meshname << endl;
  cout << "FreeFemFile = " << FreeFemFile << endl;
  cout << endl;

  cout << "[cells]" << endl;
  cout << "n_initial_cells = " << n_initial_cells << endl;
  cout << "n_steps = " << n_steps << endl;
  cout << "time_step = " << time_step << endl;
  cout << "time_death = " << time_death << endl;
  cout << "growth_rate = [" << growth_rate[0] << "," << growth_rate[1] << "]" << endl;
  cout << "alpha_birthrate = [" << alpha_birthrate[0] << "," << alpha_birthrate[1] << "]" << endl;
  cout << "alpha_YoungM = [" << alpha_YoungM[0] << "," << alpha_YoungM[1] << "]" << endl;
  cout << "alpha_PoissonNo = [" << alpha_PoissonNo[0] << "," << alpha_PoissonNo[1] << "]" << endl;
  cout << "Gcm = " << Gcm << endl;
  cout << "alpha_gcm = [" << alpha_gcm[0] << "," << alpha_gcm[1] << "]" << endl;
  cout << "adhesion_value = [" << adhesion_value[0] << "," << adhesion_value[1] << "]" << endl;
  cout << "variance_motion = " << variance_motion << endl;
  cout << "variance_phenotype = " << variance_phenotype << endl;
  cout << "variance_adhesion = " << variance_adhesion << endl;
  cout << "ic_phenotype = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << ic_phenotype[i] << ",";
  }
  cout << ic_phenotype[n_initial_cells-1] << "]" << endl;
  cout << "polarity_x = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << polarity_x[i] << ",";
  }
  cout << polarity_x[n_initial_cells-1] << "]" << endl;
  cout << "polarity_y = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << polarity_y[i] << ",";
  }
  cout << polarity_y[n_initial_cells-1] << "]" << endl;
  cout << "polarity_z = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << polarity_z[i] << ",";
  }
  cout << polarity_z[n_initial_cells-1] << "]" << endl;
  cout << "ic_follower_leader = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << ic_follower_leader[i] << ",";
  }
  cout << ic_follower_leader[n_initial_cells-1] << "]" << endl;
  cout << "follower_force = " << follower_force << endl;
  cout << "follower_denominator = " << follower_denominator << endl;
  cout << "be_multiplier = " << be_multiplier << endl;
  cout << "be_displacement = " << be_displacement << endl;
  cout << endl;

  cout << "[oxygen]" << endl;
  cout << "oxygen_response = " << oxygen_response << endl;
  cout << "threshold_death = " << threshold_death << endl;
  cout << "threshold_hypo = " << threshold_hypo << endl;
  cout << endl;

  cout << "[fibres]" << endl;
  cout << "n_initial_fibres = [" << n_initial_fibres[0] << "," << n_initial_fibres[1] << "]" << endl;
  cout << "fibre_orientation_distribution = [" << fibre_orientation_distribution[0] << "," << fibre_orientation_distribution[1] << "]" << endl;
  cout << "fibre_orientation_mean_phi = [" << fibre_orientation_mean_phi[0] << "," << fibre_orientation_mean_phi[1] << "]" << endl;
  cout << "fibre_orientation_variance_phi = [" << fibre_orientation_variance_phi[0] << "," << fibre_orientation_variance_phi[1] << "]" << endl;
  cout << "fibre_orientation_mean_theta = [" << fibre_orientation_mean_theta[0] << "," << fibre_orientation_mean_theta[1] << "]" << endl;
  cout << "fibre_orientation_variance_theta = [" << fibre_orientation_variance_theta[0] << "," << fibre_orientation_variance_theta[1]  << "]" << endl;
  cout << "fibre_length_mean = [" << fibre_length_mean[0] << "," << fibre_length_mean[1] << "]" << endl;
  cout << "fibre_length_variance = [" << fibre_length_variance[0] << "," << fibre_length_variance[1] << "]" << endl;
  cout << "fibre_make_gap = " << fibre_make_gap << endl;
  cout << "vel_adhesion = " << vel_adhesion << endl;
  cout << "vel_contact = " << vel_contact << endl;
  cout << endl;

  cout << "[vessels]" << endl;
  cout << "vessel_PoissonNo = " << vessel_PoissonNo << endl;
  cout << "vessel_YoungM = " << vessel_YoungM << endl;
  cout << "vessel_adhesion = " << vessel_adhesion << endl;
  cout << endl;


  cout << "[geo]" << endl;
  cout << "dimension = " << dimension << endl;
  cout << "boxesx = " << boxesx << endl;
  cout << "boxesy = " << boxesy << endl;
  cout << "boxesz = " << boxesz << endl;
  cout << "lattice_length_x = " << lattice_length_x << endl;
  cout << "lattice_length_y = " << lattice_length_y << endl;
  cout << "lattice_length_z = " << lattice_length_z << endl;
  cout << "ic_cell_x = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << ic_cell_x[i] << ",";
  }
  cout << ic_cell_x[n_initial_cells-1] << "]" << endl;
  cout << "ic_cell_y = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << ic_cell_y[i] << ",";
  }
  cout << ic_cell_y[n_initial_cells-1] << "]" << endl;
  cout << "ic_cell_z = [";
  for (int i=0; i<(n_initial_cells-1); i++) {
    cout << ic_cell_z[i] << ",";
  }
  cout << ic_cell_z[n_initial_cells-1] << "]" << endl;
  cout << "max_cell = " << max_cell << endl;
  cout << endl;
}






