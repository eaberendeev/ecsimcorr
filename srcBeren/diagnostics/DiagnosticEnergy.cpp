#include "Diagnostic.h"
#include <assert.h>

void DiagData::calc_energy(Mesh &mesh, const std::vector<ParticlesArray> &species){
  double3 velocity;
  for ( auto &sp : species){
    energyParticlesKinetic[sp.name] = sp.get_kinetic_energy();
    energyParticlesInject[sp.name] = sp.injectionEnergy;
    energyParticlesLost[sp.name] = sp.lostEnergy;
    energy[sp.name + "Z"] = sp.get_kinetic_energy(Axis::Z);
    energy[sp.name + "XY"] = sp.get_kinetic_energy(Axis::X,Axis::Y);
  }

  auto Jfull = Dt*mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data());
  auto Jfull2 = Dt*mesh.fieldJe.data();

  std::cout << "En-Ep norm " << (mesh.fieldEn.data() - mesh.fieldEp.data()).norm() << "\n";
  std::cout << "Je-Jp norm " << (Jfull - Jfull2).norm() << "\n";

  mesh.fieldB -= mesh.fieldBInit;
  mesh.fieldB0 -= mesh.fieldBInit;

  energyFieldE = mesh.calc_energy_field(mesh.fieldEn);
  energyFieldB = mesh.calc_energy_field(mesh.fieldB);
  double energyFieldE0 = mesh.calc_energy_field(mesh.fieldE);
  double energyFieldB0 = mesh.calc_energy_field(mesh.fieldB0);
  mesh.fieldB += mesh.fieldBInit;
  mesh.fieldB0 += mesh.fieldBInit;
  energyFieldBFull = mesh.calc_energy_field(mesh.fieldB);

  diffEB = energyFieldB + energyFieldE - energyFieldB0 - energyFieldE0;
}

void Writer::write_energies(double diffV, int timestep){
  std::stringstream ss;

  if(timestep == 0){
    ss << "Time ";
    for (auto it = diagData.energyParticlesKinetic.begin(); it != diagData.energyParticlesKinetic.end(); ++it){
      ss << "Area_" << it->first << " ";
    }
    for (auto it = diagData.energyParticlesInject.begin(); it != diagData.energyParticlesInject.end(); ++it){
      ss << "Injection_" << it->first << " ";
    }
    for (auto it = diagData.energy.begin();
         it != diagData.energy.end(); ++it) {
        ss << "Energy_" << it->first << " ";
    }
    ss << "Area_E^2 "
       << "Area_B^2 "
       << "Area_B_Full^2 "
       << "diffV "
       << "diffEB "
       << "energyError";
    ss << "\n";
  }
  
  ss << timestep*Dt << " ";
    
  for (auto it = diagData.energyParticlesKinetic.begin(); it != diagData.energyParticlesKinetic.end(); ++it){
    double energyParticles = it->second;
    ss << energyParticles << " ";
  }
  
  for (auto it = diagData.energyParticlesInject.begin(); it != diagData.energyParticlesInject.end(); ++it){
    double energyInject = it->second;
    ss << energyInject << " ";
  }
  for (auto it = diagData.energy.begin();
       it != diagData.energy.end(); ++it) {
      double energy = it->second;
      ss << energy << " ";
  }
  std::vector<double> vecEnergy = {
      diagData.energyFieldE,     diagData.energyFieldB,
      diagData.energyFieldBFull, diffV,
      diagData.diffEB,           fabs(diffV + diagData.diffEB)};

  for(uint i = 0; i< vecEnergy.size(); ++i){
  ss << vecEnergy[i] << " "; 
  }

  ss <<"\n";
      fprintf(fDiagEnergies, "%s",  ( ss.str() ).c_str() ); 
      std::cout << ss.str();
    if( timestep % TimeStepDelayDiag1D == 0){
      fflush(fDiagEnergies);
    }
}
