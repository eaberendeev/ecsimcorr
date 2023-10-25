#include "Diagnostic.h"
#include <assert.h>

void DiagData::calc_energy(const Mesh &mesh, const std::vector<ParticlesArray> &species){
  double3 velocity;
  double3 v0(0.2,0.,0.);
  for ( auto &sp : species){
    energyParticlesKinetic[sp.name] = sp.get_kinetic_energy();
    energyParticlesInject[sp.name] = sp.lostEnergy;
        // for (int ip = 0; ip < sp.size(); ip++ ) {
        //   velocity = sp.particlesData(ip).velocity;
        // }
  }
  //velocity = 0.5*(velocity + v0);

  auto Jfull = Dt*mesh.fieldJp.data() + mesh.Lmat2*(mesh.fieldE.data() + mesh.fieldEp.data());
  auto Jfull2 = Dt*mesh.fieldJe.data();

  std::cout << "En-Ep norm " << (mesh.fieldEn.data() - mesh.fieldEp.data()).norm() << "\n";
  std::cout << "Je-Jp norm " << (Jfull - Jfull2).norm() << "\n";



  energyFieldE = mesh.calc_energy_field(mesh.fieldEn);
  energyFieldB = mesh.calc_energy_field(mesh.fieldB);
  double energyFieldE0 = mesh.calc_energy_field(mesh.fieldE);
  double energyFieldB0 = mesh.calc_energy_field(mesh.fieldB0);

  Je = 0.5*Dt*mesh.fieldJe.data().dot(mesh.fieldE.data() + mesh.fieldEn.data());
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
    ss << "Area_E^2 " << "Area_B^2 " << "diffV " << "JE " << "diffEB ";
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
  std::vector<double> vecEnergy = { diagData.energyFieldE, diagData.energyFieldB};

  for(uint i = 0; i< vecEnergy.size(); ++i){
  ss << vecEnergy[i] << " "; 
  }
  ss << diffV << " "; 
  ss << diagData.Je << " "; 
  ss << diagData.diffEB  << " diff Energy " << fabs(diffV + diagData.diffEB) << " "; 

  ss <<"\n";
      fprintf(fDiagEnergies, "%s",  ( ss.str() ).c_str() ); 
      std::cout << ss.str();
    if( timestep % TimeStepDelayDiag1D == 0){
      fflush(fDiagEnergies);
    }
}
