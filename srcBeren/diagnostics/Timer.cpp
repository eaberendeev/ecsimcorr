#include "Timer.h"

void Timer::write(int timestep){
  if( times.empty() || timestep%TimeStepDelayDiag1D !=0) return;
    std::stringstream ss;

    if( firstWrite){
      ss<< "Time ";
        for (auto it = times.begin(); it != times.end(); ++it){
          ss << it->first << " ";
        }
        firstWrite = false;
         
    }
    ss << "\n" << timestep << " ";

    for (auto it = times.begin(); it != times.end(); ++it){
          ss << it->second << " ";
    } 
    
    fprintf(fTimes, "%s",  ( ss.str() ).c_str() ); 
    
    reset();

  if( timestep % TimeStepDelayDiag1D == 0 ) {
    fflush(fTimes);
  }
  
}