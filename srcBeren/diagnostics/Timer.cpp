#include "Timer.h"

void Timer::write(int timestep, int delay){
  if( times.empty() || timestep%delay !=0) return;
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

  if( timestep % delay == 0 ) {
    fflush(fTimes);
  }
  
}