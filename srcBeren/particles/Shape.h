#ifndef SHAPE_H_
#define SHAPE_H_

inline double Shape(const double& dist){
	double d = fabs(dist);
	if ( d < 1. ) return (1. - d);
	else 
	 return 0.;
}

inline double Shape2(const double& dist){
	double d = fabs(dist);
	if ( d <= 0.5 ) return (0.75 - d * d);
	else 
	    if (d < 1.5)  return ((d - 1.5) * (d-1.5) * 0.5);
	    else return 0.;
}

inline double Shape3(const double& dist){
	double d = fabs(dist);
	if(d<=1.) return (2./3. - d * d + 0.5 * d * d * d);
	else 
	  if ( d < 2 ) return ((2. - d) * (2. - d) * (2. - d) / 6.);
	  else return 0.;
}

inline double Shape4(const double& dist){
	double d = fabs(dist);
	if(d<=0.5) return ( 115. / 192 - 0.625 * d * d + 0.25 * d * d * d * d);
	else
	    if (d<=1.5 ) return (55. + 20. * d - 120. * d*d + 80.* d * d * d - 16. * d * d * d * d )/96.;
	    else 
	      if( d<2.5 )  return (5.-2.*d)*(5.-2.*d)*(5.-2.*d)*(5.-2.*d)/384.;
	      else return 0.;
} 

#endif // SHAPE_H_