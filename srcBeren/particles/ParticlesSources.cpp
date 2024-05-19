#include "Particles.h"

void ParticlesArray::add_uniform_cilinder(int numParts, double r0, double z0, double3 c){
	Particle particle;
	auto sigma = temperature;
	double3 coord, pulse;
	double rx,ry,x,y,z;

	for(int k = 0; k < numParts; k++){
		
		do {
		rx = (1 - 2*Uniform01() ) * r0;
		ry = (1 - 2*Uniform01() ) * r0;

		} while ( rx*rx + ry*ry > r0*r0);

		

		x = c.x() + rx;
		y = c.y() + ry;
		z = c.z() + z0*Uniform01(); 
		particle.coord = double3(x,y,z);

		pulse.x() = Gauss(sigma / sqrt(_mass));
		pulse.y() = Gauss(sigma / sqrt(_mass));
		pulse.z() = Gauss(sigma / sqrt(_mass));
		particle.velocity = pulse;

		add_particle(particle);
	}
    update_count_in_cell();
}

void ParticlesArray::add_uniform_line(int numParts, double3 r, double3 sizeL){
	Particle particle;
	auto sigma = temperature;
	double3 coord, pulse;
	double x, y, z;

	for(int k = 0; k < numParts; k++){

		x = r.x() + sizeL.x()*Uniform01();
		y = r.y() + sizeL.y()*Uniform01(); 
		z = r.z() + sizeL.z()*Uniform01(); 
		particle.coord = double3(x,y,z);
		do {
			pulse.x() = Gauss(sigma / sqrt(_mass));
			pulse.y() = Gauss(sigma / sqrt(_mass));
			pulse.z() = Gauss(sigma / sqrt(_mass));
		} while (fabs(pulse.x()) > 3*sigma || fabs(pulse.y()) > 3*sigma  || fabs(pulse.z()) > 3*sigma);

		particle.velocity = pulse;

		add_particle(particle);
	}
    update_count_in_cell();
}

/*
void ParticlesArray::source_uniform_from_bound(int timestep){
    double energyInj = 0;
	Particle particle;
	int i;//, kmax;
	double gama;
	double pz,x,y,z,px,py;
	double vb = velocity;
	int numInj = NumPartPerCell * (widthY / Dy) * (widthZ / Dz);
	double  particleMPW;
	bool LeaftInj =  (vb>0. && _world.region.boundType[0].x() == OPEN);
	bool RightInj =  (vb<0. && _world.region.boundType[1].y() == OPEN);

	if(!( LeaftInj || RightInj )) return;

	//kmax = size();

	for(i = 0; i < numInj; ++i){

		x = Dx*Uniform01();
		y =  widthY*(Uniform01() - 0.5);
		z =  widthZ*(Uniform01() - 0.5);

		if (y*y/(0.25*widthY*widthY) + z*z/(0.25*widthZ*widthZ) > 1. ) continue;

		y += 0.5 * Dy*_world.regionGlob.numCells.y() ;
		z += 0.5 * Dz*_world.regionGlob.numCells.z() ;
				
		px = vb / sqrt(1.0-vb*vb);

		py = Gauss(temperature);
		pz = Gauss(temperature);
		
		bool CondInj = x <= Dt*fabs(px / sqrt(1. + px*px + py*py + pz * pz) );
		
		if( !CondInj  ) continue;

				particle.coord.y() = y ;
				particle.coord.z() = z ;
				particle.pulse = double3(px,py,pz);
				
				if(LeaftInj ){
					  particle.coord.x() = x + Dx * _world.region.dampCells[0].x();
				}
				else if (RightInj ){
					  particle.coord.x() = Dx * _world.region.numCells.x() - x - Dx *  _world.region.dampCells[1].x();
					}
					else{
					  std::cout << "Error Inject Particles type" << std::endl; exit(0); 
				}
				
			
				#ifdef PARTICLE_MASS
		                   particle.mass = _mass;
				#endif	
				#ifdef PARTICLE_MPW
		                   particle.mpw = particleMPW = std::min(1.,Dt*timestep/InjectSmoothTime) * density * Dx * Dy * Dz / (NumPartPerCell);
;
				#else
						particleMPW = _mpw;
				#endif	
				gama = sqrt(_mass*_mass + dot(particle.pulse,particle.pulse) );
					
				energyInj += particleMPW*(gama - _mass);
				
				add_particle_scatter(particle);			
		
		}
			
	injectionEnergy += energyInj;
}

void ParticlesArray::source_focused_gauss(int timestep){
    double energyInj = 0;
    Particle particle;
	int i;//, kmax;
	double gama;
	double pz,x,y,z,px,py,px0,py0,pz0;
	double vb = velocity;
	double particleMPW;
	double sigma;
	double alpha = option.sourceAngle;

	int numInj = 0.5*PI* NumPartPerCell*widthY/Dy*widthZ/Dz;

	bool LeaftInj =  (vb>0. &&  _world.region.boundType[0].x() == OPEN);
	bool RightInj =  (vb<0.  && _world.region.boundType[1].x() == OPEN );

	if(!( LeaftInj || RightInj )) return;
	
	for(i = 0; i < numInj; ++i){

		x = Dx * Uniform01();

		y = 0.5*Dy*(_world.region.numCells.y()) + Gauss(0.5*widthY);
		z = 0.5*Dz*(_world.region.numCells.z()) + Gauss(0.5*widthZ);

		px0 = vb / sqrt(1.0-vb*vb);

		sigma = (Rmax-0.5*widthY)/Lmax*px0;
		
		py0 = Gauss(sigma/3.);
		pz0 = Gauss(sigma/3.);
		
		px = px0*cos(alpha) + pz0*sin(alpha); 
		pz = pz0*cos(alpha) - px0*sin(alpha);
		py = py0; 
//		r1 = sqrt((x + focus * px/pz)*(x + focus * px/pz)+(y + focus * py/pz)*(y + focus * py/pz));
		
		//x = x - focus * px/px;
		y = y - focus * py/px;
		z = z - focus * pz/px;
		
		//bool InFocus = rf < 3.*width;

		double radius_y = 0.9*Dy * (0.5*_world.region.numCells.y() - _world.region.dampCells[0].y());
		double radius_z = 0.9*Dz * (0.5*_world.region.numCells.z() - _world.region.dampCells[0].z());

		double a1 = ((y-0.5*Dy*(_world.region.numCells.y()) ) / radius_y);
		double b1 = ((z-0.5*Dz*(_world.region.numCells.z()) ) / radius_z);

		bool CondCurr = x <= Dt*fabs(px/sqrt(1. + px*px + py*py + pz * pz));
		
		bool CondInj = (a1*a1+b1*b1 <=1.) && CondCurr; //&& InFocus;
		//std::cout<<x << " " << y << " " << z << " " <<px << " "<< py << " "<< pz <<"\n";
		if( !CondInj ) continue;

				particle.coord.y() = y ;
				particle.coord.z() = z ;
					
				if(LeaftInj){
					  particle.coord.x() =  x + Dx *  _world.region.dampCells[0].x();
				}
				else if ( RightInj){
					  particle.coord.x() = Dx * _world.region.numCells.x() - x - Dx *  _world.region.dampCells[1].x();
					  py = - py;
					  pz = - pz;
					}
					else{
					  std::cout << "Error Inject Particles type" << std::endl; 
					  exit(0); 
				}
				particle.pulse = double3(px,py,pz);

				#ifdef PARTICLE_MASS
		                   particle.mass = _mass;
				#endif	
				#ifdef PARTICLE_MPW
		                   particle.mpw = particleMPW = std::min(1.,Dt*timestep/InjectSmoothTime) * density * Dx * Dy * Dz / (NumPartPerCell);
;
				#else
						particleMPW = _mpw;
				#endif	
				gama = sqrt(1. + dot(particle.pulse,particle.pulse) );
					
				energyInj += particleMPW*(gama - 1.);
				
				add_particle_scatter(particle);		
			
	}
			
	injectionEnergy += energyInj;
	if (size() > MaxSizeOfParts) {
		std::cout << "particle array is overflowing by injection!\n";
		exit(1);
	}
}

void injectionFocused(ParticlesArray& particles ,  const World& world,\
            int timestep){
                double energyInj = 0;
        
	int i, kmax;
	int count;
	double pb,gama;
	double r,r1,r2,rf,z,pr, pp, pz,x,y,px,py;
	double vb = particles.velocity;
	double width = particles.width;
	double dens = particles.density;
	double zFocus = particles.focus;
	double sigma;
	int numInj = 4. / PI * NumPartPerCell * Rfocus / Dr;
	bool LeaftInj =  (particles.name == "BeamLeft" && world.MPIconf.RankLine == 0);
	bool RightInj =  (particles.name == "BeamRight" && world.MPIconf.RankLine == world.MPIconf.SizeLine -1 );

	if(!( LeaftInj || RightInj )) return;

	count = -1;
	kmax = particles.size();

	for(i = 0; i < numInj; ++i){

		z = Dz * Uniform01();

		x = width * Uniform01();
		y = width * Uniform01();
		
		rf = sqrt(x*x+y*y);
		
		pb = vb / sqrt(1.0-vb*vb);

		sigma = (Rmax - width) / Lmax * pb;
		
		px = Gauss(sigma/3.);
		py = Gauss(sigma/3.);

		pz = pb; 
		
		r1 = sqrt((x + zFocus * px/pz)*(x + zFocus * px/pz)+(y + zFocus * py/pz)*(y + zFocus * py/pz));

		x = x - zFocus * px/pz;
		y = y - zFocus * py/pz; 
		
		r2 = sqrt(x*x+y*y);
		r = r2;
		pr = x/r * px + y/r*py; 
		pp = - y/r * px + x/r*py; 
		
		bool InFocus = rf < width;
		bool Bound1 = r1 < 0.8*Dr*(NumCellsR_glob - DampCellsR_glob);
		bool Bound2 = r2 < 0.8*Dr*(NumCellsR_glob - DampCellsR_glob);
		bool CondCurr = z <= Dt*fabs(pz/sqrt(1. + pr*pr + pp*pp + pz * pz));
		
		bool CondInj =  InFocus && Bound1 && Bound2 && CondCurr ;
		
		if( CondInj ) {

			++count;
			
			if( count % world.MPIconf.SizeDepth == world.MPIconf.RankDepth ){


				particles(kmax).r = r ;
				particles(kmax).pz = pz;
					
				if(LeaftInj){
					  particles(kmax).z = z + Dz *  world.region.dampCells_d1[0];
					  particles(kmax).pp = pp;
					  particles(kmax).pr = pr;
				}
				else if ( RightInj){
					  particles(kmax).z = Dz * world.region.numCells_d1 - z - Dz *  world.region.dampCells_d1[1];
					  particles(kmax).pp = - pp;
					  particles(kmax).pr = - pr;
					}
					else{
					  std::cout << "Error Inject Particles type" << std::endl; 
					  exit(0); 
				}
				
				particles(kmax).density = std::min(1.,Dt*timestep/InjectSmoothTime)*dens * PI * Dr * Dz * width / (NumPartPerCell);
#if MASS_OWN == 1
		                   particles(kmax).M = particles.mass;
#endif	
				gama = sqrt(1. + pr*pr + pp*pp + pz*pz);
				energyInj += particles(kmax).density*(gama - 1.);
				++kmax;	
			}
		}
	}
			
	particles.particlesData.size_ = kmax;
	particles.injectionEnergy += energyInj;
}

*/
