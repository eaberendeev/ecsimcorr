#include "magnetic_mirror.h"


Magnetic_mirrors::Magnetic_mirrors(
const Vector3d& global_size, const Vector3d& local_size,const Vector3d& offset)
:	global_size(global_size), local_size(local_size), offset(offset) {};


void Magnetic_mirrors::add_mirror(
double B_mirror, double z_mirror, double B_some, double z_some)
{
	double coil_radius = fabs(z_some - z_mirror) / sqrt((B_mirror * B_mirror) / (B_some * B_some) - 1.);
	double reduced_current = B_mirror * coil_radius * coil_radius;

    parameters.emplace_back(z_mirror, coil_radius, reduced_current);
}


Vector3d Magnetic_mirrors::field_from_mirror(
const Vector3d& r, const mirror_parameters& concrete) const
{
    const double x = r[axis::x] - concrete.z_mirror;
	const double y = r[axis::y] - 0.5 * global_size[axis::y];

	const double cr_py = concrete.coil_radius + y; 
	const double cr_my = concrete.coil_radius - y;

	const double denominator = (x * x + cr_py * cr_py) * (x * x + cr_my * cr_my);

	const double Bx = concrete.reduced_current * (concrete.coil_radius * concrete.coil_radius + x * x - y * y) / denominator;
	const double By = concrete.reduced_current * (2. * x * y ) / denominator;
	const double Bz = 0.;

	return Vector3d(Bx, By, Bz);
}


Vector3d Magnetic_mirrors::get_field(const Vector3d& r) const
{
	Vector3d B_local =  {0., 0., 0.};
	Vector3d r_global = r + offset;
	
    for (const auto& concrete : parameters) {
		B_local += field_from_mirror(r_global, concrete);
	}

	return B_local;
}


Field Magnetic_mirrors3D::field = Field::Zero(1);


Magnetic_mirrors3D::Magnetic_mirrors3D(
const Vector3d& global_size, const Vector3d& local_size, const Vector3d& offset,
double dx, int integration_steps)
:	Magnetic_mirrors(global_size, local_size, offset),
	integration_steps(integration_steps), dphi(2. * M_PI / (integration_steps - 1.)),
	dx(dx),
	nx_limit(static_cast<int>(round(local_size[axis::x] / dx)+1)), 
	ny_limit(static_cast<int>(round(local_size[axis::y] / dx)+1)),
	nz_limit(static_cast<int>(round(local_size[axis::z] / dx)+1))
	{
		Magnetic_mirrors3D::field.resize( nx_limit * ny_limit * nz_limit * 2);
		Magnetic_mirrors3D::field.setZero();
		std::cout << "Mesh for magnetic mirrors " << nx_limit << "x" << ny_limit << "x" << nz_limit <<std::endl;  
	};


void Magnetic_mirrors3D::add_mirror(
double B_mirror, double z_mirror, double B_some, double z_some)
{
	double coil_radius = fabs(z_some - z_mirror) / sqrt(pow(B_mirror / B_some, 2. / 3.) - 1.);
	double reduced_current = B_mirror * coil_radius * coil_radius / (2. * M_PI);
	parameters.emplace_back(z_mirror, coil_radius, reduced_current);

	// Presetting integrands for the Br and Bz integrals.
	auto denominator = [&coil_radius](double z, double r, double phi) { 
		return pow(z * z + coil_radius * coil_radius + r * r - 2. * coil_radius * r * cos(phi), 3. / 2.);
	};

	auto integrand_Br = [&denominator](double z, double r, double phi) {
		return z * cos(phi) / denominator(z, r, phi);
	};

	auto integrand_Bz = [&denominator, &coil_radius](double z, double r, double phi) {
		return (coil_radius - r * cos(phi)) / denominator(z, r, phi);
	};

	// Pass through every cell in the subdomain and precompute a magnetic field.
	for(int nx = 0; nx < nx_limit; ++nx) {
	for(int ny = 0; ny < ny_limit; ++ny) {
	for(int nz = 0; nz < nz_limit; ++nz) {
		
		double aa = ny * dx + offset[axis::y] - 0.5 * global_size[axis::y];
		double bb = nz * dx + offset[axis::z] - 0.5 * global_size[axis::z];
		double z  = nx * dx + offset[axis::x] - z_mirror;
		double r = sqrt(aa * aa + bb * bb);

		double Br = 0.;
		double Bz = 0.;
	
		// Calculate integrals using Simpson method.
		for(int i = 0; i < integration_steps; ++i)
		{
			// Cylindric coordinates (x, y, z) -> (x, rho, phi)	
			Br +=   integrand_Br(z, r, i * dphi)
			 + 4. * integrand_Br(z, r, (i + 0.5) * dphi)
				  + integrand_Br(z, r, (i + 1.) * dphi);
	
			Bz +=   integrand_Bz(z, r, i * dphi)
			 + 4. * integrand_Bz(z, r, (i + 0.5) * dphi)
				  + integrand_Bz(z, r, (i + 1.) * dphi);
		}

		field[vind(nx, ny, nz, 0)] += Br * reduced_current * dphi / 6.;
		field[vind(nx, ny, nz, 1)] += Bz * reduced_current * dphi / 6.;
	}}}
}

Vector3d Magnetic_mirrors3D::get_field(const Vector3d& r) const
{
	double x = r[axis::x];
	double y = r[axis::y]; // we will move on to the cylindrical coordinates at the end. 
	double z = r[axis::z]; // we will move on to the cylindrical coordinates at the end. 

	const int l = static_cast<int>(floor(x / dx));
	const int m = static_cast<int>(floor(y / dx));
	const int n = static_cast<int>(floor(z / dx));
	
	double formx[2];
	double formy[2];
	double formz[2];

	formx[0] = x / dx - l;
	formy[0] = y / dx - m;
	formz[0] = z / dx - n;
	formx[1] = 1. - formx[0];
	formy[1] = 1. - formy[0];
	formz[1] = 1. - formz[0];
					
	double Bz = 0.;
	double Br = 0.;
	
	// Using linear interpolation to get the field from the virtual grid
	Br = formx[1] * (formy[1] * (formz[1] * field[vind(l,   m,   n, 0)] + formz[0] * field[vind(l,   m,   n+1, 0)]) 
				  +  formy[0] * (formz[1] * field[vind(l,   m+1, n, 0)] + formz[0] * field[vind(l,   m+1, n+1, 0)]))   
	   + formx[0] * (formy[1] * (formz[1] * field[vind(l+1, m,   n, 0)] + formz[0] * field[vind(l+1, m,   n+1, 0)]) 
				  +  formy[0] * (formz[1] * field[vind(l+1, m+1, n, 0)] + formz[0] * field[vind(l+1, m+1, n+1, 0)]));

	Bz = formx[1] * (formy[1] * (formz[1] * field[vind(l,   m,   n, 1)] + formz[0] * field[vind(l,   m,   n+1, 1)])
				  +  formy[0] * (formz[1] * field[vind(l,   m+1, n, 1)] + formz[0] * field[vind(l,   m+1, n+1, 1)])) 
	   + formx[0] * (formy[1] * (formz[1] * field[vind(l+1, m,   n, 1)] + formz[0] * field[vind(l+1, m,   n+1, 1)]) 
				  +  formy[0] * (formz[1] * field[vind(l+1, m+1, n, 1)] + formz[0] * field[vind(l+1, m+1, n+1, 1)]));


	// Cylindric coordinates (x, y, z) -> (x, rho, phi)
	y += offset[axis::y] - 0.5 * global_size[axis::y];
	z += offset[axis::z] - 0.5 * global_size[axis::z];

	const double rho = sqrt( y * y + z * z );

	return Vector3d({ Bz, Br * y / rho, Br * z / rho });
}