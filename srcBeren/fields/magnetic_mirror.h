#ifndef MAGNETIC_MIRROR_HPP
#define MAGNETIC_MIRROR_HPP

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <vector>
#include <iostream>
//#include "util.h"

using Eigen::Vector3d;
using Eigen::Vector3i;
typedef Eigen::VectorXd Field;

//! @brief This class describes two dimensional magnetic mirrors.
class Magnetic_mirrors {
protected:
	// Axis enumarator to provide a notation. 
	enum axis { x = 0, y, z };
	
	Vector3d global_size;
	Vector3d local_size;
	Vector3d offset;

	// Storage of magnetic mirror parameters.
	struct mirror_parameters {
		
		/**
		 * @brief Construct a new mirror parameters object.
		 * 
		 * @param z_mirror 	z-coordinate of the mirror.
		 * @param coil_radius Radius of the current coil that provides magnetic field.
		 * @param reduced_current Combination of constants with current to compute a magnetic field.
		 */
		mirror_parameters(double z_mirror, double coil_radius, double reduced_current)
		:	z_mirror(z_mirror), coil_radius(coil_radius), reduced_current(reduced_current) {};
		
		const double z_mirror;
		const double coil_radius;
		const double reduced_current;
	};
  
	std::vector<mirror_parameters> parameters;

public:
	// Default destructor.
	virtual ~Magnetic_mirrors() = default;

	/**
	 * @brief Construct a new Magnetic_mirrors object.
	 * 
	 * @param global_size Global size of the computational domain [in meters].
	 * @param local_size  Size of a subdomain, including ghost cells [in meters]. 
	 * @param offset 	  Radus vector to the beginning of a subdomain [in meters]. 
	 */
	Magnetic_mirrors(const Vector3d& global_size, const Vector3d& local_size,const Vector3d& offset);

	/**
	 * @brief Adds a mirror, described by its parameters, to the list.
	 * 
	 * @param B_mirror 	Magnetic field at the mirror's center [in Tesla].
	 * @param z_mirror 	z-coordinate of the mirror [in meters].
	 * @param B_some 	Magnetic field on z-axis, at some distance from the mirror [in Tesla].
	 * @param z_some 	z-coordinate of the point where \a B_some is set [in meters].
	 */
	virtual void add_mirror(double B_mirror, double z_mirror, double B_some, double z_some);

	/**
	 * Gets a field from a magnetic system at the point \a r. 
	 * 
	 * @brief Passes through every stored mirror and calls
	 * 		the field_from_mirror, then sums up the results.
	 * 
	 * @param r Radius vector to the point where magnetic field is needed [in meters].
	 * @return Magnetic field at the point \a r [in Tesla].
	 */
	virtual Vector3d get_field(const Vector3d& r) const;

private:
	/**
	 * Gets a magnetic field from a single mirror at the point \a r.
	 * 
	 * @brief Considering the magnetic mirror as two one-dimensional
	 * 		current, it calculates the field at the desired point. 
	 * 
	 * @param r Radius vector to the point where magnetic field is needed [in meters].
	 * @param concrete Preset description of the magnetic mirror.
	 * 
	 * @return Magnetic field at the point \a r [in Tesla].
	 */
	Vector3d field_from_mirror(const Vector3d& r, const mirror_parameters& concrete) const;
};


//! @brief This class provides a three dimensional magnetic mirrors.
class Magnetic_mirrors3D : public Magnetic_mirrors {
public:
	/**
	 * @brief Construct a new Magnetic_mirrors3D object.
	 * 
	 * @param global_size Global size of the computational domain [in meters].
	 * @param local_size  Size of a subdomain, including ghost cells [in meters]. 
	 * @param offset 	  Radus vector to the beginning of a subdomain [in meters].
	 * @param dx 		  Width of the virtual grid for storing a field components [in meters].
	 * @param integrarion_steps Number of integration steps to accurately compute the integrals.
	 */
	Magnetic_mirrors3D(
	const Vector3d& global_size, const Vector3d& local_size, const Vector3d& offset,
	double dx, int integration_steps);

	/**
	 * @brief Adds a mirror to the list. It precomputes field
	 *		components Br and Bz in every cell of the grid for later use.
	 * 
	 * @param B_mirror 	   Magnetic field at the mirror [in Tesla].
	 * @param z_mirror 	   z-coordinate of the mirror [in meters].
	 * @param B_some 	   Magnetic field on the z-axis, at some distance from the mirror [in Tesla].
	 * @param z_some 	   z-coordinate of the point where the \a B_some is set [in meters].
	 */
	void add_mirror(double B_mirror, double z_mirror, double B_some, double z_some) override;

	/**
	 * Gets a field from a magnetic system at the point \a r.
	 * 
	 * @brief Interpolates magnetic field from virtual grid to the position of \a r.
	 * 
	 * @param r Radius vector to the point where magnetic field is needed [in meters].
	 * 
	 * @return Magnetic field at the point \a r [in Tesla].
	 */
	Vector3d get_field(const Vector3d& r) const override;

private: 
	static Field field;

	const int integration_steps;
	const double dphi;

	const double dx; 
	const int nx_limit, ny_limit, nz_limit;

	inline constexpr int vind(int x, int y, int z, int c) const {
		return c + 2 * (z + nz_limit * (y + ny_limit * x));
	}
};


#endif // MAGNETIC_MIRROR_HPP
