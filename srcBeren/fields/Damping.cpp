#include "Damping.h"

#include "World.h"
void Damping_Func(double& source, double i, double maxi, double& energyDamp) {
    double koeff = 0.8;
    double a, damp;
    if (i >= maxi)
        return;

    a = (1.0 - koeff) / (maxi * maxi);
    damp = a * i * i + koeff;

    energyDamp += 0.5 * source * source * (1.0 - damp * damp);
    source *= damp;
}

double damping_fields(Field3d& fieldE, Field3d& fieldB, const Domain& domain,
                               const ParametersMap& parameters){
    if (parameters.get_string("DampingType") == "None") return 0.0;

    double energyDamp = 0.0;

        if (parameters.get_string("DampingType") == "CircleXY") {
            energyDamp =
                damping_fields_circleXY(fieldE, fieldB, domain, parameters);
        } else if (parameters.get_string("DampingType") == "Rectangle") {
            energyDamp =
                damping_fields_rectangle(fieldE, fieldB, domain, parameters);
        } else{
            std::cout << "DampingType is not defined" << std::endl;
            exit(1);
        }
        return energyDamp;
    }

double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB,
                               const Domain& domain, const ParametersMap &parameters) {
    int i, j, k;
    double energyDamp;
    const double dx = domain.cell_size().x();
    const double dampSize = parameters.get_int("DampCellsX_glob") * dx;
    double dampRadius = 0.5 * dx * (fieldE.size().x() - 1);
    int dampRadiusInd = (fieldE.size().x() - 1) / 2;

    energyDamp = 0.;

    for (i = 0; i < fieldE.size().x(); i++) {
        for (j = 0; j < fieldE.size().y(); j++) {
            for (k = 0; k < fieldE.size().z(); k++) {
                double r = sqrt(dx * dx *
                                ((i - dampRadiusInd) * (i - dampRadiusInd) +
                                 (j - dampRadiusInd) * (j - dampRadiusInd)));
                if (r <= dampRadius - dampSize) {
                    continue;
                } else if (r >= dampRadius) {
                    for (int dim = 0; dim < 3; dim++) {
                        fieldE(i, j, k, dim) *= 0.8;
                        fieldB(i, j, k, dim) *= 0.8;
                    }
                } else {
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), dampRadius - r,
                                     dampSize, energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), dampRadius - r,
                                     dampSize, energyDamp);
                    }
                }
            }
        }
    }
	return energyDamp;
}

double damping_fields_rectangle(Field3d& fieldE, Field3d& fieldB, const Domain& domain,
                      const ParametersMap& parameters) {
    double energyDamp;
    int max_indx = fieldE.size().x();
    int max_indy = fieldE.size().y();
    int max_indz = fieldE.size().z();

    const int dampSizeXLeft = parameters.get_int("DampCellsX_glob");
    const int dampSizeYLeft = parameters.get_int("DampCellsY_glob");
    const int dampSizeZLeft = parameters.get_int("DampCellsZ_glob");
    const int dampSizeXRight = parameters.get_int("DampCellsX_glob",1);
    const int dampSizeYRight = parameters.get_int("DampCellsY_glob",1);
    const int dampSizeZRight = parameters.get_int("DampCellsZ_glob",1);

    energyDamp = 0.;

    if (dampSizeXLeft > 0) {
        const int dampSize = dampSizeXLeft;
        for (int i = 0; i < dampSize; i++) {
            for (int j = 0; j < max_indy; j++) {
                for (int k = 0; k < max_indz; k++) {
                    const int currentIndex = i;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }
    if (dampSizeXRight > 0) {
        const int dampSize = dampSizeXRight;

        for (int i = max_indx; i > max_indx - dampSize; i--) {
            for (int j = 0; j < max_indy; j++) {
                for (int k = 0; k < max_indz; k++) {
                     const int currentIndex = -(i - max_indx);
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }

    if (dampSizeYLeft > 0) {
        const int dampSize = dampSizeYLeft;
        for (int i = 0; i < max_indx; ++i) {
            for (int j = 0; j < dampSize; ++j) {
                for (int k = 0; k < max_indz; ++k) {
                     const int currentIndex = j;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }
    if (dampSizeYRight > 0) {
        const int dampSize = dampSizeYRight;
        for (int i = 0; i < max_indx; ++i) {
            for (int j = max_indy; j > max_indy - dampSize; --j) {
                for (int k = 0; k < max_indz; ++k) {
                     const int currentIndex = -j + max_indy;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }
    if (dampSizeZLeft > 0) {
        const int dampSize = dampSizeZLeft;

        for (int i = 0; i < max_indx; ++i) {
            for (int j = 0; j < max_indy; ++j) {
                for (int k = 0; k < dampSize; ++k) {
                     const int currentIndex = k;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }
    if (dampSizeZRight > 0) {
        const int dampSize = dampSizeZRight;

        for (int i = 0; i < max_indx; ++i) {
            for (int j = 0; j < max_indy; ++j) {
                for (int k = max_indz; k > max_indz - dampSize; --k) {
                     const int currentIndex = -k + max_indz;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }

    return energyDamp;
}
