#include "Damping.h"

#include "World.h"
#include "bounds.h"
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

double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB,
                               const Region& domain) {
    int i, j, k;
    double energyDamp, dampSize;
    dampSize = 10 * Dx;
    double dampRadius = 0.5 * Dx * (fieldE.size().x() - 1);
    int dampRadiusInd = (fieldE.size().x() - 1) / 2;

    energyDamp = 0.;

    for (i = 0; i < fieldE.size().x(); i++) {
        for (j = 0; j < fieldE.size().y(); j++) {
            for (k = 0; k < fieldE.size().z(); k++) {
                double r = sqrt(Dx * Dx *
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

double damping_fields(Field3d& fieldE, Field3d& fieldB, const Region& domain,
                      bool dampX, bool dampY, bool dampZ) {
    int i, j, k, i1, j1, k1;
    double energyDamp, dampSize;
    int max_indx = domain.numNodes.x() - 1;
    int max_indy = domain.numNodes.y() - 1;
    int max_indz = domain.numNodes.z() - 1;

    energyDamp = 0.;

    if (dampX) {
        for (i = 0; i < domain.dampCells[0].x(); i++) {
            for (j = 0; j <= max_indy; j++) {
                for (k = 0; k <= max_indz; k++) {
                    i1 = i;
                    dampSize = domain.dampCells[0].x();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), i1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), i1, dampSize,
                                     energyDamp);
                    }
                }
            }
        }

        for (i = max_indx; i > max_indx - domain.dampCells[1].x(); i--) {
            for (j = 0; j <= max_indy; j++) {
                for (k = 0; k <= max_indz; k++) {
                    i1 = -(i - max_indx);
                    dampSize = domain.dampCells[1].x();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), i1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), i1, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }

    if (dampY) {
        for (i = 0; i <= max_indx; ++i) {
            for (j = max_indy; j > max_indy - domain.dampCells[1].y(); --j) {
                for (k = 0; k <= max_indz; ++k) {
                    j1 = -j + max_indy;
                    dampSize = domain.dampCells[1].y();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, 0), j1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, 2), j1, dampSize,
                                     energyDamp);
                    }
                }
            }
            for (j = 0; j < domain.dampCells[0].y(); ++j) {
                for (k = 0; k <= max_indz; ++k) {
                    j1 = j;
                    dampSize = domain.dampCells[0].y();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), j1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, dim), j1, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }

    if (dampZ) {
        for (i = 0; i <= max_indx; ++i) {
            for (j = 0; j <= max_indy; ++j) {
                for (k = max_indz; k > max_indz - domain.dampCells[1].z();
                     --k) {
                    k1 = -k + max_indz;
                    dampSize = domain.dampCells[1].z();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, 0), k1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, 2), k1, dampSize,
                                     energyDamp);
                    }
                }
                for (k = 0; k < domain.dampCells[0].z(); ++k) {
                    k1 = k;
                    dampSize = domain.dampCells[0].z();
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, 0), k1, dampSize,
                                     energyDamp);
                        Damping_Func(fieldB(i, j, k, 2), k1, dampSize,
                                     energyDamp);
                    }
                }
            }
        }
    }
    return energyDamp;
}
