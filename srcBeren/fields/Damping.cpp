#include "Damping.h"

#include <algorithm>

#include "World.h"

void Damping_Func(double& source, double i, double maxi, double& energyDamp, double coeff) {
    double a, damp;
    if (i >= maxi)
        return;

    a = (1.0 - coeff) / (maxi * maxi);
    damp = a * i * i + coeff;

    energyDamp += 0.5 * source * source * (1.0 - damp * damp);
    source *= damp;
}

double damping_fields(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const nlohmann::json& config) {
    std::string damping_type = config.value("DampingType", "None");
    if (damping_type == "None")
        return 0.0;

    const double coeff = config.value("DampingCoeff", 0.8);
    double energyDamp = 0.0;

    if (damping_type == "CircleXY") {
        energyDamp = damping_fields_circleXY(fieldE, fieldB, domain, config, coeff);
    } else if (damping_type == "Rectangle") {
        energyDamp = damping_fields_rectangle(fieldE, fieldB, config, coeff);
    } else {
        std::cout << "DampingType is not defined" << std::endl;
        exit(1);
    }
    return energyDamp;
}

double damping_fields_circleXY(Field3d& fieldE, Field3d& fieldB, const Domain& domain, const json& config,
                               double coeff) {
    const auto sizes = fieldE.sizes();
    int i, j, k;
    double energyDamp;
    const double dx = domain.cell_size().x();
    const auto dampX = config.value("DampCellsX_glob", json::array({0, 0}));
    const int dampLeft = (dampX.is_array() && dampX.size() > 0) ? dampX[0].get<int>() : 0;
    const int dampRight = (dampX.is_array() && dampX.size() > 1) ? dampX[1].get<int>() : 0;
    const double dampSize = std::max(dampLeft, dampRight) * dx;
    double dampRadius = 0.5 * dx * (sizes.x() - 1);
    int dampRadiusInd = (sizes.x() - 1) / 2;

    energyDamp = 0.;

    for (i = 0; i < sizes.x(); i++) {
        for (j = 0; j < sizes.y(); j++) {
            for (k = 0; k < sizes.z(); k++) {
                double r = sqrt(
                    dx * dx * ((i - dampRadiusInd) * (i - dampRadiusInd) + (j - dampRadiusInd) * (j - dampRadiusInd)));
                if (r <= dampRadius - dampSize) {
                    continue;
                } else if (r >= dampRadius) {
                    for (int dim = 0; dim < 3; dim++) {
                        fieldE(i, j, k, dim) *= coeff;
                        fieldB(i, j, k, dim) *= coeff;
                    }
                } else {
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), dampRadius - r, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), dampRadius - r, dampSize, energyDamp, coeff);
                    }
                }
            }
        }
    }
    return energyDamp;
}

double damping_fields_rectangle(Field3d& fieldE, Field3d& fieldB, const json& config, double coeff) {
    double energyDamp;
    const auto sizes = fieldE.sizes();

    int max_indx = sizes.x();
    int max_indy = sizes.y();
    int max_indz = sizes.z();

    const int dampSizeXLeft = config.at("DampCellsX_glob")[0];
    const int dampSizeYLeft = config.at("DampCellsY_glob")[0];
    const int dampSizeZLeft = config.at("DampCellsZ_glob")[0];
    const int dampSizeXRight = config.at("DampCellsX_glob")[1];
    const int dampSizeYRight = config.at("DampCellsY_glob")[1];
    const int dampSizeZRight = config.at("DampCellsZ_glob")[1];

    energyDamp = 0.;

    if (dampSizeXLeft > 0) {
        const int dampSize = dampSizeXLeft;
        for (int i = 0; i < dampSize; i++) {
            for (int j = 0; j < max_indy; j++) {
                for (int k = 0; k < max_indz; k++) {
                    const int currentIndex = i;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                    }
                }
            }
        }
    }
    if (dampSizeXRight > 0) {
        const int dampSize = dampSizeXRight;

        // currentIndex — расстояние от правой кромки слоя в ячейках
        for (int i = max_indx - dampSize; i < max_indx; i++) {
            for (int j = 0; j < max_indy; j++) {
                for (int k = 0; k < max_indz; k++) {
                    const int currentIndex = max_indx - 1 - i;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
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
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                    }
                }
            }
        }
    }
    if (dampSizeYRight > 0) {
        const int dampSize = dampSizeYRight;
        for (int i = 0; i < max_indx; ++i) {
            // currentIndex — расстояние от правой кромки слоя в ячейках
            for (int j = max_indy - dampSize; j < max_indy; ++j) {
                for (int k = 0; k < max_indz; ++k) {
                    const int currentIndex = max_indy - 1 - j;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
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
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                    }
                }
            }
        }
    }
    if (dampSizeZRight > 0) {
        const int dampSize = dampSizeZRight;

        for (int i = 0; i < max_indx; ++i) {
            for (int j = 0; j < max_indy; ++j) {
            // currentIndex — расстояние от правой кромки слоя в ячейках
            for (int k = max_indz - dampSize; k < max_indz; ++k) {
                const int currentIndex = max_indz - 1 - k;
                    for (int dim = 0; dim < 3; dim++) {
                        Damping_Func(fieldE(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                        Damping_Func(fieldB(i, j, k, dim), currentIndex, dampSize, energyDamp, coeff);
                    }
                }
            }
        }
    }

    return energyDamp;
}
