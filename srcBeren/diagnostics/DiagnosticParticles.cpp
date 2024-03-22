#include "Diagnostic.h"

void Writer::write_particles2D(int timestep) {
    for (auto& sp : _species) {
        sp.get_Pr();
        _mesh.make_periodic_border_with_add(sp.Pxx);
        _mesh.make_periodic_border_with_add(sp.Pyy);
        _mesh.make_periodic_border_with_add(sp.Pzz);

        std::string fnameDense = ".//Particles//" + sp.name + "//Diag2D//Dens";
        std::string fnamePxx = ".//Particles//" + sp.name + "//Diag2D//Pxx";
        std::string fnamePyy = ".//Particles//" + sp.name + "//Diag2D//Pyy";
        std::string fnamePzz = ".//Particles//" + sp.name + "//Diag2D//Pzz";
        std::string fnameJ = ".//Particles//" + sp.name + "//Diag2D//Current";
        for (auto coordX : diagData.params.sliceFieldsPlaneX) {
            write_array2D_planeX(sp.densityOnGrid, coordX, fnameDense,
                                 timestep);
            write_array2D_planeX(sp.Pxx, coordX, fnamePxx, timestep);
            write_array2D_planeX(sp.Pyy, coordX, fnamePyy, timestep);
            write_array2D_planeX(sp.Pzz, coordX, fnamePzz, timestep);
            write_fields2D_planeX(sp.currentOnGrid, coordX, fnameJ, timestep);
        }
        for (auto coordY : diagData.params.sliceFieldsPlaneY) {
            write_array2D_planeY(sp.densityOnGrid, coordY, fnameDense,
                                 timestep);
            write_array2D_planeY(sp.Pxx, coordY, fnamePxx, timestep);
            write_array2D_planeY(sp.Pyy, coordY, fnamePyy, timestep);
            write_array2D_planeY(sp.Pzz, coordY, fnamePzz, timestep);
            write_fields2D_planeY(sp.currentOnGrid, coordY, fnameJ, timestep);
        }
        for (auto coordZ : diagData.params.sliceFieldsPlaneZ) {
            write_array2D_planeZ(sp.densityOnGrid, coordZ, fnameDense,
                                 timestep);
            write_array2D_planeZ(sp.Pxx, coordZ, fnamePxx, timestep);
            write_array2D_planeZ(sp.Pyy, coordZ, fnamePyy, timestep);
            write_array2D_planeZ(sp.Pzz, coordZ, fnamePzz, timestep);
            write_fields2D_planeZ(sp.currentOnGrid, coordZ, fnameJ, timestep);
        }
        write_array2D_planeZ_avg(sp.densityOnGrid, fnameDense, timestep);
        write_array2D_planeZ_avg(sp.Pxx, fnamePxx, timestep);
        write_array2D_planeZ_avg(sp.Pyy, fnamePyy, timestep);
        write_array2D_planeZ_avg(sp.Pzz, fnamePzz, timestep);
        write_fields2D_AvgPlaneZ(sp.currentOnGrid, fnameJ, timestep);
    }
}

void Writer::write_array2D_planeX(const Array3D<double>& data, double coordX,
                                  const std::string& fname,
                                  const int& timestep) {
    int globIndex = _mesh.get_node_from_coordX(coordX);
    int index = _world.region.get_index_loc(globIndex);

    char filenameCh[100];
    float info;
    int indx;

    int size_y = data.size().y();
    int size_z = data.size().z();
    int size1 = size_y;
    int size2 = size_z;

    float* floatData = new float[size1 * size2];

    sprintf(filenameCh, (fname + "PlaneX_%04d_%04d").c_str(), globIndex,
            timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fdata2D(filename, std::ios::out | std::ios::binary);

    int i = index;
    for (auto j = 0; j < size_y; j++) {
        for (auto k = 0; k < size_z; k++) {
            indx = j * size_z + k;
            floatData[indx] = float(data(i, j, k));
        }
    }

    info = float(size1);
    fdata2D.write((char*) &info, sizeof(info));
    info = float(size2);
    fdata2D.write((char*) &info, sizeof(info));

    fdata2D.write((char*) floatData, size1 * size2 * sizeof(float));

    delete[] floatData;
}

void Writer::write_array2D_planeZ(const Array3D<double>& data, double coordZ,
                                  const std::string& fname,
                                  const int& timestep) {
    char filenameCh[100];
    float info;

    int indx;

    int size_x = data.size().x();
    int size_y = data.size().y();
    int size1 = size_x;
    int size2 = size_y;

    float* floatData = new float[size1 * size2];
    int globIndex = _mesh.get_node_from_coordZ(coordZ);

    sprintf(filenameCh, (fname + "PlaneZ_%04d_%04d").c_str(), globIndex,
            timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fdata2D(filename, std::ios::out | std::ios::binary);

    int k = _mesh.get_node_from_coordZ(coordZ);
    for (auto i = 0; i < size_x; i++) {
        for (auto j = 0; j < size_y; j++) {
            indx = i * size_y + j;
            floatData[indx] = float(data(i, j, k));
        }
    }

    info = float(size1);
    fdata2D.write((char*) &info, sizeof(info));
    info = float(size2);
    fdata2D.write((char*) &info, sizeof(info));

    fdata2D.write((char*) floatData, size1 * size2 * sizeof(float));

    delete[] floatData;
}

void Writer::write_array2D_planeZ_avg(const Array3D<double>& data,
                                      const std::string& fname,
                                      const int& timestep) {
    char filenameCh[100];
    float info;

    int indx;

    int size_x = data.size().x();
    int size_y = data.size().y();
    int size_z = data.size().z();
    int size1 = size_x;
    int size2 = size_y;

    float* floatData = new float[size1 * size2];

    sprintf(filenameCh, (fname + "PlaneAvgZ_%04d").c_str(),
            timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fdata2D(filename, std::ios::out | std::ios::binary);

    for (auto i = 0; i < size_x; i++) {
        for (auto j = 0; j < size_y; j++) {
            indx = i * size_y + j;
            floatData[indx] = 0.;
        }
    }

    for (auto i = 0; i < size_x; i++) {
        for (auto j = 0; j < size_y; j++) {
            for (auto k = 1; k < size_z - 2; k++) {
                indx = i * size_y + j;
                floatData[indx] += float(data(i, j, k) / (size_z - ADD_NODES));
            }
        }
    }

    info = float(size1);
    fdata2D.write((char*) &info, sizeof(info));
    info = float(size2);
    fdata2D.write((char*) &info, sizeof(info));

    fdata2D.write((char*) floatData, size1 * size2 * sizeof(float));

    delete[] floatData;
}

void Writer::write_array2D_planeY(const Array3D<double>& data, double coordY,
                                  const std::string& fname,
                                  const int& timestep) {
    char filenameCh[100];
    float info;

    int indx;

    int size_x = data.size().x();
    int size_z = data.size().z();
    int size1 = size_x;
    int size2 = size_z;

    float* floatData = new float[size1 * size2];
    int globIndex = _mesh.get_node_from_coordY(coordY);

    sprintf(filenameCh, (fname + "PlaneY_%04d_%04d").c_str(), globIndex,
            timestep / TimeStepDelayDiag2D);
    std::string filename(filenameCh);
    std::ofstream fdata2D(filename, std::ios::out | std::ios::binary);

    int j = _mesh.get_node_from_coordY(coordY);
    for (auto i = 0; i < size_x; i++) {
        for (auto k = 0; k < size_z; k++) {
            indx = i * size_z + k;
            floatData[indx] = float(data(i, j, k));
        }
    }

    info = float(size1);
    fdata2D.write((char*) &info, sizeof(info));
    info = float(size2);
    fdata2D.write((char*) &info, sizeof(info));

    fdata2D.write((char*) floatData, size1 * size2 * sizeof(float));

    delete[] floatData;
}
