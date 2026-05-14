const int NumProcs = 1;
const int NumAreas = 1;
const double Dx = 0.5;
const double Dy = 0.5;
const double Dz = 0.5;
const double Dt = 1.5;
const double Tau = 6000;
const int NumCellsX_glob = 300;
const int NumCellsY_glob = 300;
const int NumCellsZ_glob = 5;
const int NumCellsXmax_glob = 303;
const int NumCellsYmax_glob = 303;
const int DampCellsX_glob[2] = {20, 20};
const int DampCellsY_glob[2] = {20, 20};
const int DampCellsZ_glob[2] = {0, 0};
const int PlasmaCellsY_glob = 1;
const int PlasmaCellsX_glob = 300;
const int PlasmaCellsZ_glob = 1;
const int NumOfPartSpecies = 2;
const int NumPartPerLine = 1;
const int NumPartPerCell = 2000;
const int MaxTimeStep = 88;
const int RecTimeStep = 400;
const int StartTimeStep = 0;
const int TimeStepDelayDiag1D = 1;
const int TimeStepDelayDiag2D = 20;
const double BUniform[3] = {0, 0, 0.2};
const double BCoil[76] = {25,    -13.75, 70,    3.0,  -12.5, 70,    -3.0,  -11.25, 70,    3.0,   -10.0, 70,   -3.0,
                          -8.75, 70,     3.0,   -7.5, 70,    -3.0,  -6.25, 70,     3.0,   -5.0,  70,    -3.0, -3.75,
                          70,    3.0,    -2.5,  70,   -3.0,  -1.25, 70,    3.0,    0.0,   70,    -3.0,  1.25, 70,
                          3,     2.5,    70,    -3,   3.75,  70,    3,     5.0,    70,    -3,    6.25,  70,   3,
                          7.5,   70,     -3,    8.75, 70,    3,     10.0,  70,     -3,    11.25, 70,    3,    12.5,
                          70,    -3,     13.75, 70,   3,     15.0,  70,    -3,     16.25, 70,    3};
const int BoundTypeY_glob[2] = {2, 2};
const int BoundTypeZ_glob[2] = {2, 2};
const int BoundTypeX_glob[2] = {2, 2};
const double MC2 = 512.0;
const double n0 = 10000000000000.0;
const int isVelDiag = 0;
const int PxMax = 800;
const int PpMax = 800;
const int MaxSizeOfParts = 1200001.0;
const double PI = 3.141592653589793;
const double k_particles_reservation = -1.0;
