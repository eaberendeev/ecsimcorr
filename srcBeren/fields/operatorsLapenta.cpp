// #include "Mesh.h"
// #include "World.h"
// #include "Shape.h"
// #include "SolverFDTD.h"
// #include "SolverFDTD_PML.h"
//     void triplet_dxdy(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     void triplet_dxdz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     void triplet_dydz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     void triplet_dxdx(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     void triplet_dydy(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     void triplet_dzdz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val);
//     Operator& nodes_dxdx(Operator &oper, int dimN, int dim);
//     Operator& nodes_dydy(Operator &oper, int dimN, int dim);
//     Operator& nodes_dzdz(Operator &oper, int dimN, int dim);
//     Operator& nodes_dxdz(Operator &oper, int dimN, int dim);
//     Operator& nodes_dxdy(Operator &oper, int dimN, int dim);
//     Operator& nodes_dydz(Operator &oper, int dimN, int dim);
//     Operator& nodes_dx(Operator &oper, int dimN, int dim);
//     Operator& nodes_dy(Operator &oper, int dimN, int dim);
//     Operator& nodes_dz(Operator &oper, int dimN, int dim);
//     Operator& cells_dx(Operator &oper, int dimC, int dim);
//     Operator& cells_dy(Operator &oper, int dimC, int dim);
//     Operator& cells_dz(Operator &oper, int dimC, int dim);

// // matrix ColMajor
// // (0,0) (0,1) (0,2)
// // (1,0) (1,1) (1,2)
// // (2,0) (2,1) (2,2)

// // E(i,j,k), first E[0,0,0]
// // B (i+1/2,j+1/2,k+1/2), first B[-1/2,-1/2,-1/2]

// // Derivatives in the form of triplets
// void Mesh::triplet_dxdy(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {
//   val *= 1. / (4.*Dx*Dy);
//   trips.push_back(Trip(ind_row,vind(i-1,j-1,k,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i+1,j-1,k,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i-1,j+1,k,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i+1,j+1,k,dimN),val));
// }

// void Mesh::triplet_dxdz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {
//   val *= 1. / (4.*Dx*Dz);
//   trips.push_back(Trip(ind_row,vind(i-1,j,k-1,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i+1,j,k-1,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i-1,j,k+1,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i+1,j,k+1,dimN),val));
// }

// void Mesh::triplet_dydz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {  
//   val *= 1. / (4.*Dy*Dz);
//   trips.push_back(Trip(ind_row,vind(i,j-1,k-1,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i,j+1,k-1,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i,j-1,k+1,dimN),-val));         
//   trips.push_back(Trip(ind_row,vind(i,j+1,k+1,dimN),val));
// }

// void Mesh::triplet_dxdx(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {
//   val *= 1. / (Dx*Dx);
//   trips.push_back(Trip(ind_row,vind(i-1,j,k,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i  ,j,k,dimN),-2*val));         
//   trips.push_back(Trip(ind_row,vind(i+1,j,k,dimN),val));
// }

// void Mesh::triplet_dydy(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {
//   val *= 1. / (Dy*Dy);  
//   trips.push_back(Trip(ind_row,vind(i,j-1,k,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i,j  ,k,dimN),-2*val));         
//   trips.push_back(Trip(ind_row,vind(i,j+1,k,dimN),val));
// }

// void Mesh::triplet_dzdz(std::vector<Trip> &trips, int ind_row, int i, int j, int k,
//                 int dimN, double val)
// {
//   val *= 1. / (Dz*Dz);
//   trips.push_back(Trip(ind_row,vind(i,j,k-1,dimN),val));         
//   trips.push_back(Trip(ind_row,vind(i,j,k  ,dimN),-2*val));         
//   trips.push_back(Trip(ind_row,vind(i,j,k+1,dimN),val));
// }

// // Grid derivatives at nodes
// // Makes the derivative of the oper at the grid nodes 
// // and returns the result at the cell centers

// // oper - operator on which the derivative is taken. 
// // dimN - dimension (for nodes) from which the derivative is taken.
// // dim -  dimendion in which the result will be placed
// Operator& Mesh::nodes_dxdx(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*27);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dimN);

//         val = 1./4.;
//         triplet_dxdx(trips,ind_row,i,j,k,dimN,val);
//         val = 1./8;
//         triplet_dxdx(trips,ind_row,i,j-1,k,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j,k-1,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j+1,k,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j,k+1,dimN,val);
//         val = 1./16.;
//         triplet_dxdx(trips,ind_row,i,j-1,k-1,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j+1,k-1,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j-1,k+1,dimN,val);
//         triplet_dxdx(trips,ind_row,i,j+1,k+1,dimN,val);
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dydy(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*27);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dimN);

//         val = 1./4.;
//         triplet_dydy(trips,ind_row,i,  j,  k,  dimN,val);
//         val = 1./8;
//         triplet_dydy(trips,ind_row,i-1,j,  k,  dimN,val);
//         triplet_dydy(trips,ind_row,i  ,j,  k-1,dimN,val);
//         triplet_dydy(trips,ind_row,i+1,j,  k,  dimN,val);
//         triplet_dydy(trips,ind_row,i  ,j,  k+1,dimN,val);
//         val = 1./16.;
//         triplet_dydy(trips,ind_row,i-1,j,  k-1,dimN,val);
//         triplet_dydy(trips,ind_row,i+1,j,  k-1,dimN,val);
//         triplet_dydy(trips,ind_row,i-1,j,  k+1,dimN,val);
//         triplet_dydy(trips,ind_row,i+1,j,  k+1,dimN,val);
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dzdz(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*27);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dimN);

//         val = 1./4.;
//         triplet_dzdz(trips,ind_row,i,  j,  k,  dimN,val);
//         val = 1./8;
//         triplet_dzdz(trips,ind_row,i,  j-1,k,  dimN,val);
//         triplet_dzdz(trips,ind_row,i-1,j,  k,  dimN,val);
//         triplet_dzdz(trips,ind_row,i,  j+1,k,  dimN,val);
//         triplet_dzdz(trips,ind_row,i+1,j,  k,  dimN,val);
//         val = 1./16.;
//         triplet_dzdz(trips,ind_row,i-1,j-1,k  ,dimN,val);
//         triplet_dzdz(trips,ind_row,i-1,j+1,k  ,dimN,val);
//         triplet_dzdz(trips,ind_row,i+1,j-1,k  ,dimN,val);
//         triplet_dzdz(trips,ind_row,i+1,j+1,k  ,dimN,val);
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dxdy(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*12);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dim);
//         val = 1./2.;
//         triplet_dxdy(trips,ind_row,i,  j,  k,  dimN,val);
//         val = 1./4.;
//         triplet_dxdy(trips,ind_row,i,  j,  k-1,dimN,val);
//         triplet_dxdy(trips,ind_row,i,  j,  k+1,dimN,val);   
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dxdz(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*12);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dim);
//         val = 1./2.;
//         triplet_dxdz(trips,ind_row,i,  j,  k,  dimN,val);
//         val = 1./4.;
//         triplet_dxdz(trips,ind_row,i,  j-1,k,  dimN,val);
//         triplet_dxdz(trips,ind_row,i,  j+1,k,  dimN,val);   
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dydz(Operator &oper, int dimN, int dim)
// { 
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*12);
//   double val;

//   for(int i = 1; i < size.x()-1; i++){
//     for(int j = 1; j < size.y()-1; j++){
//       for(int k = 1; k < size.z()-1; k++){
//         auto ind_row = vind(i,j,k,dim);
//         val = 1./2.;
//         triplet_dydz(trips,ind_row,i,  j,  k,  dimN,val);
//         val = 1./4.;
//         triplet_dydz(trips,ind_row,i-1,j,  k,dimN,val);
//         triplet_dydz(trips,ind_row,i+1,j,  k,dimN,val);   
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dx(Operator &oper, int dimN, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dx; 
//   for(int i = 1; i < size.x(); i++){
//     for(int j = 1; j < size.y(); j++){
//       for(int k = 1; k < size.z(); k++){    
//         auto ind_row = vind(i,j,k,dim);
//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k-1,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k-1,dimN),-val));

//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k  ,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k-1,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k-1,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k  ,dimN),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }
// Operator& Mesh::nodes_dy(Operator &oper, int dimN, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dy; 
//   for(int i = 1; i < size.x(); i++){
//     for(int j = 1; j < size.y(); j++){
//       for(int k = 1; k < size.z(); k++){    
//         auto ind_row = vind(i,j,k,dim);
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k-1,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k-1,dimN),-val));

//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k-1,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k-1,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k  ,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k  ,dimN),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::nodes_dz(Operator &oper, int dimN, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dz; 
//   for(int i = 1; i < size.x(); i++){
//     for(int j = 1; j < size.y(); j++){
//       for(int k = 1; k < size.z(); k++){    
//         auto ind_row = vind(i,j,k,dim);
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j-1,k-1,dimN),-val));

//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j-1,k-1,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i-1,j  ,k-1,dimN),-val));
        
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k  ,dimN),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k-1,dimN),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// void Mesh::stencil_curlE(){
//   auto rows = curlE.rows();
//   auto cols = curlE.cols();
//   Operator Ex(rows,cols);
//   Operator Ey(rows,cols);
//   Operator Ez(rows,cols);

//   curlE =  nodes_dy(Ez,2,0) - nodes_dz(Ey,1,0);
//   curlE += nodes_dz(Ex,0,1) - nodes_dx(Ez,2,1);
//   //curlE += nodes_dx(Ey,1,2) - nodes_dy(Ex,0,2); 

// }

// void Mesh::stencil_curlB(){
//   auto rows = curlB.rows();
//   auto cols = curlB.cols();
//   Operator Bx(rows,cols);
//   Operator By(rows,cols);
//   Operator Bz(rows,cols);

//   //curlB =  cells_dy(Bz,2,0) - cells_dz(By,1,0);
//   //curlB += cells_dz(Bx,0,1) - cells_dx(Bz,2,1);
//   curlB = cells_dx(By,1,2) - cells_dy(Bx,0,2); 

// }


// // Grid derivatives at cell centres
// // Makes the derivative of the oper at the cell centres
// // and returns the result at thegrid nodes

// // oper - operator on which the derivative is taken. 
// // dimN - dimension (for cell centres) from which the derivative is taken.
// // dim -  dimendion in which the result will be placed
// Operator& Mesh::cells_dx(Operator &oper, int dimC, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dx; 
//   for(int i = 0; i < size.x()-1; i++){
//     for(int j = 0; j < size.y()-1; j++){
//       for(int k = 0; k < size.z()-1; k++){    
//         auto ind_row = vind(i,j,k,dim);
        
//         trips.push_back(Trip(ind_row,vind(i+1,j  ,k  ,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k  ,dimC),-val));

//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k  ,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i,  j+1,k  ,dimC),-val));

//         trips.push_back(Trip(ind_row,vind(i+1,j,  k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k+1,dimC),-val));

//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i,  j+1,k+1,dimC),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

// Operator& Mesh::cells_dy(Operator &oper, int dimC, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dy; 
//   for(int i = 0; i < size.x()-1; i++){
//     for(int j = 0; j < size.y()-1; j++){
//       for(int k = 0; k < size.z()-1; k++){    
//         auto ind_row = vind(i,j,k,dim);

//         trips.push_back(Trip(ind_row,vind(i,  j+1,k  ,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i,  j,  k  ,dimC),-val));
        
//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k  ,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i+1,j  ,k  ,dimC),-val));        
        
//         trips.push_back(Trip(ind_row,vind(i,  j+1,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k+1,dimC),-val));

//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i+1,j,  k+1,dimC),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }


// Operator& Mesh::cells_dz(Operator &oper, int dimC, int dim)
// { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize*8);
//   double val;

//   val = 0.25 / Dz; 
//   for(int i = 0; i < size.x()-1; i++){
//     for(int j = 0; j < size.y()-1; j++){
//       for(int k = 0; k < size.z()-1; k++){    
//         auto ind_row = vind(i,j,k,dim);
//         trips.push_back(Trip(ind_row,vind(i,  j  ,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j  ,k,  dimC),-val));

//         trips.push_back(Trip(ind_row,vind(i+1,j  ,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i+1,j  ,k,  dimC),-val));
        
//         trips.push_back(Trip(ind_row,vind(i,  j+1,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i  ,j+1,k,  dimC),-val));
        
//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k+1,dimC),val));
//         trips.push_back(Trip(ind_row,vind(i+1,j+1,k,  dimC),-val));
//       }
//     }
//   }
//   oper.setFromTriplets(trips.begin(), trips.end());
//   return oper;
// }

//  void Mesh::stencil_Imat()
//  { 
//  // !!!!! needs bound condition and if cases!!!!!!
//   std::vector<Trip> trips;
//   const auto size = fieldE.size();
//   int totalSize = size.x()*size.y()*size.z();
//   trips.reserve(totalSize);
    
//   for(int i = 0; i < size.x(); i++){
//     for(int j = 0; j < size.y(); j++){
//       for(int k = 0; k < size.z(); k++){

//         // i,j,k
//         trips.push_back(Trip(vind(i,j,k,0),vind(i,j,k,0),1.0));

//         // i,j,k
//         trips.push_back(Trip(vind(i,j,k,1),vind(i,j,k,1),1.0));
          
//         // i,j,k
//         trips.push_back(Trip(vind(i,j,k,2),vind(i,j,k,2),1.0));
          
//       }
//     }
//   }
//   Imat.setFromTriplets(trips.begin(), trips.end());
// }