#ifndef READ_H_
#define READ_H_
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
void read_params_to_string(const std::string& target,const std::string& filename,std::vector< std::vector<std::string> > &params);


#endif 
