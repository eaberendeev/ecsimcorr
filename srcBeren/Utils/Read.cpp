#include "Read.h"

//функция разделяющая string'овую переменную по делимитору delim
void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

void read_params_to_string(const std::string& target,const std::string& filename,std::vector< std::vector<std::string> > &params){
    std::ifstream fParams;
    std::string line;
    std::vector<std::string> strvec;
    int currentParam = -1;
	
	params.clear();

    fParams.open(filename, std::ifstream::in);
    getline(fParams, line);


    while (fParams.good() ){
        strvec = split(line, ' ');
        if(strvec[0] == target) {
            ++currentParam;
            params.resize(currentParam+1);
        }
        if (currentParam >= 0)           
        	params[currentParam].push_back(line);

        getline(fParams,line);
    }

    fParams.close();

}
