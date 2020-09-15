#include "Pileup.hpp"


namespace bluntifier{


const char Pileup::space = '_';


Pileup::Pileup():
        id_map(true)
{}


void Pileup::to_string(string& s){
    s.clear();

    for (size_t j=0; j<matrix[0].size(); j++){
        for (size_t i=0; i<matrix.size(); i++){
            s += matrix[i][j];
        }
        s += '\n';
    }
}


}