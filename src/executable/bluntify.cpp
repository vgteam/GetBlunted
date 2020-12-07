#include "Bluntifier.hpp"

using bluntifier::Bluntifier;


int main(int argc, char **argv){
    string gfa_path;

    if (argc == 1){
        throw runtime_error("No input gfa path provided");
    }
    else if (argc == 2){
        gfa_path = argv[1];
    }
    else{
        throw runtime_error("Too many arguments. Specify 1 input gfa path.");
    }

    Bluntifier bluntifier(gfa_path);
    bluntifier.bluntify();

    return 0;
}

