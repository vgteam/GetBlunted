#include "Cigar.hpp"

using std::stol;


namespace bluntifier {

///
/// M   I   D   N   S   H   P   =   X
/// 0   1   2   3   4   5   6   7   8
///

                                            // 0   1   2   3   4   5   6   7   8   9
const array<uint8_t, 128> Cigar::cigar_code = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 0
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 10
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 20
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 30
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 40
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 50
                                               10, 7,  10, 10, 10, 10, 10, 10, 2,  10,  // 60    =61, D68
                                               10, 10, 5,  1,  10, 10, 10, 0,  3,  10,  // 70    H72, I73, M77, N78
                                               6,  10, 10, 4,  10, 10, 10, 10, 8,  10,  // 80    P80, S83, X88
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 90
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 100
                                               10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 110
                                               10, 10, 10, 10, 10, 10, 10, 10};         // 120


const array<bool, 9> Cigar::is_ref_move = {true,      //MATCH      0
                                           false,     //INS        1
                                           true,      //DEL        2
                                           true,      //REF_SKIP   3
                                           false,     //SOFT_CLIP  4
                                           false,     //HARD_CLIP  5
                                           false,     //PAD        6
                                           true,      //EQUAL      7
                                           true};     //MISMATCH   8


const array<bool, 9> Cigar::is_query_move = {true,     //MATCH      0
                                             true,     //INS        1
                                             false,    //DEL        2
                                             false,    //REF_SKIP   3
                                             true,     //SOFT_CLIP  4
                                             false,    //HARD_CLIP  5
                                             false,    //PAD        6
                                             true,     //EQUAL      7
                                             true};    //MISMATCH   8


Cigar::Cigar(const string& s) {
    string token;

    for (auto& c: s) {
        if (isdigit(c)) {
            token += c;
        } else {
            operations.emplace_back(stol(token), c);
            token.clear();
        }
    }
}

}