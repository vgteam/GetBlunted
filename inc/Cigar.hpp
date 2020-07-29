#ifndef BLUNTIFIER_CIGAR_HPP
#define BLUNTIFIER_CIGAR_HPP

#include <string>
#include <utility>
#include <vector>
#include <array>

using std::string;
using std::pair;
using std::vector;
using std::array;


namespace bluntifier{

class Cigar {
public:
    /// Attributes ///
    vector <pair<uint32_t, char>> operations;

    // Map from all possible chars to 0-8
    static const array<uint8_t,128> cigar_code;

    // Check if a cigar operation consumes a base in the reference or query
    static const array<bool,9> is_ref_move;
    static const array<bool,9> is_query_move;

    /// Methods ///
    Cigar(const string& s);
};


}

#endif //BLUNTIFIER_CIGAR_HPP
