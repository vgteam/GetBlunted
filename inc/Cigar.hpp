#ifndef BLUNTIFIER_CIGAR_HPP
#define BLUNTIFIER_CIGAR_HPP

#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <array>

using std::ostream;
using std::vector;
using std::string;
using std::pair;
using std::array;


namespace bluntifier{


class CigarOperation{
public:
    /// Attributes ///
    char type;
    uint32_t length;

    /// Methods ///
    CigarOperation(uint32_t length, char type);
};


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

ostream& operator<<(ostream& o, bluntifier::Cigar& c);

#endif //BLUNTIFIER_CIGAR_HPP
