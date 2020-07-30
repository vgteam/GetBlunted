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


class Cigar{
public:
    /// Attributes ///
    uint8_t code;
    uint32_t length;

    /// Methods ///
    Cigar(uint32_t length, char type);
    char type() const;
};


class AlignmentIterator{
public:
    uint64_t query_index;
    uint64_t ref_index;
    uint64_t cigar_index;
    uint64_t intra_cigar_index;
    bool first_step;

    /// Methods ///
    AlignmentIterator();
    AlignmentIterator(uint64_t query_index, uint64_t ref_index);
};


class Alignment {
public:
    /// Attributes ///
    vector <Cigar> operations;

    // Map from all possible cigar chars to 0-8
    static const array<uint8_t,128> cigar_code;

    // Map back to the human readable characters
    static const array<char,9> cigar_type;

    // Check if a cigar operation consumes a base in the reference or query
    static const array<bool,9> is_ref_move;
    static const array<bool,9> is_query_move;

    /// Methods ///
    Alignment(const string& s);

    // Return the sequence lengths of the {query,ref} in a pair
    void compute_lengths(pair<uint64_t,uint64_t>& lengths);

    // Given a pair of start indexes {query_start,ref_start} walk through the alignment, updating the indexes in the
    // pair to retrace the coordinates of the path in the alignment matrix that produces this cigar. Return true as long
    // as the iteration is still progressing
    //
    // Intended usage:
    //     while(step_through_alignment(alignment_iterator)){
    //         ref_sequence[alignment_iterator.ref_index];
    //         query_sequence[alignment_iterator.query_index];
    //     }
    //
    bool step_through_alignment(AlignmentIterator& iterator);

    // Make a cool looking alignment string
    string generate_formatted_alignment_string(const string& ref_sequence, const string& query_sequence);
};


}

ostream& operator<<(ostream& o, bluntifier::Alignment& c);
ostream& operator<<(ostream& o, bluntifier::AlignmentIterator& c);

#endif //BLUNTIFIER_CIGAR_HPP
