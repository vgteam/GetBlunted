#include "Cigar.hpp"

using std::runtime_error;
using std::to_string;
using std::stol;


namespace bluntifier {

///
/// M   I   D   N   S   H   P   =   X
/// 0   1   2   3   4   5   6   7   8
///


                                                // 0   1   2   3   4   5   6   7   8   9
const array<uint8_t, 128> Alignment::cigar_code = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 0
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


const array<bool, 9> Alignment::is_ref_move = {true,      //MATCH      0
                                           false,     //INS        1
                                           true,      //DEL        2
                                           true,      //REF_SKIP   3
                                           false,     //SOFT_CLIP  4
                                           false,     //HARD_CLIP  5
                                           false,     //PAD        6
                                           true,      //EQUAL      7
                                           true};     //MISMATCH   8


const array<bool, 9> Alignment::is_query_move = {true,     //MATCH      0
                                             true,     //INS        1
                                             false,    //DEL        2
                                             false,    //REF_SKIP   3
                                             true,     //SOFT_CLIP  4
                                             false,    //HARD_CLIP  5
                                             false,    //PAD        6
                                             true,     //EQUAL      7
                                             true};    //MISMATCH   8

const array<char, 9> Alignment::cigar_type = {'M','I','D','N','S','H','P','=','X'};


Cigar::Cigar(uint32_t length, char type):
        code(Alignment::cigar_code[type]),
        length(length)
{
    if (code > 8){
        throw runtime_error("ERROR: unrecognized cigar character: " + string(1,type) + " has ASCII value: " + to_string(int(code)));
    }
}

char Cigar::type() const{
    return Alignment::cigar_type[code];
}

AlignmentIterator::AlignmentIterator(uint64_t ref_index, uint64_t query_index):
    query_index(query_index),
    ref_index(ref_index),
    cigar_index(0),
    intra_cigar_index(0),
    first_step(true)
{}


AlignmentIterator::AlignmentIterator():
    query_index(0),
    ref_index(0),
    cigar_index(0),
    intra_cigar_index(0),
    first_step(true)
{}


void AlignmentIterator::next_cigar(){
    cigar_index++;
    intra_cigar_index = 0;
}


Alignment::Alignment(const string& s) {
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


void Alignment::compute_lengths(pair<size_t,size_t>& lengths) const{
    // Reset lengths
    lengths.first = 0;
    lengths.second = 0;

    // Count up the cigar operations
    for (const auto& c: operations){
        const uint8_t code = Alignment::cigar_code[c.type()];

        // Assume the left side (source) node is treated as the "reference" in the cigar
        if (Alignment::is_ref_move[code]){
            // Increment by the length of the cigar operation
            lengths.first += c.length;
        }

        // Assume the right side (sink) node is treated as the "query" in the cigar
        if (Alignment::is_query_move[code]){
            // Increment by the length of the cigar operation
            lengths.second += c.length;
        }
    }
}


uint64_t Alignment::compute_common_length(){
    uint64_t n_matches = 0;

    // Count up the cigar operations
    for (const auto& c: operations){
        const uint8_t code = Alignment::cigar_code[c.type()];

        if (Alignment::is_ref_move[code] and Alignment::is_query_move[code]){
            n_matches += c.length;
        }
    }

    return n_matches;
}


bool Alignment::step_through_alignment(AlignmentIterator& iterator){
    // Don't do anything on the first step
    if (iterator.first_step){
        iterator.first_step = false;

        // Do a quick check to see if this is actually a non-overlap (0M). If so, don't iterate
        if (operations.empty() or (operations.size() == 1 and operations[0].length == 0)){
            return false;
        }
        else {
            return true;
        }
    }

    bool is_last_cigar = (iterator.cigar_index >= operations.size() - 1);
    bool is_last_step_in_cigar = (iterator.intra_cigar_index >= operations[iterator.cigar_index].length - 1);
    bool done = is_last_cigar and is_last_step_in_cigar;

    // Exit early if this is the end my friend
    if (done){
        return false;
    }

    // Move up to the next cigar tuple if necessary
    if (is_last_step_in_cigar){
        iterator.next_cigar();
    }
    // Or just increment the intra_cigar_index
    else{
        iterator.intra_cigar_index++;
    }

    // Ref
    if (is_ref_move[operations[iterator.cigar_index].code]){
        iterator.ref_index++;
    }

    // Query
    if (is_query_move[operations[iterator.cigar_index].code]){
        iterator.query_index++;
    }


    return true;
}


vector<Cigar> Alignment::explicitize_mismatches(
        const string& query_sequence,
        const string& ref_sequence,
        uint64_t ref_start_index,
        uint64_t query_start_index){

    // TODO: rewrite this function without copying? Use insert operations instead

    AlignmentIterator iterator(query_start_index, ref_start_index);
    vector<Cigar> explicit_operations;

    while (step_through_alignment(iterator)) {
        if (operations[iterator.cigar_index].type() == 'M'){
            char ref_base = ref_sequence[iterator.ref_index];
            char query_base = query_sequence[iterator.query_index];

            if (ref_base == query_base){
                // If the last operation was already a Match (=), then extend its length
                if (not explicit_operations.empty() and explicit_operations.back().type() == '='){
                    explicit_operations.back().length++;
                }
                // Otherwise make a new operation of type =
                else{
                    explicit_operations.emplace_back(1,'=');
                }
            }
            else{
                // If the last operation was already a Mismatch (X), then extend its length
                if (not explicit_operations.empty() and explicit_operations.back().type() == 'X'){
                    explicit_operations.back().length++;
                }
                // Otherwise make a new operation of type X
                else{
                    explicit_operations.emplace_back(1,'X');
                }
            }
        }
        else{
            // Just copy any non-match operations
            explicit_operations.emplace_back(operations[iterator.cigar_index]);

            // Skip iterating the remaining coordinates in this cigar if its not a Match operation
            iterator.next_cigar();
        }
    }

    return explicit_operations;
}


vector<Cigar> Alignment::explicitize_mismatches(
        const HandleGraph& graph,
        const edge_t& edge,
        uint64_t ref_start_index,
        uint64_t query_start_index) {

    // TODO: rewrite this function without copying? Use insert operations instead

    AlignmentIterator iterator(query_start_index, ref_start_index);
    vector<Cigar> explicit_operations;

    while (step_through_alignment(iterator)) {
        if (operations[iterator.cigar_index].type() == 'M'){
            char query_base = graph.get_base(edge.second, iterator.query_index);
            char ref_base = graph.get_base(edge.first, iterator.ref_index);

            if (ref_base == query_base){
                // If the last operation was already a Match (=), then extend its length
                if (not explicit_operations.empty() and explicit_operations.back().type() == '='){
                    explicit_operations.back().length++;
                }
                    // Otherwise make a new operation of type =
                else{
                    explicit_operations.emplace_back(1,'=');
                }
            }
            else{
                // If the last operation was already a Mismatch (X), then extend its length
                if (not explicit_operations.empty() and explicit_operations.back().type() == 'X'){
                    explicit_operations.back().length++;
                }
                    // Otherwise make a new operation of type X
                else{
                    explicit_operations.emplace_back(1,'X');
                }
            }
        }
        else{
            // Just copy any non-match operations
            explicit_operations.emplace_back(operations[iterator.cigar_index]);

            // Skip iterating the remaining coordinates in this cigar if its not a Match operation
            iterator.next_cigar();
        }
    }

    return explicit_operations;
}


string Alignment::create_formatted_alignment_string(
        const HandleGraph& graph,
        const edge_t& edge,
        uint64_t ref_start_index,
        uint64_t query_start_index
        ) {

    AlignmentIterator iterator(ref_start_index, query_start_index);

    string aligned_ref;
    string aligned_query;
    string alignment_symbols;

    while (step_through_alignment(iterator)) {
        char ref_base = graph.get_base(edge.first, iterator.ref_index);
        char query_base = graph.get_base(edge.second, iterator.query_index);

        uint8_t code = operations[iterator.cigar_index].code;

        if (Alignment::is_ref_move[code]) {
            aligned_ref += ref_base;
        } else {
            aligned_ref += "-";
        }

        if (Alignment::is_query_move[code]) {
            aligned_query += query_base;
        } else {
            aligned_query += "-";
        }

        if (Alignment::is_ref_move[code] and Alignment::is_query_move[code] and ref_base == query_base) {
            alignment_symbols += '|';
        } else {
            if (cigar_type[code] == 'M') {
                alignment_symbols += "*";
            }
            else{
                alignment_symbols += " ";
            }
        }
    }

    return aligned_ref + '\n' + alignment_symbols + '\n' + aligned_query;
}


string Alignment::create_formatted_alignment_string(
        const string& ref_sequence,
        const string& query_sequence,
        uint64_t ref_start_index,
        uint64_t query_start_index
        ) {

    AlignmentIterator iterator(ref_start_index, query_start_index);

    string aligned_ref;
    string aligned_query;
    string alignment_symbols;

    while (step_through_alignment(iterator)) {
        char ref_base = ref_sequence[iterator.ref_index];
        char query_base = query_sequence[iterator.query_index];

        uint8_t code = operations[iterator.cigar_index].code;

        if (Alignment::is_ref_move[code]) {
            aligned_ref += ref_base;
        } else {
            aligned_ref += "-";
        }

        if (Alignment::is_query_move[code]) {
            aligned_query += query_base;
        } else {
            aligned_query += "-";
        }

        if (Alignment::is_ref_move[code] and Alignment::is_query_move[code] and ref_base == query_base) {
            alignment_symbols += '|';
        } else {
            if (cigar_type[code] == 'M') {
                alignment_symbols += "*";
            }
            else{
                alignment_symbols += " ";
            }
        }

    }

    return aligned_ref + '\n' + alignment_symbols + '\n' + aligned_query;
}


}

ostream& operator<<(ostream& o, bluntifier::Alignment& a){
    for (const auto& cigar: a.operations) {
        o << '(' << cigar.length << ',' << cigar.type() << ')' << ',';
    }

    return o;
}


ostream& operator<<(ostream& o, bluntifier::AlignmentIterator& a){
    o << "cigar_index:\t\t" << a.cigar_index;
    o << "\nintra_cigar_index:\t" << a.intra_cigar_index;
    o << "\nref_index:\t\t" << a.ref_index;
    o << "\nquery_index:\t\t" << a.query_index;
    o << "\nfirst_step:\t\t" << a.first_step;

    return o;
}
