#include "IncrementalIdMap.hpp"
#include "utility.hpp"
#include "Cigar.hpp"

using std::to_string;
using std::ifstream;
using std::cerr;

using bluntifier::parent_path;
using bluntifier::join_paths;
using bluntifier::AlignmentIterator;
using bluntifier::Alignment;
using bluntifier::Cigar;


int main(){
    string ref = "ACGT";
    string match = "ACGT";
    string singleton_match = "A";
    string mismatch = "ACAT";
    string ins = "ACAGT";
    string del = "ACT";

    string match_vs_ref = "4M";
    string singleton_match_vs_ref = "1M";
    string mismatch_vs_ref = "4M";
    string mismatch_vs_ref_explicit = "2=1X1=";
    string ins_vs_ref = "2M1I2M";
    string del_vs_ref = "2M1D1M";

    vector<string> test_names = {
            "match_vs_ref",
            "singleton_match_vs_ref",
            "mismatch_vs_ref",
            "mismatch_vs_ref_explicit",
            "ins_vs_ref",
            "del_vs_ref",
    };

    vector<Cigar> cigar_match_vs_ref = {Cigar(4,'M')};

    vector<Cigar> cigar_singleton_match_vs_ref = {Cigar(1,'M')};

    vector<Cigar> cigar_mismatch_vs_ref = {Cigar(4,'M')};

    vector<Cigar> cigar_mismatch_vs_ref_explicit = {Cigar(2,'='),
                                                    Cigar(1,'X'),
                                                    Cigar(1,'=')};

    vector<Cigar> cigar_ins_vs_ref = {Cigar(2,'M'),
                                      Cigar(1,'I'),
                                      Cigar(2,'M')};

    vector<Cigar> cigar_del_vs_ref = {Cigar(2,'M'),
                                      Cigar(1,'D'),
                                      Cigar(1,'M')};

    vector<uint64_t> query_index_match_vs_ref = {0,1,2,3};
    vector<uint64_t> query_index_singleton_match_vs_ref = {0};
    vector<uint64_t> query_index_mismatch_vs_ref = {0,1,2,3};
    vector<uint64_t> query_index_mismatch_vs_ref_explicit = {0,1,2,3};
    vector<uint64_t> query_index_ins_vs_ref = {0,1,2,3,4};
    vector<uint64_t> query_index_del_vs_ref = {0,1,1,2};

    vector<uint64_t> ref_index_match_vs_ref = {0,1,2,3};
    vector<uint64_t> ref_index_singleton_match_vs_ref = {0};
    vector<uint64_t> ref_index_mismatch_vs_ref = {0,1,2,3};
    vector<uint64_t> ref_index_mismatch_vs_ref_explicit = {0,1,2,3};
    vector<uint64_t> ref_index_ins_vs_ref = {0,1,1,2,3};
    vector<uint64_t> ref_index_del_vs_ref = {0,1,2,3};


    vector<string> cigar_strings = {
            match_vs_ref,
            singleton_match_vs_ref,
            mismatch_vs_ref,
            mismatch_vs_ref_explicit,
            ins_vs_ref,
            del_vs_ref
    };

    vector <vector <Cigar> > cigars = {
            cigar_match_vs_ref,
            cigar_singleton_match_vs_ref,
            cigar_mismatch_vs_ref,
            cigar_mismatch_vs_ref_explicit,
            cigar_ins_vs_ref,
            cigar_del_vs_ref
    };

    vector<string> query_sequences = {
            match,
            singleton_match,
            mismatch,
            mismatch,
            ins,
            del
    };

    vector <vector <uint64_t> > query_indexes = {
            query_index_match_vs_ref,
            query_index_singleton_match_vs_ref,
            query_index_mismatch_vs_ref,
            query_index_mismatch_vs_ref_explicit,
            query_index_ins_vs_ref,
            query_index_del_vs_ref
    };

    vector <vector <uint64_t> > ref_indexes = {
            ref_index_match_vs_ref,
            ref_index_singleton_match_vs_ref,
            ref_index_mismatch_vs_ref,
            ref_index_mismatch_vs_ref_explicit,
            ref_index_ins_vs_ref,
            ref_index_del_vs_ref
    };

    for (size_t i=0; i < cigar_strings.size(); i++){
        cerr << '\n' << test_names[i] << '\n' << '\n';

        Alignment alignment(cigar_strings[i]);
        AlignmentIterator iterator;

        size_t n = 0;
        while(alignment.step_through_alignment(iterator)){
            uint8_t cigar_code = alignment.operations[iterator.cigar_index].code;
            uint8_t cigar_length = alignment.operations[iterator.cigar_index].length;

            cerr << alignment.operations[iterator.cigar_index].type() << '\n';
            cerr << iterator << "\n\n";

            if (iterator.query_index != query_indexes[i][n]){
                throw runtime_error("FAIL: unexpected query index: " + to_string(cigar_code) +
                " should be " + to_string(query_indexes[i][n]));
            }

            if (iterator.ref_index != ref_indexes[i][n]){
                throw runtime_error("FAIL: unexpected ref index: " + to_string(cigar_length) +
                " should be " + to_string(ref_indexes[i][n]));
            }

            if (cigar_code != cigars[i][iterator.cigar_index].code){
                throw runtime_error("FAIL: unexpected cigar code: " + to_string(cigar_code) +
                " should be " + to_string(cigars[i][iterator.cigar_index].code));
            }

            if (cigar_length != cigars[i][iterator.cigar_index].length){
                throw runtime_error("FAIL: unexpected cigar length: " + to_string(cigar_length) +
                " should be " + to_string(cigars[i][iterator.cigar_index].length));
            }


            n++;
        }
        cerr << '\n';
    }

    for (size_t i=0; i < cigar_strings.size(); i++) {
        cerr << '\n' << test_names[i] << '\n' << '\n';

        Alignment alignment(cigar_strings[i]);

        cerr << alignment.create_formatted_alignment_string(ref, query_sequences[i]) << '\n';
    }

    cerr << '\n';

    Alignment non_explicit(mismatch_vs_ref);
    cerr << non_explicit << '\n';
    non_explicit.explicitize_mismatches(ref, mismatch);
    cerr << non_explicit << '\n';

    for (size_t i=0; i<cigar_mismatch_vs_ref_explicit.size(); i++){
        if (i>non_explicit.operations.size()){
            throw runtime_error("FAIL: explitized matches don't agree with truth set");
        }
        if (non_explicit.operations[i].code != cigar_mismatch_vs_ref_explicit[i].code){
            throw runtime_error("FAIL: explitized matches don't agree with truth set");
        }
        if (non_explicit.operations[i].length != cigar_mismatch_vs_ref_explicit[i].length){
            throw runtime_error("FAIL: explitized matches don't agree with truth set");
        }
    }

    cerr << "\nPASS\n";
    return 0;
}



