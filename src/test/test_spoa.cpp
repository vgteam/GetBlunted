#include "spoa/alignment_engine.hpp"
#include "spoa/graph.hpp"
#include <string>

using spoa::AlignmentEngine;
using spoa::AlignmentType;
using spoa::Graph;
using std::string;


int main(){

    std::vector<std::string> sequences = {
            "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
            "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
            "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
            "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
            "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
            "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    auto alignment_engine = spoa::AlignmentEngine::Create(
            spoa::AlignmentType::kNW, 3, -5, -3);  // linear gaps

    spoa::Graph graph{};

    for (const auto& it : sequences) {
        auto alignment = alignment_engine->Align(it, graph);
        graph.AddAlignment(alignment, it);
    }

    auto consensus = graph.GenerateConsensus();

    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
              << consensus << std::endl;

    auto msa = graph.GenerateMultipleSequenceAlignment();

    for (const auto& it : msa) {
        std::cerr << it << std::endl;
    }



    return 0;
}