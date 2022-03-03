#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include "abpoa.h"
#include "bdsg/hash_graph.hpp"
#include "unchop.hpp"
#include "sdsl/bits.hpp"

using namespace std;
using namespace bdsg;
using namespace bluntifier;
/**
 * copied directly from https://github.com/yangao07/abPOA/blob/main/example.c
 */

// for nt
// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nt4_table[256] = {
    0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
    'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
};

int main(){

    int i, j, n_seqs = 10;
    char seqs[10][100] = {
        "CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
        "CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
        "CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
        "CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
        "CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
        "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
        "CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT"
    };
    
    vector<string> seq_names{
        "seq0",
        "seq1",
        "seq2",
        "seq3",
        "seq4",
        "seq5",
        "seq6",
        "seq7",
        "seq8",
        "seq9"
    };
    
    // convert to c style names, sigh...
    char** c_seq_names = (char**) malloc(sizeof(char*) * seq_names.size());
    for (size_t i = 0; i < seq_names.size(); ++i) {
        c_seq_names[i] = (char*) malloc(sizeof(char) * (seq_names[i].size() + 1));
        for (size_t j = 0; j < seq_names[i].size(); ++j) {
            c_seq_names[i][j] = seq_names[i][j];
        }
        c_seq_names[i][seq_names[i].size()] = '\0';
    }
    
    
    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;
    
    abpoa_post_set_para(abpt);
    
    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = strlen(seqs[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j)
            bseqs[i][j] = nt4_table[(int)seqs[i][j]];
    }
    
    // 1. output to stdout
    fprintf(stdout, "=== output to stdout ===\n");
    
    // perform abpoa-msa
    abpoa_msa(ab, abpt, n_seqs, c_seq_names, seq_lens, bseqs, stdout, NULL, NULL, NULL, NULL, NULL, NULL);
    
    fprintf(stdout, "=== testing my conversion algorithm ===\n");
    
    unordered_map<int, handle_t> abpoa_node_to_handle;
    HashGraph graph;
    abpoa_graph_t* ab_graph = ab->abg;
    
    // get all of the nodes
    // note: 0 and 1 are reserved for the source and sink nodes (which are only
    // for internal use, no associated sequence)
    for (int j = 2; j < ab_graph->node_n; ++j) {
        string seq(1, nt256_table[ab_graph->node[j].base]);
        handle_t handle = graph.create_handle(seq);
        abpoa_node_to_handle[j] = handle;
    }
    // get all of the edges
    for (int j = 2; j < ab_graph->node_n; ++j) {
        for (int k = 0; k < ab_graph->node[j].out_edge_n; ++k) {
            int next_id = ab_graph->node[j].out_id[k];
            if (next_id == ABPOA_SRC_NODE_ID || next_id == ABPOA_SINK_NODE_ID) {
                continue;
            }
            graph.create_edge(abpoa_node_to_handle[j],
                              abpoa_node_to_handle[next_id]);
        }
    }
    
    
    // figure out the correspondence between abpoa's integer ids and the path
    // handles that we added during the init step
    abpoa_seq_t* ab_seqs = ab->abs;
    vector<path_handle_t> read_id_to_path;
    for (uint64_t j = 0; j < ab_seqs->n_seq; ++j) {
        
        assert(ab_seqs->name[j].l != 0);
        
        string path_name(ab_seqs->name[j].s);
        read_id_to_path.push_back(graph.create_path_handle(path_name));
    }
        
    // Kahn's algorithm to compute topological order
    vector<int> in_degree(ab_graph->node_n);
    vector<int> stack;
    for (int j = 0; j < ab_graph->node_n; ++j) {
        in_degree[j] = ab_graph->node[j].in_edge_n;
        if (in_degree[j] == 0) {
            stack.push_back(j);
        }
    }
    vector<int> topological_order;
    topological_order.reserve(in_degree.size());
    while (!stack.empty()) {
        int here = stack.back();
        stack.pop_back();
        topological_order.push_back(here);
        
        for (int j = 0; j < ab_graph->node[here].out_edge_n; ++j) {
            int next = ab_graph->node[here].out_id[j];
            --in_degree[next];
            if (in_degree[next] == 0) {
                stack.push_back(next);
            }
        }
    }
    
    // add paths for the input sequences
    for (int nid : topological_order) {
        if (nid == ABPOA_SRC_NODE_ID || nid == ABPOA_SINK_NODE_ID) {
            continue;
        }
        // add node id to read path
        int b = 0;
        for (size_t j = 0; j < ab_graph->node[nid].read_ids_n; ++j) {
            uint64_t num = ab_graph->node[nid].read_ids[j];
            uint64_t tmp;
            while (num) {
                tmp = num & -num;
                int read_id = sdsl::bits::hi(tmp);
                graph.append_step(read_id_to_path[b+read_id], abpoa_node_to_handle.at(nid));
                num ^= tmp;
            }
            b += 64;
        }
        
//        // TODO: the abPOA function actually has a much more complicated
//        // indexing scheme, hopefully this is enough
//        cerr << "checking " << ab_graph->node[nid].read_ids_n << " reads on abpoa " << nid << ", handle " << graph.get_id(abpoa_node_to_handle.at(nid)) << endl;
//        for (int j = 0; j < ab_graph->node[nid].read_ids_n; ++j) {
//            uint64_t read_id = ab_graph->node[nid].read_ids[j];
//            cout << "\tgot read ID " << read_id << endl;
//            path_handle_t path_handle = read_id_to_path.at(read_id);
//            graph.append_step(path_handle, abpoa_node_to_handle.at(nid));
//        }
    }
    
    unchop(&graph);
    
    
    vector<handle_t> handles;
    graph.for_each_handle([&](const handle_t& h) {
        handles.push_back(h);
    });
    sort(handles.begin(), handles.end(), [&](handle_t a, handle_t b) {
        return graph.get_id(a) < graph.get_id(b);
    });
    
    for (auto h : handles) {
        cout << graph.get_id(h) << " " << graph.get_sequence(h) << endl;
        for (bool left : {true, false}) {
            graph.follow_edges(h, left, [&](const handle_t& n) {
                cout << '\t';
                if (!left) {
                    cout << "-> ";
                }
                cout << graph.get_id(n);
                if (left) {
                    cout << " <-";
                }
                cout << endl;
            });
        }
    }
    
    
    graph.for_each_path_handle([&](const path_handle_t& p) {
        cout << graph.get_path_name(p) << endl;
        bool first = true;
        for (auto h : graph.scan_path(p)) {
            if (!first){
                cout << ", ";
            }
            cout << graph.get_id(h);
            first = false;
        }
        cout << endl;
    });
    
    
    
    // 2. variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n=0;
    uint8_t **msa_seq; int msa_l=0;
    
    // perform abpoa-msa
    ab->abs->n_seq = 0; // To re-use ab, n_seq needs to be set as 0
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);
    
    fprintf(stdout, "=== output to variables ===\n");
    for (i = 0; i < cons_n; ++i) {
        fprintf(stdout, ">Consensus_sequence\n");
        for (j = 0; j < cons_l[i]; ++j)
            fprintf(stdout, "%c", nt256_table[cons_seq[i][j]]);
        fprintf(stdout, "\n");
    }
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < n_seqs; ++i) {
        for (j = 0; j < msa_l; ++j) {
            fprintf(stdout, "%c", nt256_table[msa_seq[i][j]]);
        }
        fprintf(stdout, "\n");
    }
    
    if (cons_n) {
        for (i = 0; i < cons_n; ++i) {
            free(cons_seq[i]); free(cons_cov[i]);
        } free(cons_seq); free(cons_cov); free(cons_l);
    }
    if (msa_l) {
        for (i = 0; i < n_seqs; ++i) free(msa_seq[i]); free(msa_seq);
    }
    
    /* generate DOT partial order graph plot */
    abpt->out_pog = strdup("example.png"); // dump parital order graph to file
    if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);
    
    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
    abpoa_free(ab); abpoa_free_para(abpt);

    return 0;
}
