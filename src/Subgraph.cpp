#include "Subgraph.hpp"


namespace bluntifier {

PathInfo::PathInfo(path_handle_t path_handle, uint32_t spoa_id, bool biclique_side) :
        path_handle(path_handle),
        spoa_id(spoa_id),
        biclique_side(biclique_side) {}


}