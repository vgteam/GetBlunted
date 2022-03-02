#include "Subgraph.hpp"


namespace bluntifier {

PathInfo::PathInfo(path_handle_t path_handle, uint32_t poa_id, bool biclique_side) :
        path_handle(path_handle),
        poa_id(poa_id),
        biclique_side(biclique_side) {}


}
