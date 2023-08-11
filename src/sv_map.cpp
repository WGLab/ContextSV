#include "sv_map.h"

void SVMap::addSVCall(std::string chr, int start, int end, int sv_type)
{
    key_type key(chr, start, end, sv_type);

    // Check if the SV call has already been added
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        this->sv_calls[key]++;
    } else {
        this->sv_calls[key] = 1;
    }
}

// map_type SVMap::getSVCalls()
// {
//     return this->sv_calls;
// }
