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

void SVMap::addSVCalls(SVMap sv_calls)
{
    // Iterate over the SV calls
    for (auto const& sv_call : sv_calls.sv_calls) {
        // Add the SV call to the map
        this->addSVCall(std::get<0>(sv_call.first), std::get<1>(sv_call.first), std::get<2>(sv_call.first), std::get<3>(sv_call.first));
    }
}
