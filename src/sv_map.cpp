#include "sv_map.h"

void SVMap::addSVCall(std::string chr, int start, int end, int sv_type)
{
    std::tuple<std::string, int, int, int> key(chr, start, end, sv_type);

    // Check if the SV call has already been added
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        this->sv_calls[key]++;
    } else {
        this->sv_calls[key] = 1;
    }
}

std::map<std::tuple<std::string, int, int, int>, int> SVMap::getSVCalls()
{
    return this->sv_calls;
}
