// Region data (chr, start, end)

#ifndef REGION_H
#define REGION_H

/// @cond
#include <string>
/// @endcond

struct Region {
    std::string chr;
    int start;
    int end;
    
    // Return the region as a string, ignore the start and end positions if they
    // are -1
    std::string toString() {
        if (this->start == -1 && this->end == -1) {
            return this->chr;
        } else {
            return this->chr + ":" + std::to_string(this->start) + "-" + std::to_string(this->end);
        }
    }

    Region(std::string chr, int start, int end):
        chr(chr), start(start), end(end) {}

    Region():
        chr(""), start(-1), end(-1) {}
};

#endif // REGION_H
