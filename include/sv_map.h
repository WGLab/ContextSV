#ifndef SVMAP_H
#define SVMAP_H

#include <string>
#include <map>

class SVMap {
    private:
        // Map of SV calls
        // Key:   chr, start, end, sv_type
        // Value: read count
        std::map<std::tuple<std::string, int, int, int>, int> sv_calls;

    public:
        void addSVCall(std::string, int start, int end, int sv_type);
        std::map<std::tuple<std::string, int, int, int>, int> getSVCalls();
};

#endif // SVMAP_H
