#ifndef SVMAP_H
#define SVMAP_H

#include <string>
#include <map>

class SVMap {
    private:
        // Map aliases
        typedef std::map<std::tuple<std::string, int, int, int>, int> map_type;

        // Key aliases
        typedef std::tuple<std::string, int, int, int> key_type;

        // Map of SV calls
        // Key:   chr, start, end, sv_type
        // Value: read count
        map_type sv_calls;

    public:
        void addSVCall(std::string, int start, int end, int sv_type);
        //map_type getSVCalls();
};

#endif // SVMAP_H
