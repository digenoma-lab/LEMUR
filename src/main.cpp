#include "merge_bedmethyl/merger.hpp"
#include "merge_bedmethyl/options.hpp"

int main(int argc, char* argv[]) {
    merge_bedmethyl::ParsedArgs args;
    const int parse_rc = merge_bedmethyl::parse_arguments(argc, argv, args);
    if (parse_rc != 0) return parse_rc;
    if (args.show_help) return 0;
    return merge_bedmethyl::run_merge(args);
}
