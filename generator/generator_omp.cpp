/* Copyright (C) 2009-2010 The Trustees of Indiana University.             */
/*                                                                         */
/* Use, modification and distribution is subject to the Boost Software     */
/* License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at */
/* http://www.boost.org/LICENSE_1_0.txt)                                   */
/*                                                                         */
/*  Authors: Jeremiah Willcock                                             */
/*           Andrew Lumsdaine                                              */

#include <iostream>
#include <filesystem>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstdio>

#include <unistd.h>
#include <omp.h>

#include "fcntl.h"
#include "sys/mman.h"
#include "fmt/format.h"
#include "cxxopts.hpp"

#include "make_graph.h"
#include "utils.h"

using namespace std;
namespace fs = std::filesystem;

void write_to_stdout(packed_edge* result, size_t m) {
    for (size_t i = 0; i < m; i++) {
        int64_t v0 = get_v0_from_edge(result + i);
        int64_t v1 = get_v1_from_edge(result + i);
        cout << fmt::format("{:10} {:10}", v0, v1) << endl;
    }
}

template<typename VertexType>
void write_to_file_binary(fs::path path, packed_edge* result, size_t nedges) {
    size_t vsize = sizeof(VertexType);
    size_t esize = 2ull * vsize;
    size_t fsize = esize * nedges;
    int fd = open64(path.c_str(), O_WRONLY | O_CREAT, 0664);
    fs::resize_file(path, fsize);
    VertexType *file = static_cast<VertexType*>(mmap(NULL, fsize, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0));
    if(file == MAP_FAILED) {
        cout << "mmap failed." << endl;
        cout << std::strerror(errno) << endl;
        exit(errno);
    }
    // parallel write to file.
    const size_t LOGN_BLOCK_SZ = 22;
    const size_t BLOCK_SZ = 1 << LOGN_BLOCK_SZ;

    #pragma omp parallel for schedule(static, BLOCK_SZ)
    for(size_t i=0; i<nedges; i++) {
        *(file + 2*i) = static_cast<VertexType>(get_v0_from_edge(result + i));
        *(file + 2*i + 1) = static_cast<VertexType>(get_v1_from_edge(result + i));
    }
    munmap(file, fsize);
    close(fd);
}

void write_to_file_text(fs::path path, packed_edge* result, size_t nedges) {
    FILE *file = fopen(path.c_str(), "w");
    for (size_t i = 0; i < nedges; i++) {
        int64_t v[2];
        v[0] = get_v0_from_edge(result + i);
        v[1] = get_v1_from_edge(result + i);
        fputs(fmt::format("{} {}\n", v[0], v[1]).c_str(), file);
    }
    fclose(file);
}

int main(int argc, char* argv[]) {
    cxxopts::Options options("KronGenerator", "Generate Kron Graph with 2^n vertices and m*2^n edges");
    options
        .positional_help("[n] [m] [output_file]")
        .show_positional_help();
    options.add_options()
        ("n,log_numverts", "log2(#vertices)", cxxopts::value<int>()->default_value("16"))
        ("m,nedges_per_verts", "#edges per vertex", cxxopts::value<int>()->default_value("16"))
        ("o,path", "output file path, {n} is the wildcards and pass to the fmt::format, replacement rule: \n"
                    "{0}: log_numverts\n{1}: nedges_per_verts\n"
                    "{2}: data format, 'txt' for text, 'bin' for binary\n"
                    "{3}: file number, necessary when graph is large than filesize",
                    cxxopts::value<string>()->default_value("/data/Kron/Kron{0}-{1}/block-{3:02}.{2}"))
        ("b,log_blocksize", "max number of edges be generated in a iteration, must fit in memory, also the max single file size",
                    cxxopts::value<int>()->default_value("30"))
        ("f,format", "output format (0: stdout, 1:binary, 2:text)", cxxopts::value<int>()->default_value("0"))
        ("S,short", "use 32bit int as vertex ID in binary format")
        ("v,info", "show debug messages")
        ("seed1", "user seed 1", cxxopts::value<uint64_t>()->default_value("1"))
        ("seed2", "user seed 2", cxxopts::value<uint64_t>()->default_value("2"))
        ("h,help", "print usage")
        ;
    options.parse_positional({"n", "m", "o"});
    auto opt = options.parse(argc, argv);
    if(opt.count("help")) {
        cout << options.help() << endl;
        exit(0);
    } else {
        cout << fmt::format("Running in arguments:\n{}", opt.arguments_string()) << endl;
    }

    string path_format = opt["path"].as<string>();
    int log_numverts = opt["n"].as<int>();
    int64_t nedges_per_verts = opt["m"].as<int>();
    int64_t desired_nedges = nedges_per_verts << log_numverts;
    int64_t seed1 = opt["seed1"].as<uint64_t>();
    int64_t seed2 = opt["seed2"].as<uint64_t>();
    int64_t block_size = 1l << opt["log_blocksize"].as<int>();
    int64_t nblocks = desired_nedges / block_size + (desired_nedges % block_size != 0);

    uint_fast32_t seed[5];
    make_mrg_seed(seed1, seed2, seed);

    for(int64_t i=0; i < nblocks; i++) {
        int64_t start_edge = i * block_size;
        int64_t end_edge = min((i+1)*block_size, desired_nedges);
        size_t nblock_edges = static_cast<size_t>(end_edge - start_edge);
        packed_edge* edges = (packed_edge*)xmalloc(block_size * sizeof(packed_edge));

        if(opt["info"].as<bool>()) {
            cout << fmt::format("Generating block {}, range [{}, {})", i, start_edge, end_edge) << endl;
        }

        /* Start of graph generation timing */
        double time_taken = omp_get_wtime();
        generate_kronecker_range(seed, log_numverts, start_edge, end_edge, edges);
        time_taken = omp_get_wtime() - time_taken;
        /* End of graph generation timing */

        if(opt["info"].as<bool>()) {
            cout << fmt::format("{} edges generated in {}s ({} Medges/s)", nblock_edges, time_taken, 1e-6 * nblock_edges / time_taken  ) << endl;
        }

        time_taken = omp_get_wtime();
        fs::path path;
        switch (opt["format"].as<int>()) {
        case 0: // stdout
            write_to_stdout(edges, nblock_edges);
            break;
        case 1: // binary
            path = fmt::format(path_format, log_numverts, nedges_per_verts, "bin", i);
            fs::create_directories(path.parent_path());
            if(opt["short"].as<bool>()){
                write_to_file_binary<uint32_t>(path, edges, nblock_edges);
            } else {
                write_to_file_binary<int64_t>(path, edges, nblock_edges);
            }
            break;
        case 2: // text
            path = fmt::format(path_format, log_numverts, nedges_per_verts, "txt", i);
            fs::create_directories(path.parent_path());
            write_to_file_text(path, edges, nblock_edges);
            break;
        default:
            cout << "wrong format." << endl;
            cout << options.help() << endl;
            break;
        }
        time_taken = omp_get_wtime() - time_taken;
        /* End of graph writing timing */

        if(opt["info"].as<bool>()) {
            cout << fmt::format("{} edges written in {}s ({} Medges/s)", nblock_edges, time_taken, 1e-6 * nblock_edges / time_taken  ) << endl;
        }
    }

    return 0;
}
