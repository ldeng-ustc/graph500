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
#include "fcntl.h"
#include <omp.h>

#include "fmt/format.h"
#include "cxxopts.hpp"

#include "make_graph.h"

using namespace std;
namespace fs = std::filesystem;
int main(int argc, char* argv[]) {
    cxxopts::Options options("KronGenerator", "Generate Kron Graph with 2^n vertices and m*2^n edges");
    options
        .positional_help("[n] [m] [output_file]")
        .show_positional_help();
    options.add_options()
        ("n,log_numverts", "log2(#vertices)", cxxopts::value<int>()->default_value("16"))
        ("m,nedges_per_verts", "#edges per vertex", cxxopts::value<int>()->default_value("16"))
        ("o,outdir", "output directory", cxxopts::value<string>()->default_value("./"))
        ("file", "output file {0}, {1}, {2} will be replaced to n, m, and format (txt or bin)", 
                    cxxopts::value<string>()->default_value("Kron{0}-{1}.{2}"))
        ("f,format", "output format (0: stdout, 1:binary, 2:text)", cxxopts::value<int>()->default_value("0"))
        ("S,short", "use 32bit int as vertex ID in binary format")
        ("seed1", "user seed 1", cxxopts::value<uint64_t>()->default_value("1"))
        ("seed2", "user seed 2", cxxopts::value<uint64_t>()->default_value("2"))
        ("h,help", "print usage")
        ;
    options.parse_positional({"n", "m", "o"});
    auto opt = options.parse(argc, argv);
    if(opt.count("help")) {
        cout << options.help() << endl;
    } else {
        cout << fmt::format("Running in arguments:\n{}", opt.arguments_string()) << endl;
    }

    int log_numverts = opt["n"].as<int>();
    int64_t nedges_per_verts = opt["m"].as<int>();
    int64_t desired_nedges = nedges_per_verts << log_numverts;
    int64_t seed1 = opt["seed1"].as<uint64_t>();
    int64_t seed2 = opt["seed2"].as<uint64_t>();


    int64_t actual_nedges;
    packed_edge* result;

    /* Start of graph generation timing */
    double time_taken = omp_get_wtime();

    make_graph(log_numverts, desired_nedges, seed1, seed2, &actual_nedges, &result);

    time_taken = omp_get_wtime() - time_taken;
    /* End of graph generation timing */

    cout << fmt::format("{} edges generated in {}s ({} Medges/s)", actual_nedges, time_taken, 1e-6 * actual_nedges / time_taken  ) << endl;

    string format_string = opt["file"].as<string>();
    string filename;
    fs::path outdir = opt["outdir"].as<string>();
    fs::path fullpath;
    switch (opt["format"].as<int>()) {
    case 0: // stdout
        for (int i = 0; i < actual_nedges; i++) {
            int64_t v0 = get_v0_from_edge(result + i);
            int64_t v1 = get_v1_from_edge(result + i);
            cout << fmt::format("{:10}  {:10}", v0, v1) << endl;
        }
        break;
    case 1: // binary 
    {
        fs::create_directories(outdir);
        filename = fmt::format(format_string, log_numverts, nedges_per_verts, "bin");
        fullpath = outdir / filename;

        bool short_vid = opt["short"].as<bool>();
        size_t vsize = short_vid ? sizeof(uint32_t) : sizeof(uint64_t);
        size_t esize = 2ull * vsize;
        size_t fsize = esize * static_cast<size_t>(actual_nedges);

        int fd = open64(fullpath.c_str(), O_WRONLY | O_CREAT, 0664);
        fs::resize_file(fullpath, fsize);

        // parallel write to file.
        const size_t LOGN_BLOCK_SZ = 22;
        const size_t BLOCK_SZ = 1 << LOGN_BLOCK_SZ;
        size_t nblocks = max(actual_nedges >> LOGN_BLOCK_SZ, 1l);

        #pragma omp parallel for
        for(size_t b = 0; b < nblocks; b++) {
            size_t begin_edge = b * BLOCK_SZ;
            size_t end_edge = min( (b + 1) * BLOCK_SZ, static_cast<size_t>(actual_nedges));
            size_t pos = begin_edge * esize;
            //cout << fmt::format("writing block {:4}, start position: {:14}", b, pos) << endl;
            
            const int BUFFER_SIZE = 1024 * 1024;
            auto buffer = make_unique<uint64_t[]>(BUFFER_SIZE);
            int32_t *buffer32 = (int32_t*)(buffer.get());
            size_t n = short_vid? 2*BUFFER_SIZE : BUFFER_SIZE;
            size_t p = 0;
            for (size_t i = begin_edge; i < end_edge; i++) {
                int64_t v[2];
                v[0] = get_v0_from_edge(result + i);
                v[1] = get_v1_from_edge(result + i);
                if(opt["short"].as<bool>()) {
                    uint32_t v_short[2];
                    v_short[0] = static_cast<uint32_t>(v[0]);
                    v_short[1] = static_cast<uint32_t>(v[1]);
                    buffer32[p++] = v_short[0];
                    buffer32[p++] = v_short[1];
                } else {
                    buffer[p++] = v[0];
                    buffer[p++] = v[1];
                }
                if(p == n) {
                    pos += pwrite64(fd, buffer.get(), BUFFER_SIZE * sizeof(int64_t), pos);
                    p = 0;
                }
            }
            if(p != 0){
                p = short_vid ? p/2 : p;
                pos += pwrite64(fd, buffer.get(), p * sizeof(int64_t), pos);
            }
        }
        break;
    }
    case 2: // text
    {
        FILE *file = fopen(fullpath.c_str(), "w");
        fs::create_directories(outdir);
        filename = fmt::format(format_string, log_numverts, nedges_per_verts, "txt");
        fullpath = outdir / filename;
        for (int i = 0; i < actual_nedges; i++) {
            int64_t v[2];
            v[0] = get_v0_from_edge(result + i);
            v[1] = get_v1_from_edge(result + i);
            fputs(fmt::format("{} {}\n", v[0], v[1]).c_str(), file);
        }
        break;
    }
    default:
        cout << "wrong format." << endl;
        cout << options.help() << endl;
        break;
    }

    free(result);

    return 0;
}
