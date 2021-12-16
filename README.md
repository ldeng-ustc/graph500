# Graph500-3.0.0

## Generator

This fork add an openmp generator to generate Kron graph.

usage:

```sh
cd generator
make
./generator_omp 20 16 ./Kron-26
```

or use `-h` to see help as below:

```sh
Generate Kron Graph with 2^n vertices and m*2^n edges
Usage:
  KronGenerator [OPTION...] [n] [m] [output_file]

  -n, --log_numverts arg      log2(#vertices) (default: 16)
  -m, --nedges_per_verts arg  #edges per vertex (default: 16)
  -o, --path arg              output file path, {n} is the wildcards and
                              pass to the fmt::format, replacement rule:
                              {0}: log_numverts
                              {1}: nedges_per_verts
                              {2}: data format, 'txt' for text, 'bin' for
                              binary
                              {3}: file number, necessary when graph is
                              large than filesize (default:
                              /data/Kron/Kron{0}-{1}/block-{3:02}.{2})
  -b, --log_blocksize arg     max number of edges be generated in a
                              iteration, must fit in memory, also the max
                              single file size (default: 30)
  -f, --format arg            output format (0: stdout, 1:binary, 2:text)
                              (default: 0)
  -S, --short                 use 32bit int as vertex ID in binary format
  -v, --info                  show debug messages
      --seed1 arg             user seed 1 (default: 1)
      --seed2 arg             user seed 2 (default: 2)
  -h, --help                  print usage
```

## Original Readme


Compiling should be pretty straightforward as long as you have a valid MPI-3 library loaded in your PATH.
There is no more OpenMP,Sequential and XMT versions of benchmark.

On single node you can run MPI code with reasonable performance.

To build binaries change directory to src and execute make.
If you are lucky four binaries would be built, two of which are of interest:

graph500_reference_bfs runs BFS kernel (and skips weight generation)
graph500_reference_bfs_sssp runs both BFS and SSSP kernels

Both binaries require one integer parameter which is scale of the graph.
Validation can be deactivated by specifying SKIP_VALIDATION=1 as an environment variable.
bfs_sssp binary would skip BFS part if SKIP_BFS=1 is present in your environment.

If you want to store/read generated graph from/to file use environment variables TMPFILE=<filename> and also REUSEFILE=1 to keep the file.
It's advised to use bfs_sssp binary to generate graph files as it generates both files of edges and weights (filename.weights)
bfs binary would only use/write edges file. And once bfs_sssp cant open weights file it would generate both files even if edges files is present.

N.B:

Current settings assume you are using powers of 2: total number of cores and number of cores per node.
It's possible to have non-power of two of nodes if you comment macro defined in common.h SIZE_MUST_BE_POWER_OF_TWO.
Be aware normally that will drop your performance by more then 20%.

If you want to use non-power of two processes per node, you should add -DPROCS_PER_NODE_NOT_POWER_OF_TWO to CFLAGS in src/Makefile,
this one will enable SIZE_MUST_BE_POWER_OF_TWO automatically.

For GreenGraph500 runs define compile time macro ENERGYLOOP_BFS or ENERGYLOOP_SSSP to measure energy in long enough loop to measure average power.

