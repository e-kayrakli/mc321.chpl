### Files

- `mc321.c`: Base implementation in C.
- `mc321.chpl`: The experimental Chapel version that can run on GPUs.
- `mc321.out`: 10M photon output of the base implementation
- `mc321.out.chpl`: 10M photon output of the first Chapel port that's not
  parallel
- `mc321.out.chpl.gpu`: 10M photon output of the GPU-enabled Chapel version

### How to compile

```
> gcc -lm -O3 mc321.c   # for the C version
> chpl --fast mc321.chpl -suseCudaRng=true  # for the Chapel version
```

### Notes

- The Chapel version wants to use the custom RNG by default. However, there are
  some issues launching the kernel with that. Instead, `-suseCudaRng=true` can
  be used to use CUDA's RNG, which works fine.

### Performance

Below are some results. The Chapel version uses 1M GPU threads. All runs are
using 1B photons.

RTX A2000 is a low-end GPU, whereas A100 is an HPC-grade GPU. There are more
powerful ones that I can test with as well.

|                         |  Sequential C  | Chapel w/ RTX A2000 | Chapel w/ A100
|-------------------------|----------------|---------------------|---------------
|Time (s)                 |  83.57         | 6.67                | 0.80
|Throughput (MPhotons/s)  |  11.97         | 149.94              | 1254.21
|Relative Performance     |  1.0x          | 12.53x              | 104.78


- Chapel 2.2 (compiled with only `--fast`)
- gcc 13.2 (compiled with `-O3 -lm`)

#### Pushing the limits

On A100, I wanted to see how high the throughput can get and played around a bit
with number of photons and threads. With 100B photons and 100M threads, I get
the following, which was the highest throughput I've seen so far:

```
Number of photons : 100000000000
MPhotons/s : 1384.64
Elapsed time : 72.2208
```

This is about 115x better throughput than the sequential C version.


