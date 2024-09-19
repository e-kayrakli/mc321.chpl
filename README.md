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

Just using `time`, I see that the Chapel version runs in 32 seconds compared to
C's 69 seconds using 1B photons. The GPU-enabled Chapel version is right now the
first working prototype and can benefit from several optimizations.


- CPU: i5-11400 (6 cores)
- GPU: RTX A2000
- Chapel 2.2 (compiled with only `--fast`)
- gcc 13.2 (compiled with `-O3 -lm`)
