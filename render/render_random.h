/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 - 2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __SURFACE_RENDER_RANDOM_H__
#define __SURFACE_RENDER_RANDOM_H__

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { unsigned long long state;  unsigned long long inc; } pcg32_random_t;

inline unsigned int pcg32_random_r(pcg32_random_t* rng)
{
    unsigned long long oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    unsigned int xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    unsigned int rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

/// Draw single precision pseudo random number in [0, 1]
inline float random_float(pcg32_random_t* rng)
{
  // If you want double precision random number in [0, 1],
  // http://mumble.net/~campbell/tmp/random_real.c

  // divide by 2^32
  return (float)((double)pcg32_random_r(rng) / (double)4294967296.0);
}


#endif // __SURFACE_RENDER_RANDOM_H__
