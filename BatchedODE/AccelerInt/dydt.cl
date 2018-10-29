#include "species_rates.oclh"
void dydt (__private __ValueType const t,
           __global __ValueType const * __restrict__ param,
           __global __ValueType const * __restrict__ y,
           __global __ValueType * __restrict__ dy,
           __global __ValueType * __restrict__ rwk)
{
    species_rates(0, param, y, dy, rwk);
}
