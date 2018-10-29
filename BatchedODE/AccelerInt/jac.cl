#include "jacobian.oclh"
void jacob (__private __ValueType const t,
            __global __ValueType const * __restrict__ param,
            __global __ValueType const * __restrict__ y,
            __global __ValueType * __restrict__ jac,
            __global __ValueType * __restrict__ rwk)
{
   jacobian(0, param, y, jac, rwk);
}
