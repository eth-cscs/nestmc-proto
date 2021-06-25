#ifndef ARB_TYPES_H_
#define ARB_TYPES_H_

#ifdef __cplusplus
#include <cstdint>
using std::uint32_t;
#else
#include <stdint.h>
#endif

typedef double   arb_value_type;
typedef float    arb_weight_type;
typedef int      arb_index_type;
typedef uint32_t arb_size_type;

#endif // ARB_TYPES_H_
