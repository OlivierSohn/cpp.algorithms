#ifdef __APPLE__
#include "TargetConditionals.h"
#endif

#if defined ( WIN32 )
#define __func__ __FUNCTION__
#endif

namespace imajuscule
{
#define ERR_LOG(x,type) do{ \
LG(ERR, "%s : ' %s '\n   |in %s\n   |(%s,%d)", #type, #x, __func__, __FILE__, __LINE__ ); \
logStack(); \
} while(0)

// "soft" Assert
#define C(x) if(!(x)) {ERR_LOG(x,check);} else do{}while(0)

/*
"A" is used to replace "assert(x)" with "Assert(x)"

with the advantage that :
- a comprehensive error message is logged, and the stack is available in debugger (not always the case with assert)
- in release, code execution is not broken but an error message is logged
*/

// use __builtin_trap to stop execution and allow debugging
#ifndef NDEBUG
#define ASSERT__THROW do{ __builtin_trap(); throw std::logic_error("assertion failed");}while(0)
#else
#define ASSERT__THROW do{}while(0)
#endif

#define ASSERT_ERR_LOG(x) ERR_LOG(x,assert)
#define ASSERT_ERR(x) do{ASSERT_ERR_LOG(x); ASSERT__THROW;}while(0)

#ifndef NDEBUG
#define if_A(x) if(unlikely(!(x))) ASSERT_ERR(x); else
#else
#define if_A(x)
#endif

#ifndef NDEBUG
#define Assert(x) if_A(x) do {} while ( 0 )
#else
#define Assert(x) do {} while ( 0 )
#endif

}
