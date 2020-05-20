#define _WIN64
#include <stdint.h>
/* ==================== GNU C and possibly other UNIX compilers ===================== */
#if !(defined(WIN32) || defined(_WIN64)) || defined(__GNUC__)

	#if defined(__GNUC__) || defined(__linux__)
		#define VOLATILE __volatile__
		#define ASM __asm__
	#else
		/* if we're neither compiling with gcc or under linux, we can hope
		 * the following lines work, they probably won't */
		#define ASM asm
		#define VOLATILE
	#endif

	#define myInt64 unsigned long long
	#define INT32 unsigned int

/* ======================== WIN32 ======================= */
#else

	#include <bitset>
	#include <array>

	#include <intrin.h>
	#pragma intrinsic(__rdtsc, __cpuid)
	//#pragma intrinsic(__cpuid)


	#define myInt64 uint64_t
	#define INT32 uint32_t

#endif

/* This is the RDTSC timer.
 * RDTSC is an instruction on several Intel and compatible CPUs that Reads the
 * Time Stamp Counter. The Intel manuals contain more information.
 */


#define COUNTER_LO(a) ((a).int32.lo)
#define COUNTER_HI(a) ((a).int32.hi)
#define COUNTER_VAL(a) ((a).int64)

#define COUNTER(a) \
	((unsigned long long)COUNTER_VAL(a))

#define COUNTER_DIFF(a,b) \
	(COUNTER(a)-COUNTER(b))

/* ==================== GNU C and possibly other UNIX compilers ===================== */
#if !(defined(WIN32) || defined(_WIN64)) || defined(__GNUC__)

	typedef union
	{       myInt64 int64;
		    struct {INT32 lo, hi;} int32;
	} tsc_counter;

  #define RDTSC(cpu_c) \
	  ASM VOLATILE ("rdtsc" : "=a" ((cpu_c).int32.lo), "=d"((cpu_c).int32.hi))
	#define CPUID() \
		ASM VOLATILE ("cpuid" : : "a" (0) : "bx", "cx", "dx" )

/* ======================== WIN32 ======================= */

#elif defined(_WIN64)

	typedef union
	{
	myInt64 int64;
	struct { INT32 lo, hi; } int32;
	} tsc_counter;

#define RDTSC(cpu_c) \
	(cpu_c).int64 = __rdtsc()
#if 0
#define CPUID(void)               \
	{                          \ 
		int cpu_test[4];\
		int* cpu_ptr = cpu_test;\
		__cpuid(cpu_ptr, 0);   \
		//__cpuid(cpui_info.data(), 0);  \
	}
#endif


#else

	typedef union
	{       myInt64 int64;
			struct {INT32 lo, hi;} int32;
	} tsc_counter;

	#define RDTSC(cpu_c)   \
	{       __asm{rdtsc}    \
			__asm{mov (cpu_c).int32.lo,eax}  \
			__asm{mov (cpu_c).int32.hi,edx}  \
	}

	#define CPUID() \
	{ \
		__asm mov eax, 0 \
		__asm cpuid \
	}

#endif


void init_tsc() {
	; // no need to initialize anything for x86
}

inline myInt64 start_tsc(void) {
    tsc_counter start;
	{                          
		int cpu_test[4];
		__cpuid(cpu_test, 0); 
		
	} 
    RDTSC(start);
    return COUNTER_VAL(start);
}

inline myInt64 stop_tsc(myInt64 start) {
	tsc_counter end;
	RDTSC(end);
	{
		int cpu_test[4];
		__cpuid(cpu_test, 0);

	}
	return COUNTER_VAL(end) - start;
}