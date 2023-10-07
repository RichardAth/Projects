/* code 'borrowed' from YAFU */
/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may
benefit from your work.

                       --bbuhrow@gmail.com 7/28/10
----------------------------------------------------------------------*/

#include "pch.h"


int CLSIZE;
char HAS_SSE41;
#if defined(WIN32)
char sysname[MAX_COMPUTERNAME_LENGTH + 1];
unsigned long sysname_sz;
#else
char sysname[256];
int sysname_sz;
#endif

#ifdef _MSC_VER

/* Core aware timing on Windows, courtesy of Brian Gladman */

#if defined( _WIN64 )

#define current_processor_number GetCurrentProcessorNumber

#else

unsigned long current_processor_number(void)
{
    __asm
    {
        mov     eax, 1
        cpuid
        shr     ebx, 24
        mov     eax, ebx
    }
}

#endif

static int lock_thread_to_core(void)
{
    DWORD_PTR afp, afs;

    if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
    {
        afp &= (DWORD_PTR)(1LL << current_processor_number());
        if (SetThreadAffinityMask(GetCurrentThread(), afp))
            return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}

static int unlock_thread_from_core(void)
{
    DWORD_PTR afp, afs;

    if (GetProcessAffinityMask(GetCurrentProcess(), &afp, &afs))
    {
        if (SetThreadAffinityMask(GetCurrentThread(), afp))
            return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}

double cycles_per_second = 0.0;
double ticks_per_second = 0.0;
double cycles_per_tick = 0.0;

static unsigned long long measure_processor_speed(void)
{
    unsigned long long cycles;

    lock_thread_to_core();
    cycles = __rdtsc();
    Sleep(100);
    cycles = __rdtsc() - cycles;
    unlock_thread_from_core();
    cycles_per_second = 10.0 * (double)cycles;

    if (ticks_per_second == 0.0)
    {
        LARGE_INTEGER ll;
        QueryPerformanceFrequency(&ll);
        ticks_per_second = (double)ll.QuadPart;
        cycles_per_tick = cycles_per_second / ticks_per_second;
    }
    return cycles;
}


#else

uint64 measure_processor_speed(void)
{
    uint64 cycles;
    struct timeval start, stop;
    double t_time;
    TIME_DIFF* difference;

    gettimeofday(&start, NULL);

    cycles = yafu_read_clock();
    do
    {
        gettimeofday(&stop, NULL);
        difference = my_difftime(&start, &stop);
        t_time = ((double)difference->secs +
            (double)difference->usecs / 1000000);
        std::free(difference);
    } while (t_time < 0.1);
    cycles = yafu_read_clock() - cycles;

    return cycles;                  /* return cycles per second  */
}

#endif

#if defined(GCC_ASM32X)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
        ASM_G volatile(					\
            "movl %%ebx, %%esi   \n\t"		\
            "cpuid               \n\t"		\
            "movl %%ebx, %1      \n\t"		\
            "movl %%esi, %%ebx   \n\t"		\
            :"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
            :"0"(code) : "%esi")
#define CPUID2(code1, code2, a, b, c, d) 			\
        ASM_G volatile(					\
            "movl %%ebx, %%esi   \n\t"		\
            "cpuid               \n\t"		\
            "movl %%ebx, %1      \n\t"		\
            "movl %%esi, %%ebx   \n\t"		\
            :"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
            :"0"(code1), "2"(code2) : "%esi")

#elif defined(GCC_ASM64X)
#define HAS_CPUID
#define CPUID(code, a, b, c, d) 			\
        ASM_G volatile(					\
            "movq %%rbx, %%rsi   \n\t"		\
            "cpuid               \n\t"		\
            "movl %%ebx, %1      \n\t"		\
            "movq %%rsi, %%rbx   \n\t"		\
            :"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
            :"0"(code) : "%rsi")
#define CPUID2(code1, code2, a, b, c, d)		\
        ASM_G volatile(					\
            "movq %%rbx, %%rsi   \n\t"		\
            "cpuid               \n\t"		\
            "movl %%ebx, %1      \n\t"		\
            "movq %%rsi, %%rbx   \n\t"		\
            :"=a"(a), "=m"(b), "=c"(c), "=d"(d) 	\
            :"0"(code1), "2"(code2) : "%rsi")

#elif defined(_MSC_VER)
#include <intrin.h>
#define HAS_CPUID
#define CPUID(code, a, b, c, d)	\
    {	int _z[4]; \
        __cpuid(_z, code); \
        a = _z[0]; \
        b = _z[1]; \
        c = _z[2]; \
        d = _z[3]; \
    }
#define CPUID2(code1, code2, a, b, c, d) \
    {	int _z[4]; \
        __cpuidex(_z, code1, code2); \
        a = _z[0]; \
        b = _z[1]; \
        c = _z[2]; \
        d = _z[3]; \
    }
#endif


// http://msdn.microsoft.com/en-us/library/hskdteyh.aspx
// cpuid.cpp 
// processor: x86, x64
// Use the __cpuid intrinsic to get information about a CPU
// modified for c compliers and use of CPUID macros 
//		- brb, 10/26/10

const char* szFeatures[] =
{
    "x87 FPU On Chip",
    "Virtual-8086 Mode Enhancement",
    "Debugging Extensions",
    "Page Size Extensions",
    "Time Stamp Counter",
    "RDMSR and WRMSR Support",
    "Physical Address Extensions",
    "Machine Check Exception",
    "CMPXCHG8B Instruction",
    "APIC On Chip",
    "Unknown1",
    "SYSENTER and SYSEXIT",
    "Memory Type Range Registers",
    "PTE Global Bit",
    "Machine Check Architecture",
    "Conditional Move/Compare Instruction",
    "Page Attribute Table",
    "36-bit Page Size Extension",
    "Processor Serial Number",
    "CFLUSH Extension",
    "Unknown2",
    "Debug Store",
    "Thermal Monitor and Clock Ctrl",
    "MMX Technology",
    "FXSAVE/FXRSTOR",
    "SSE Extensions",
    "SSE2 Extensions",
    "Self Snoop",
    "Multithreading Technology",
    "Thermal Monitor",
    "Unknown4",
    "Pending Break Enable"
};



/* return cpuid , cache size, etc  */
static int extended_cpuid(char* CPUidstr, int* cachelinesize, char* bSSE41Extensions,
    int do_print)
{
    char CPUString[0x20];
    char CPUBrandString[0x40] = { 0 };
    int CPUInfo[4] = { -1 };
    int nSteppingID = 0;
    int nModel = 0;
    int nFamily = 0;
    int nProcessorType = 0;
    int nExtendedmodel = 0;
    int nExtendedfamily = 0;
    int nBrandIndex = 0;
    int nCLFLUSHcachelinesize = 0;
    int nLogicalProcessors = 0;
    int nAPICPhysicalID = 0;
    int nFeatureInfo = 0;
    int nCacheLineSize = 0;
    int nL2Associativity = 0;
    int nCacheSizeK = 0;
    int nPhysicalAddress = 0;
    int nVirtualAddress = 0;
    int nRet = 0;    /* return value */

    int nCores = 0;
    int nCacheType = 0;
    int nCacheLevel = 0;
    int nMaxThread = 0;
    int nSysLineSize = 0;
    int nPhysicalLinePartitions = 0;
    int nWaysAssociativity = 0;
    int nNumberSets = 0;

    unsigned    nIds, nExIds, i;

    char    bSSE3Instructions = 0;
    char    bMONITOR_MWAIT = 0;
    char    bCPLQualifiedDebugStore = 0;
    char    bVirtualMachineExtensions = 0;
    char    bEnhancedIntelSpeedStepTechnology = 0;
    char    bThermalMonitor2 = 0;
    char    bSupplementalSSE3 = 0;
    char    bL1ContextID = 0;
    char    bCMPXCHG16B = 0;
    char    bxTPRUpdateControl = 0;
    char    bPerfDebugCapabilityMSR = 0;
    //char    bSSE41Extensions = 0;
    char    bSSE42Extensions = 0;
    char    bPOPCNT = 0;

    char    bMultithreading = 0;

    char    bLAHF_SAHFAvailable = 0;
    char    bCmpLegacy = 0;
    char    bSVM = 0;
    char    bExtApicSpace = 0;
    char    bAltMovCr8 = 0;
    char    bLZCNT = 0;
    char    bSSE4A = 0;
    char    bMisalignedSSE = 0;
    char    bPREFETCH = 0;
    char    bSKINITandDEV = 0;
    char    bSYSCALL_SYSRETAvailable = 0;
    char    bExecuteDisableBitAvailable = 0;
    char    bMMXExtensions = 0;
    char    bFFXSR = 0;
    char    b1GBSupport = 0;
    char    bRDTSCP = 0;
    char    b64Available = 0;
    char    b3DNowExt = 0;
    char    b3DNow = 0;
    char    bNestedPaging = 0;
    char    bLBRVisualization = 0;
    char    bFP128 = 0;
    char    bMOVOptimization = 0;

    char    bSelfInit = 0;
    char    bFullyAssociative = 0;

    *bSSE41Extensions = 0;

    // __cpuid with an InfoType argument of 0 returns the number of
    // valid Ids in CPUInfo[0] and the CPU identification string in
    // the other three array elements. The CPU identification string is
    // not in linear order. The code below arranges the information 
    // in a human readable form.
    CPUID(0, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    //__cpuid(CPUInfo, 0);
    nIds = CPUInfo[0];
    std::memset(CPUString, 0, sizeof(CPUString));
    *((int*)CPUString) = CPUInfo[1];
    *((int*)(CPUString + 4)) = CPUInfo[3];
    *((int*)(CPUString + 8)) = CPUInfo[2];

    // Get the information associated with each valid Id
    for (i = 0; i <= nIds; ++i)
    {
        CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuid(CPUInfo, i);

        if (do_print > 1)
        {
            printf("\nFor InfoType %d\n", i);
            printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
            printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
            printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
            printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
        }

        // Interpret CPU feature information.
        if (i == 1)
        {
            nSteppingID = CPUInfo[0] & 0xf;
            nModel = (CPUInfo[0] >> 4) & 0xf;
            nFamily = (CPUInfo[0] >> 8) & 0xf;
            nProcessorType = (CPUInfo[0] >> 12) & 0x3;
            nExtendedmodel = (CPUInfo[0] >> 16) & 0xf;
            nExtendedfamily = (CPUInfo[0] >> 20) & 0xff;
            nBrandIndex = CPUInfo[1] & 0xff;
            *cachelinesize = nCLFLUSHcachelinesize = ((CPUInfo[1] >> 8) & 0xff) * 8;
            nLogicalProcessors = ((CPUInfo[1] >> 16) & 0xff);
            nAPICPhysicalID = (CPUInfo[1] >> 24) & 0xff;
            bSSE3Instructions = (CPUInfo[2] & 0x1) || 0;
            bMONITOR_MWAIT = (CPUInfo[2] & 0x8) || 0;
            bCPLQualifiedDebugStore = (CPUInfo[2] & 0x10) || 0;
            bVirtualMachineExtensions = (CPUInfo[2] & 0x20) || 0;
            bEnhancedIntelSpeedStepTechnology = (CPUInfo[2] & 0x80) || 0;
            bThermalMonitor2 = (CPUInfo[2] & 0x100) || 0;
            bSupplementalSSE3 = (CPUInfo[2] & 0x200) || 0;
            bL1ContextID = (CPUInfo[2] & 0x300) || 0;
            bCMPXCHG16B = (CPUInfo[2] & 0x2000) || 0;
            bxTPRUpdateControl = (CPUInfo[2] & 0x4000) || 0;
            bPerfDebugCapabilityMSR = (CPUInfo[2] & 0x8000) || 0;
            *bSSE41Extensions = (CPUInfo[2] & 0x80000) || 0;
            bSSE42Extensions = (CPUInfo[2] & 0x100000) || 0;
            bPOPCNT = (CPUInfo[2] & 0x800000) || 0;
            nFeatureInfo = CPUInfo[3];
            bMultithreading = (nFeatureInfo & (1 << 28)) || 0;
        }
    }

    // Calling __cpuid with 0x80000000 as the InfoType argument
    // gets the number of valid extended IDs.
    CPUID(0x80000000, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
    //__cpuid(CPUInfo, 0x80000000);
    nExIds = CPUInfo[0];
    std::memset(CPUBrandString, 0, sizeof(CPUBrandString));

    // Get the information associated with each extended ID.
    for (i = 0x80000000; i <= nExIds; ++i)
    {
        CPUID(i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuid(CPUInfo, i);
        if (do_print > 1)
        {
            printf("\nFor InfoType %x\n", i);
            printf("CPUInfo[0] = 0x%x\n", CPUInfo[0]);
            printf("CPUInfo[1] = 0x%x\n", CPUInfo[1]);
            printf("CPUInfo[2] = 0x%x\n", CPUInfo[2]);
            printf("CPUInfo[3] = 0x%x\n", CPUInfo[3]);
        }

        if (i == 0x80000001)
        {
            bLAHF_SAHFAvailable = (CPUInfo[2] & 0x1) || 0;
            bCmpLegacy = (CPUInfo[2] & 0x2) || 0;
            bSVM = (CPUInfo[2] & 0x4) || 0;
            bExtApicSpace = (CPUInfo[2] & 0x8) || 0;
            bAltMovCr8 = (CPUInfo[2] & 0x10) || 0;
            bLZCNT = (CPUInfo[2] & 0x20) || 0;
            bSSE4A = (CPUInfo[2] & 0x40) || 0;
            bMisalignedSSE = (CPUInfo[2] & 0x80) || 0;
            bPREFETCH = (CPUInfo[2] & 0x100) || 0;
            bSKINITandDEV = (CPUInfo[2] & 0x1000) || 0;
            bSYSCALL_SYSRETAvailable = (CPUInfo[3] & 0x800) || 0;
            bExecuteDisableBitAvailable = (CPUInfo[3] & 0x10000) || 0;
            bMMXExtensions = (CPUInfo[3] & 0x40000) || 0;
            bFFXSR = (CPUInfo[3] & 0x200000) || 0;
            b1GBSupport = (CPUInfo[3] & 0x400000) || 0;
            bRDTSCP = (CPUInfo[3] & 0x8000000) || 0;
            b64Available = (CPUInfo[3] & 0x20000000) || 0;
            b3DNowExt = (CPUInfo[3] & 0x40000000) || 0;
            b3DNow = (CPUInfo[3] & 0x80000000) || 0;
        }

        // Interpret CPU brand string and cache information.
        if (i == 0x80000002)
            std::memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            std::memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            std::memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000006)
        {
            nCacheLineSize = CPUInfo[2] & 0xff;
            nL2Associativity = (CPUInfo[2] >> 12) & 0xf;
            nCacheSizeK = (CPUInfo[2] >> 16) & 0xffff;
        }
        else if (i == 0x80000008)
        {
            nPhysicalAddress = CPUInfo[0] & 0xff;
            nVirtualAddress = (CPUInfo[0] >> 8) & 0xff;
        }
        else if (i == 0x8000000A)
        {
            bNestedPaging = (CPUInfo[3] & 0x1) || 0;
            bLBRVisualization = (CPUInfo[3] & 0x2) || 0;
        }
        else if (i == 0x8000001A)
        {
            bFP128 = (CPUInfo[0] & 0x1) || 0;
            bMOVOptimization = (CPUInfo[0] & 0x2) || 0;
        }
    }
    ptrdiff_t ix = 0;
    while (ix < 64 &&  std::isblank(CPUBrandString[ix]))
        ix++;        /* find 1st non-blank character */
    strcpy_s(CPUidstr, 64 - ix, CPUBrandString + ix);
    // Display all the information in user-friendly format.
    if (do_print > 0)
        printf("\n\nCPU String: %s\n", CPUString);

    if (nIds >= 1)
    {
        if (do_print > 0)
        {
            if (nSteppingID)
                printf("Stepping ID = %d\n", nSteppingID);
            if (nModel)
                printf("Model = %d\n", nModel);
            if (nFamily)
                printf("Family = %d\n", nFamily);
            if (nProcessorType)
                printf("Processor Type = %d\n", nProcessorType);
            if (nExtendedmodel)
                printf("Extended model = %d\n", nExtendedmodel);
            if (nExtendedfamily)
                printf("Extended family = %d\n", nExtendedfamily);
            if (nBrandIndex)
                printf("Brand Index = %d\n", nBrandIndex);
            if (nCLFLUSHcachelinesize)
                printf("CLFLUSH cache line size = %d\n",
                    nCLFLUSHcachelinesize);
            if (bMultithreading && (nLogicalProcessors > 0))
                printf("Logical Processor Count = %d\n", nLogicalProcessors);
            if (nAPICPhysicalID)
                printf("APIC Physical ID = %d\n", nAPICPhysicalID);

            if (nFeatureInfo || bSSE3Instructions ||
                bMONITOR_MWAIT || bCPLQualifiedDebugStore ||
                bVirtualMachineExtensions || bEnhancedIntelSpeedStepTechnology ||
                bThermalMonitor2 || bSupplementalSSE3 || bL1ContextID ||
                bCMPXCHG16B || bxTPRUpdateControl || bPerfDebugCapabilityMSR ||
                *bSSE41Extensions || bSSE42Extensions || bPOPCNT ||
                bLAHF_SAHFAvailable || bCmpLegacy || bSVM ||
                bExtApicSpace || bAltMovCr8 ||
                bLZCNT || bSSE4A || bMisalignedSSE ||
                bPREFETCH || bSKINITandDEV || bSYSCALL_SYSRETAvailable ||
                bExecuteDisableBitAvailable || bMMXExtensions || bFFXSR || b1GBSupport ||
                bRDTSCP || b64Available || b3DNowExt || b3DNow || bNestedPaging ||
                bLBRVisualization || bFP128 || bMOVOptimization)
            {
                printf("\nThe following features are supported:\n");

                if (bSSE3Instructions)
                    printf("\tSSE3\n");
                if (bMONITOR_MWAIT)
                    printf("\tMONITOR/MWAIT\n");
                if (bCPLQualifiedDebugStore)
                    printf("\tCPL Qualified Debug Store\n");
                if (bVirtualMachineExtensions)
                    printf("\tVirtual Machine Extensions\n");
                if (bEnhancedIntelSpeedStepTechnology)
                    printf("\tEnhanced Intel SpeedStep Technology\n");
                if (bThermalMonitor2)
                    printf("\tThermal Monitor 2\n");
                if (bSupplementalSSE3)
                    printf("\tSupplemental Streaming SIMD Extensions 3\n");
                if (bL1ContextID)
                    printf("\tL1 Context ID\n");
                if (bCMPXCHG16B)
                    printf("\tCMPXCHG16B Instruction\n");
                if (bxTPRUpdateControl)
                    printf("\txTPR Update Control\n");
                if (bPerfDebugCapabilityMSR)
                    printf("\tPerf\\Debug Capability MSR\n");
                if (*bSSE41Extensions)
                    printf("\tSSE4.1 Extensions\n");
                if (bSSE42Extensions)
                    printf("\tSSE4.2 Extensions\n");
                if (bPOPCNT)
                    printf("\tPPOPCNT Instruction\n");

                i = 0;
                nIds = 1;
                while (i < (sizeof(szFeatures) / sizeof(const char*)))
                {
                    if (nFeatureInfo & nIds)
                    {
                        printf("\t");
                        printf("%s", szFeatures[i]);
                        printf("\n");
                    }

                    nIds <<= 1;
                    ++i;
                }
                if (bLAHF_SAHFAvailable)
                    printf("\tLAHF/SAHF in 64-bit mode\n");
                if (bCmpLegacy)
                    printf("\tCore multi-processing legacy mode\n");
                if (bSVM)
                    printf("\tSecure Virtual Machine\n");
                if (bExtApicSpace)
                    printf("\tExtended APIC Register Space\n");
                if (bAltMovCr8)
                    printf("\tAltMovCr8\n");
                if (bLZCNT)
                    printf("\tLZCNT instruction\n");
                if (bSSE4A)
                    printf("\tSSE4A (EXTRQ, INSERTQ, MOVNTSD, MOVNTSS)\n");
                if (bMisalignedSSE)
                    printf("\tMisaligned SSE mode\n");
                if (bPREFETCH)
                    printf("\tPREFETCH and PREFETCHW Instructions\n");
                if (bSKINITandDEV)
                    printf("\tSKINIT and DEV support\n");
                if (bSYSCALL_SYSRETAvailable)
                    printf("\tSYSCALL/SYSRET in 64-bit mode\n");
                if (bExecuteDisableBitAvailable)
                    printf("\tExecute Disable Bit\n");
                if (bMMXExtensions)
                    printf("\tExtensions to MMX Instructions\n");
                if (bFFXSR)
                    printf("\tFFXSR\n");
                if (b1GBSupport)
                    printf("\t1GB page support\n");
                if (bRDTSCP)
                    printf("\tRDTSCP instruction\n");
                if (b64Available)
                    printf("\t64 bit Technology\n");
                if (b3DNowExt)
                    printf("\t3Dnow Ext\n");
                if (b3DNow)
                    printf("\t3Dnow! instructions\n");
                if (bNestedPaging)
                    printf("\tNested Paging\n");
                if (bLBRVisualization)
                    printf("\tLBR Visualization\n");
                if (bFP128)
                    printf("\tFP128 optimization\n");
                if (bMOVOptimization)
                    printf("\tMOVU Optimization\n");
            }
        }
    }

    if (nExIds >= 0x80000004 && do_print > 0)
        printf("\nCPU Brand String: %s\n", CPUBrandString);

    if (nExIds >= 0x80000006 && do_print > 0)
    {
        printf("Cache Line Size = %d\n", nCacheLineSize);
        printf("L2 Associativity = %d\n", nL2Associativity);
        printf("Cache Size = %dK\n", nCacheSizeK);
    }


    for (i = 0;; i++)
    {
        CPUID2(0x4, i, CPUInfo[0], CPUInfo[1], CPUInfo[2], CPUInfo[3]);
        //__cpuidex(CPUInfo, 0x4, i);
        if (!(CPUInfo[0] & 0xf0))
            break;

        if (i == 0) {
            nCores = CPUInfo[0] >> 26;
            if (do_print > 0)
                printf("\n\nNumber of Cores = %d\n", nCores + 1);
        }

        nCacheType = (CPUInfo[0] & 0x1f);
        nCacheLevel = (CPUInfo[0] & 0xe0) >> 5;
        bSelfInit = (CPUInfo[0] & 0x100) >> 8;
        bFullyAssociative = (CPUInfo[0] & 0x200) >> 9;
        nMaxThread = (CPUInfo[0] & 0x03ffc000) >> 14;
        nSysLineSize = (CPUInfo[1] & 0x0fff);
        nPhysicalLinePartitions = (CPUInfo[1] & 0x03ff000) >> 12;
        nWaysAssociativity = (CPUInfo[1]) >> 22;
        nNumberSets = CPUInfo[2];

        if (do_print > 0) {
            printf("\n");

            printf("ECX Index %d\n", i);
            switch (nCacheType)
            {
            case 0:
                printf("   Type: Null\n");
                break;
            case 1:
                printf("   Type: Data Cache\n");
                break;
            case 2:
                printf("   Type: Instruction Cache\n");
                break;
            case 3:
                printf("   Type: Unified Cache\n");
                break;
            default:
                printf("   Type: Unknown\n");
            }

            printf("   Level = %d\n", nCacheLevel + 1);
            if (bSelfInit) {
                printf("   Self Initializing\n");
            }
            else {
                printf("   Not Self Initializing\n");
            }
            if (bFullyAssociative) {
                printf("   Is Fully Associatve\n");
            }
            else {
                printf("   Is Not Fully Associatve\n");
            }
            printf("   Max Threads = %d\n", nMaxThread + 1);
            printf("   System Line Size = %d\n", nSysLineSize + 1);
            printf("   Physical Line Partions = %d\n", nPhysicalLinePartitions + 1);
            printf("   Ways of Associativity = %d\n", nWaysAssociativity + 1);
            printf("   Number of Sets = %d\n", nNumberSets + 1);
        }
    }

    return  nRet;
}

// function containing system commands to get the computer name, CPU speed, etc
void get_computer_info(char* CPUidstr, double &MEAS_CPU_FREQUENCY)
{

    //figure out cpu freq. 0.1 seconds won't be very accurate 
    MEAS_CPU_FREQUENCY = measure_processor_speed() / 1.0e5;

#ifdef __APPLE__
    // something in extended cpuid causes a segfault on mac builds.
    // just disable it for now - this information is not critical for
    // program operation.
    std::strcpy(idstr, "N/A");
    CLSIZE = 0;
    L1CACHE = DEFAULT_L1_CACHE_SIZE;
    L2CACHE = DEFAULT_L2_CACHE_SIZE;
    HAS_SSE41 = 0;

#else

    // run an extended cpuid command to get the cache line size, and
    // optionally print a bunch of info to the screen

    extended_cpuid(CPUidstr, &CLSIZE, &HAS_SSE41, verbose);
#endif

#if defined(WIN32)

    sysname_sz = MAX_COMPUTERNAME_LENGTH + 1;
    GetComputerNameA(sysname, &sysname_sz);

#else

    ret = gethostname(sysname, sizeof(sysname) / sizeof(*sysname));
    sysname[(sizeof(sysname) - 1) / sizeof(*sysname)] = 0;	// null terminate
    if (ret != 0)
    {
        printf("error occured when getting host name\n");
        std::strcpy(sysname, "N/A");
    }
    sysname_sz = strlen(sysname);

#endif

    return;
}
