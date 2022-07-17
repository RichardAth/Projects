#pragma once
/*  This file was created by Configure. */
#ifndef __PARICFG_H__
#define __PARICFG_H__
//#define UNIX
#define GPHELP "perl -S gphelp -detex"
#define GPDATADIR "/usr/local/share/pari"
#define SHELL_Q '\"'

#define PARIVERSION "GP/PARI CALCULATOR Version 2.13.3 (released)"
#define PARIINFO "amd64 running Windows, 64-bit version"
#define PARI_VERSION_CODE 134403
#define PARI_VERSION(a,b,c) (((a) << 16) + ((b) << 8) + (c))
#define PARI_VERSION_SHIFT 8
#define PARI_VCSVERSION ""
#define PARI_MT_ENGINE "single"

#define PARI_DOUBLE_FORMAT -
 //#define GCC_VERSION "gcc version 11.3.0 (GCC)"
#define ASMINLINE

/*  Location of GNU gzip program (enables reading of .Z and .gz files). */
#define GNUZCAT
#define ZCAT "/usr/bin/gzip -dc"


#undef UNIX
#define GNUZCAT
#undef ZCAT
#define ZCAT "gzip.exe -dc"

#define LONG_IS_64BIT
#define HAS_EXP2
#define HAS_LOG2
#define HAS_ISATTY
#define HAS_ALARM
#define HAS_SYSTEM
//#define USE_GETRUSAGE 1
//#define USE_GETTIMEOFDAY 1
//#define USE_TIMES 1
//#define HAS_SIGACTION
//#define HAS_WAITPID
#define HAS_GETENV
#define HAS_SETSID
//#define HAS_DLOPEN
#define STACK_CHECK
#define HAS_VSNPRINTF
//#define HAS_TIOCGWINSZ
#define HAS_STRFTIME
//#define HAS_STAT
//#define HAS_MMAP
#endif
