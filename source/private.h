/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "../include/public.h"

#define IMJ_WITH_DEBUG_STACK
#ifdef _WIN32 // dbg stack
#  include <process.h>
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <Windows.h>
#  include "dbghelp.h"
#  pragma comment(lib, "Dbghelp.lib")
#  include <strsafe.h>      // for StringCchPrintfW
#  include <cstdlib>
#  include <stdlib.h>
#elif __has_include(<execinfo.h>)
#  include <execinfo.h>
#else
# undef IMJ_WITH_DEBUG_STACK
#endif

#if !defined (_MSC_VER)
#  if defined _WIN32 // mingw
#    include <direct.h>
#  endif

#  include <dirent.h>
#  include <sys/types.h>
#  include <cxxabi.h>
#endif

#include <unistd.h>
#include <sys/stat.h>
