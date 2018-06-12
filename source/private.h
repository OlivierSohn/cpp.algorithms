/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "../include/public.h"

#ifdef _WIN32
#  include <process.h>
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <Windows.h>
#  include "dbghelp.h"
#  pragma comment(lib, "Dbghelp.lib")
#  include <strsafe.h>      // for StringCchPrintfW
#  include <cstdlib>
#  include <stdio.h>
#  include <stdlib.h>
#else
#  include <execinfo.h>
#endif

#if !defined (_MSC_VER)
#include <dirent.h>
#include <sys/types.h>
//#include <pwd.h>
#include <cxxabi.h>
#endif

#include <unistd.h>
#include <sys/stat.h>
