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
#  include <cstdlib>
#  include <stdio.h>
#  include <stdlib.h>
#  include <unistd.h>
#else
#  include <dirent.h>
#  include <sys/types.h>
#  include <unistd.h>
#  include <pwd.h>
#  include <execinfo.h>
#  include <cxxabi.h>
#endif
