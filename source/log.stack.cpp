

namespace imajuscule
{
#ifdef _WIN32 // dbg stack
    static void CreateMiniDump(EXCEPTION_POINTERS* pep)
    {
        // Open the file

        LG(ERR, "writing minidump...");
        HANDLE hFile = CreateFile(L"MiniDump.dmp", GENERIC_READ | GENERIC_WRITE,
            0, nullptr, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);

        if ( (hFile != nullptr) && (hFile != INVALID_HANDLE_VALUE) )
        {
            // Create the minidump

            MINIDUMP_EXCEPTION_INFORMATION mdei;

            mdei.ThreadId = GetCurrentThreadId();
            mdei.ExceptionPointers = pep;
            mdei.ClientPointers = FALSE;

            MINIDUMP_TYPE mdt = MiniDumpNormal;

            BOOL rv = MiniDumpWriteDump(GetCurrentProcess(), GetCurrentProcessId(),
                hFile, mdt, (pep != 0) ? &mdei : 0, 0, 0);

            if ( !rv )
                LG(ERR,"MiniDumpWriteDump failed. Error: %u", GetLastError());

            // Close the file

            CloseHandle(hFile);

        } else
        {
            LG(ERR, "CreateFile failed. Error: %u", GetLastError());
        }

    }
#endif

    void logStack()
{
#ifdef _WIN32 // dbg stack
    CreateMiniDump(0);
    return;
#else

    const char * lineS =  "    ----------------------------------------------------------------------------------------------";
    const char * Header = "                                           STACK BEGIN";
    const char * Footer = "                                           STACK END";
    const char * pipesS = "------------------------------------------------------------------------------------------------------";
  const char* endl = "\n";

    std::cout << endl
    << lineS << endl
    << Header << endl
    << pipesS << endl
    << endl;

    constexpr auto n_remove = 1; // for this function
    for(auto & t : debugging::getProgramStack(n_remove)) {
        std::cout << t << endl;
    }

    std::cout << endl
    << pipesS << endl
    << Footer << endl
    << lineS << endl
  << std::endl; // std::endl at the very end because it flushes.

#endif

	}
}
