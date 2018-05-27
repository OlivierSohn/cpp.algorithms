

namespace imajuscule
{
#ifdef _WIN32
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
#ifdef _WIN32
    CreateMiniDump(0);
    return;
#else
    
    const char * lineS =  "    ----------------------------------------------------------------------------------------------";
    const char * Header = "                                           STACK BEGIN";
    const char * Footer = "                                           STACK END";
    const char * pipesS = "------------------------------------------------------------------------------------------------------";

    std::cout << std::endl
    << lineS << std::endl
    << Header << std::endl
    << pipesS << std::endl
    << std::endl;
    
    constexpr auto n_remove = 1; // for this function
    for(auto & t : debugging::getProgramStack(n_remove)) {
        std::cout << t << std::endl;
    }
    
    std::cout << std::endl
    << pipesS << std::endl
    << Footer << std::endl
    << lineS << std::endl
    << std::endl;
    
#endif

	}
}
