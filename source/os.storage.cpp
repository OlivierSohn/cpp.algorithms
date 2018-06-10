
using namespace imajuscule;
using namespace StorageStuff;

namespace imajuscule {
    DirectoryPath DirectoryPath::referentiablesPath;
    DirectoryPath DirectoryPath::capturePath;
    
    bool split_path(std::string const & str, DirectoryPath & dir, FileName & filename)
    {
        using namespace std;
        auto pos = str.find_last_of("/");
        if(pos == string::npos) {
            return false;
        }
        dir = {std::string(str.begin(), str.begin()+pos)};
        filename = std::string(str.begin()+pos+1, str.end());
        return true;
    }
    
    std::string getParentDirectory(std::string path)
    {
        using namespace std;
        auto pos = path.find_last_of("/");
        if(pos == string::npos) {
            throw logic_error("no parent");
        }
        path.resize(pos);

        return std::move(path);
    }
}

std::set<std::string> WritableStorage::g_openedForWrite;

void DirectoryPath::setReferentiablesDir(DirectoryPath const & d) {
    referentiablesPath = d;
}
bool DirectoryPath::getReferentiablesDir(DirectoryPath & p) {
    p = referentiablesPath;
    return !referentiablesPath.empty();
}

void DirectoryPath::setCaptureImageDir(DirectoryPath const & d) {
    capturePath = d;
}
bool DirectoryPath::getCaptureImageDir(DirectoryPath & p) {
    p = capturePath;
    return !capturePath.empty();
}

ReadableStorage::ReadableStorage(DirectoryPath const &d, FileName const &f) :
m_pFile(nullptr),
m_bufferReadPos(0),
m_directoryPath(d),
m_filename(f)
{}


ReadableStorage::~ReadableStorage()
{
    CloseFile();
}

void ReadableStorage::CloseFile()
{
    if (m_pFile) {
        fclose((FILE*)m_pFile);
        m_pFile = 0;
    }
    m_bufferReadPos = 0;
}

void WritableStorage::Finalize()
{
    FlushData();
    CloseFile();
}


eResult ReadableStorage::OpenFileForOperation(const std::string & sFilePath, FileMode op)
{
    //LG(INFO, "WritableStorage::OpenFileForOperation( %s, %s) begin", sFilePath.c_str(), FileOperationToString(op));
    
    CloseFile();
    
    m_pFile = fopen(sFilePath.c_str(), op == FileMode::READ ? "rb" : "wb");
    
    if ( unlikely(!m_pFile)) {
        LG(ERR, "WritableStorage::OpenFileForOperation : fopen failed : %d", errno);
        return ILE_BAD_PARAMETER;
    }
    if (op == FileMode::READ) {
        m_bufferReadPos = 0;
        ReadToBuffer();
    }
    return ILE_SUCCESS;
}

eResult ReadableStorage::OpenForRead()
{
    std::string filePath;
    
    for( auto const & directory_name : m_directoryPath)
    {
        filePath.append(directory_name);
        filePath.append("/");
    }
    
    filePath.append(m_filename);
    ReplaceStringInPlace(filePath, "//", "/" );
    eResult ret = OpenFileForOperation(filePath, FileMode::READ);
    if ( unlikely(ret != ILE_SUCCESS))
    {
        LG(ERR, "WritableStorage::OpenForRead : OpenFileForOperation returned %d", ret);
    }
    
    return ret;
}
eResult WritableStorage::OpenForWrite()
{
    eResult ret = ILE_SUCCESS;
    for( auto const & directory_name : m_directoryPath)
    {
        m_filePath.append( directory_name );
        m_filePath.append("/");
        if (!dirExists(m_filePath))
        {
            ret = makeDir(m_filePath);
            if ( unlikely(ILE_SUCCESS != ret))
            {
                LG(ERR, "WritableStorage::OpenForWrite : WritableStorage::makeDir(%s) error : %d", m_filePath.c_str(), ret);
                return ret;
            }
        }
    }
    
    m_filePath.append(m_filename);
    
    auto it2 = g_openedForWrite.find(m_filePath);
    if(it2 != g_openedForWrite.end())
    {
        return ILE_RECURSIVITY;
    }
    g_openedForWrite.insert(m_filePath);
    
    ret = OpenFileForOperation(m_filePath, FileMode::WRITE);
    if ( unlikely(ret != ILE_SUCCESS))
    {
        LG(ERR, "WritableStorage::OpenForWrite : OpenFileForOperation returned %d", ret);
    }
    
    return ret;
}

eResult WritableStorage::Save()
{
    m_filePath.clear();
    
    eResult ret = doSaveBegin();
    if (ret != ILE_SUCCESS)
    {
        if( unlikely(ret != ILE_RECURSIVITY)) {
            LG(ERR, "WritableStorage::Save : doSaveBegin returned %d", ret);
        }
        return ret;
    }
    
    ret = doSave();
    if (unlikely(ret != ILE_SUCCESS))
    {
        LG(ERR, "WritableStorage::Save : doSave returned %d", ret);
        return ret;
    }
    
    doSaveEnd();
    
    return ILE_SUCCESS;
}

eResult WritableStorage::doSaveBegin()
{
    {
        eResult ret = OpenForWrite();
        if ( ret != ILE_SUCCESS )
        {
            if ( unlikely(ret != ILE_RECURSIVITY) )
                LG(ERR, "WritableStorage::SaveBegin : OpenForWrite returned %d", ret);
            return ret;
        }
    }
    
    // to reserve the header space (it will be overwritten in ::SaveEnd())
    DoUpdateFileHeader();
    
    return ILE_SUCCESS;
}
eResult WritableStorage::doSave()
{
    return ILE_SUCCESS;
}
void WritableStorage::doSaveEnd()
{
    UpdateFileHeader();    
}

void WritableStorage::UpdateFileHeader()
{
    if(!m_pFile) {
        assert(0);
        return;
    }
    // write the data for this frame
    FlushMyBuffer();
    
    // now that the data has been written, we can modify the file position
    fpos_t curPos;
    if (!fgetpos((FILE*)m_pFile, &curPos))
    {
        rewind((FILE*)m_pFile);
        
        DoUpdateFileHeader();
        
        int ret = FlushData();
        if (unlikely(ret))
        {
            LG(ERR, "WritableStorage::UpdateFileHeader : FlushData returned %d", ret );
        }
        
        if (likely(!fsetpos((FILE*)m_pFile, &curPos)))
        {
            
        }
        else
        {
            LG(ERR, "WritableStorage::UpdateFileHeader : fsetpos failed : %d", errno);
            Assert(0);
        }
    }
    else
    {
        LG(ERR, "WritableStorage::UpdateFileHeader : fgetpos failed : %d", errno);
        Assert(0);
    }
}

void WritableStorage::FlushMyBuffer()
{
    size_t count = m_writeBuffer.size();
    if (count == 0)
        return;
#ifdef _WIN32
    _fwrite_nolock(m_writeBuffer.data(), 1, count, (FILE*)m_pFile);
#else
    fwrite(m_writeBuffer.data(), 1, count, (FILE*)m_pFile);
#endif
    
    m_writeBuffer.clear();
}

void ReadableStorage::ReadToBuffer()
{
#ifdef _WIN32
    auto result = _fread_nolock(m_freadBuffer, 1, SIZE_READ_BUFFER, (FILE*)m_pFile);
#else
    auto result = fread(m_freadBuffer, 1, SIZE_READ_BUFFER, (FILE*)m_pFile);
#endif
  // we don't check the result and it's okay
  //if (result != SIZE_READ_BUFFER) {fputs ("Reading error",stderr); exit (3);}
}

void ReadableStorage::ReadData(void * p, size_t size, size_t count)
{
    //LG(INFO, "WritableStorage::ReadData(%x, %d, %d)", p, size, count);
    
    size_t total = size * count;
    
    do
    {
        //LG(INFO, "WritableStorage::ReadData m_bufferReadPos = %d", m_bufferReadPos);
        
        size_t max = m_bufferReadPos + total;
        
        //LG(INFO, "WritableStorage::ReadData max = %d", max);
        
        long secondRead = max - SIZE_READ_BUFFER;
        
        //LG(INFO, "WritableStorage::ReadData secondRead = %d", secondRead);
        
        if (secondRead > 0)
        {
            //LG(INFO, "WritableStorage::ReadData secondRead > 0");
            
            long i = SIZE_READ_BUFFER - m_bufferReadPos;
            
            memcpy(p, &m_freadBuffer[m_bufferReadPos], i);
            
            m_bufferReadPos = 0;
            
            //LG(INFO, "WritableStorage::ReadData before ReadToBuffer");
            ReadToBuffer();
            //LG(INFO, "WritableStorage::ReadData after ReadToBuffer");
            
            total -= i;
            p = (char*)p + i;
        }
        else
        {
            //LG(INFO, "WritableStorage::ReadData secondRead < 0");
            
            memcpy(p, &m_freadBuffer[m_bufferReadPos], total);
            m_bufferReadPos = max;
            
            break;
        }
    }
    while(total>0);
    
    //LG(INFO, "WritableStorage::ReadData end");
}

int WritableStorage::WriteData(void const * p, size_t size, size_t count)
{
    size_t add = size*count;
    m_writeBuffer.insert(m_writeBuffer.end(), (unsigned char*)p, ((unsigned char*)p) + add);
    return add;
}

int WritableStorage::FlushData()
{
    FlushMyBuffer();
    
#ifdef _WIN32
    return _fflush_nolock((FILE*)m_pFile);
#else
    return fflush((FILE*)m_pFile);
#endif
}

namespace imajuscule {
    
    DirectoryPath DirectoryPath::root() {
        return Platform::user_path() + "grid3d";
    }
    
    namespace StorageStuff {
    const char * FileOperationToString(WritableStorage::FileMode op)
    {
        switch (op)
        {
            case WritableStorage::FileMode::WRITE:
                return "FileMode::WRITE";
                break;
            case WritableStorage::FileMode::READ:
                return "FileMode::READ";
                break;
            default:
                return "UNKNOWN";
                break;
        }
    }
    
#ifdef _WIN32
    void string_cast(const wchar_t* pSource, unsigned int codePage, std::string & oCast)
    {
        Assert(pSource != 0);
        oCast.clear();
        size_t sourceLength = std::wcslen(pSource);
        if (likely(sourceLength > 0))
        {
            int length = ::WideCharToMultiByte(codePage, 0, pSource, sourceLength, nullptr, 0, nullptr, nullptr);
            if (likely(length != 0))
            {
                StackVector<char> buffer(length);
                ::WideCharToMultiByte(codePage, 0, pSource, sourceLength, &buffer[0], length, nullptr, nullptr);
                oCast.assign(buffer.begin(), buffer.end());
            }
        }
    }
#endif
    
    bool dirExists(const std::string & path)
    {
        //LG(INFO, "dirExists(%s)", (path.c_str() ? path.c_str() : "nullptr"));
        bool bExists = false;
        
        struct stat info;
        bExists = ((stat(path.c_str(), &info) == 0) && (info.st_mode & S_IFDIR));
        
        //LG(INFO, "dirExists(%s) returns %s", (path.c_str() ? path.c_str() : "nullptr"), (bExists ? "true" : "false"));
        return bExists;
    }
    
    bool fileExists(const std::string & path)
    {
        //LG(INFO, "fileExists(%s)", (path.c_str() ? path.c_str() : "nullptr"));
        bool bExists = false;
        
        struct stat info;
        bExists = (stat(path.c_str(), &info) == 0) && !(info.st_mode & S_IFDIR);
        
        //LG(INFO, "fileExists(%s) returns %s", (path.c_str() ? path.c_str() : "nullptr"), (bExists ? "true" : "false"));
        return bExists;
    }
    
    bool fileCreationDate(const std::string & path, std::string & oDate)
    {
        //LG(INFO, "fileCreationDate(%s)", (path.c_str() ? path.c_str() : "nullptr"));
        
        oDate.clear();
        
        bool bExists = false;
        
        struct stat info;
        bExists = (stat(path.c_str(), &info) == 0) && !(info.st_mode & S_IFDIR);
        if (likely(bExists))
        {
            tm * time = gmtime((const time_t*)&(info.st_mtime));
            FormatDate(time, oDate);
        }
        else
        {
            LG(ERR, "fileCreationDate : file does not exist");
            oDate.assign("../../.. ..:..:..");
        }
        
        //LG(INFO, "fileCreationDate(%s) returns %s", (path.c_str() ? path.c_str() : "nullptr"), (bExists ? "true" : "false"));
        return bExists;
    }
    
        eResult removeDir(const std::string & path)
        {
            if(rmdir(path.c_str())) {
                LG(ERR, "rmdir error %d", errno);
                return ILE_ERROR;
            }
            return ILE_SUCCESS;
        }
        eResult removeFile(const std::string & path)
        {
            if(unlink(path.c_str())) {
                LG(ERR, "unlink error %d", errno);
                return ILE_ERROR;
            }
            return ILE_SUCCESS;
        }
        
    eResult makeDir(const std::string & path)
    {
        //LG(INFO, "makeDir(%s)", (path.c_str() ? path.c_str() : "nullptr"));
        eResult res = ILE_SUCCESS;
        
#ifdef _WIN32
        std::wstring swName = std::wstring(path.begin(), path.end());
        const wchar_t * pwStr = swName.c_str();
        //LG(INFO, "makeDir : before CreateDirectory");
        if (!CreateDirectory(pwStr, nullptr))
        {
            DWORD dwErr = GetLastError();
            if (unlikely(dwErr != ERROR_ALREADY_EXISTS))
            {
                LG(ERR, "makeDir : CreateDirectory error %x", dwErr);
                res = ILE_ERROR;
            }
            else
            {
                //LG(INFO, "makeDir : directory already exists");
            }
        }
#else
        
        int ret;
        ret = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (ret != 0)
        {
            if( unlikely(errno != EEXIST))
            {
                LG(ERR, "makeDir : CreateDirectory error %x", errno);
                res = ILE_ERROR;
            }
            else
            {
                //LG(INFO, "makeDir : directory already exists");
            }
        }
        
#endif
        
        //LG(INFO, "makeDir(%s) returns %d", (path.c_str() ? path.c_str() : "nullptr"), res);
        return res;
    }
    
        
        std::vector< std::string > listFilenames( const DirectoryPath & path ) {
            return listFilenames(path.toString());
        }
        std::vector< std::string > listFilenames( const std::string & path )
        {
            //LG(INFO, "listFilenames(%s)", (path.c_str() ? path.c_str() : "nullptr"));
            std::vector<std::string> filenames;
            
#ifdef _WIN32
            WIN32_FIND_DATA ffd;
            TCHAR szDir[MAX_PATH];
            size_t length_of_arg;
            HANDLE hFind = INVALID_HANDLE_VALUE;
            DWORD dwError = 0;
            
            // Check that the input path plus 3 is not longer than MAX_PATH.
            // Three characters are for the "\*" plus nullptr appended below.
            
            TCHAR tstrTo[MAX_PATH*2];
            const int nMax = sizeof(tstrTo) / sizeof(tstrTo[0]);
            int tstrLen;
#ifdef UNICODE
            tstrLen = MultiByteToWideChar(CP_ACP, 0, path.c_str(), strlen(path.c_str()), nullptr, 0);
            if ( unlikely(tstrLen >= nMax) ) {
                LG(ERR, "listFilenames : string %s is tool long", path.c_str());
                Assert(0);
                return filenames;
            }
            tstrTo[tstrLen] = 0;
            MultiByteToWideChar(CP_ACP, 0, path.c_str(), strlen(path.c_str()), tstrTo, tstrLen);
#else
            int err = strcpy_s( tstrTo, nMax, path.c_str() );
            if ( err != 0 )
            {
                LG(ERR, "listFilenames : strcpy_s error %d", err);
                Assert(0);
                return filenames;
            }
            tstrLen = strlen( tstrTo );
#endif
            
            HRESULT hr=StringCchLength(tstrTo, MAX_PATH, &length_of_arg);
            
            if (unlikely(FAILED(hr)))
            {
                LG(ERR, "listFilenames : StringCchLength failed (%x)", hr);
            }
            else if (unlikely(length_of_arg > (MAX_PATH - 3)))
            {
                // can fix this by using unicode version of FindFirstFile and prepending \\?\ to the path
                LG(ERR, "listFilenames : Directory path is too long");
            }
            else
            {
                // Prepare string for use with FindFile functions.  First, copy the
                // string to a buffer, then append '\*' to the directory name.
                
                StringCchCopy(szDir, MAX_PATH, tstrTo);
                StringCchCat(szDir, MAX_PATH, TEXT("\\*"));
                
                // Find the first file in the directory.
                
                hFind = FindFirstFile(szDir, &ffd);
                
                if (unlikely(INVALID_HANDLE_VALUE == hFind))
                {
                    LG(INFO, "listFilenames : FindFirstFile returned INVALID_HANDLE_VALUE");
                }
                else
                {
                    // List all the files in the directory with some info about them.
                    
                    do
                    {
                        if (!(ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
                        {
                            std::string cast;
                            string_cast(ffd.cFileName, CP_ACP, cast);
                            filenames.push_back(cast);
                        }
                    } while (FindNextFile(hFind, &ffd) != 0);
                    
                    dwError = GetLastError();
                    if (unlikely(dwError != ERROR_NO_MORE_FILES))
                    {
                        LG(ERR, "listFilenames : FindNextFile returned %d", dwError);
                    }
                }
                FindClose(hFind);
            }
            
#else
            DIR           *d;
            struct dirent *dir;
            d = opendir(path.c_str());
            if (likely(d))
            {
                while ((dir = readdir(d)) != nullptr)
                {
                    if (dir->d_type == DT_REG)
                    {
                        filenames.push_back(dir->d_name);
                    }
                }
                
                closedir(d);
            }
#endif
            
            //LG(INFO, "listFilenames(%s) found %d files and returns %s", (path.c_str() ? path.c_str() : "nullptr"), filenames.size(), (bExists ? "true" : "false"));
            return filenames;
        }

        
        std::vector< std::string > listDirectories( const DirectoryPath & path ) {
            return listDirectories(path.toString());
        }
        std::vector< std::string > listDirectories( const std::string & path )
        {
            std::vector<std::string> dirnames;
            
#ifdef _WIN32
            throw std::logic_error("not implemented on windows");
#else
            DIR           *d;
            struct dirent *dir;
            d = opendir(path.c_str());
            if (likely(d))
            {
                while ((dir = readdir(d)) != nullptr)
                {
                    if(!strcmp(dir->d_name, ".")) {
                        continue;
                    }
                    if (dir->d_type == DT_DIR)
                    {
                        dirnames.push_back(dir->d_name);
                    }
                }
                
                closedir(d);
            }
#endif
            
            return dirnames;
        }
        

    bool isGUID(std::string const & str)
    {
        bool bIsGUID = true;
        
        int count = -1;
        int idx_parenthesis_open = -1;
        int idx_parenthesis_close = -1;
        for(char const &c:str)
        {
            count++;
            
            switch(c)
            {
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case 'A':
                case 'B':
                case 'C':
                case 'D':
                case 'E':
                case 'F':
                case 'a':
                case 'b':
                case 'c':
                case 'd':
                case 'e':
                case 'f':
                    break;
                case '-':
                    break;
                case '{':
                    if(idx_parenthesis_open != -1) {
                        bIsGUID = false;
                    }
                    else {
                        idx_parenthesis_open = count;
                    }
                    break;
                case '}':
                    if(idx_parenthesis_close != -1) {
                        bIsGUID = false;
                    }
                    else {
                        idx_parenthesis_close = count;
                    }
                    break;
                default:
                    bIsGUID = false;
                    break;
            }
            if (!bIsGUID) {
                break;
            }
        }
        
        if(bIsGUID) {
            if( idx_parenthesis_close != -1 ) {
                if( idx_parenthesis_open != -1 ) {
                    if( ( idx_parenthesis_open != 0 ) || ( idx_parenthesis_close != count ) ) {
                        bIsGUID = false;
                    }
                } else {
                    bIsGUID = false;
                }
            } else if( idx_parenthesis_open != -1 ) {
                bIsGUID = false;
            }
        }
        return bIsGUID;
    }
    
    } // namespace StorageStuff
} // namespace imajuscule

DirectoryPath::DirectoryPath(const std::string & sInput) {
    set(sInput);
}

DirectoryPath::DirectoryPath(const char * sInput) {
    set(std::string(sInput));
}
void DirectoryPath::set(const std::string & sInput)
{
    std::istringstream f(sInput);
    std::string s;
    while (getline(f, s, '/')) {
        vec.push_back(s);
    }
}

std::string DirectoryPath::toString() const
{
    std::string ret;
    for( auto & st : vec )
    {
        ret.append(st);
        ret.append("/");
    }
    
    if(ret.size() > 1) {
        // remove trailing '/'
        ret.pop_back();
    }
    
    return ret;
}
