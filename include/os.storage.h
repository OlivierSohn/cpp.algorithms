
namespace imajuscule {

enum eResult
{
  ILE_SUCCESS = 0,
  ILE_RECURSIVITY,
  ILE_ERROR,
  ILE_NOT_IMPLEMENTED,
  ILE_BAD_PARAMETER,
  ILE_OBJECT_INVALID
};

#define SIZE_READ_BUFFER 2048
namespace StorageStuff {
#if defined (_MSC_VER)
void string_cast(const wchar_t* pSource, unsigned int codePage /*= CP_ACP*/, std::string & oString);
#endif

bool dirExists(const std::string & path);
bool fileExists(const std::string & path);
bool fileCreationDate(const std::string & path, std::string & oDate);
eResult makeDir(const std::string & path);
eResult removeDir(const std::string & path);
eResult removeFile(const std::string & path);

bool isGUID(std::string const & str);
}

class WritableStorage;
class DirectoryPath {
  std::filesystem::path m_path;
  
  static DirectoryPath referentiablesPath;
  static DirectoryPath capturePath;
public:
  DirectoryPath( const std::filesystem::path & path );
    
  const std::filesystem::path & path() const {return m_path;}
  
  static DirectoryPath root();
  static void setReferentiablesDir(DirectoryPath const &);
  static bool getReferentiablesDir(DirectoryPath &);
  static void setCaptureImageDir(DirectoryPath const &);
  static bool getCaptureImageDir(DirectoryPath &);
  
  bool empty() const { return m_path.empty(); }
  bool isFile() const { return std::filesystem::is_regular_file(m_path); }
  bool isDir() const { return std::filesystem::is_directory(m_path); }
    
  std::string toString() const;
  void set(const std::string & path);
  
  DirectoryPath operator + ( const DirectoryPath & other ) const {
    return { m_path / other.m_path };
  }
  void operator += ( const DirectoryPath & other ) {
    m_path = m_path / other.m_path;
  }
  DirectoryPath & append(const char * p) {
    m_path = m_path / std::filesystem::path{p};
    return *this;
  }
};

using FileName = std::string;

bool split_path(std::string const & str, DirectoryPath & dir, FileName & filename);
std::string getParentDirectory(std::string path);

class WritableStorage;

struct ReadableStorage {
  friend class WritableStorage;
  
  auto const & path() {
    return m_path;
  }
  
  enum class FileMode
  {
    WRITE,
    READ
  };
  
protected:
  ReadableStorage(std::filesystem::path const & p);
  virtual ~ReadableStorage();
  
  eResult OpenForRead();
  void CloseFile();
  
  // returns the size read
  size_t ReadData(void * p, size_t size, size_t count);
  
  
protected:
  void* m_pFile = nullptr;
private:
  unsigned char m_freadBuffer[SIZE_READ_BUFFER];
  size_t m_bufferReadPos;
  
protected:
  std::filesystem::path m_path;
  
private:
  eResult OpenFileForOperation(const std::string & sFilePath, FileMode);
  void ReadToBuffer();
};

class WritableStorage : public ReadableStorage
{
public:
  
  eResult Save();
  
protected:
  
  WritableStorage(std::filesystem::path const & p) : ReadableStorage(p) {}
  
  ~WritableStorage() {
    std::lock_guard l(g_mutex_openedForWrite);

    g_openedForWrite.erase(m_filePath);
  }
  
  eResult OpenForWrite();
  
  virtual int WriteData(void const * p, size_t size, size_t count);
  
  void Finalize();
  
  void UpdateFileHeader();
  
  // child classes should call this method directly only the first time the header is written.
  // for subsequent header writes they should call instead UpdateFileHeader that will call this method at the appropriate moment
  // and then restore the file position to the position it had before writing the header
  virtual void DoUpdateFileHeader() { Assert(0); }
  
  eResult doSaveBegin();
private:
  bool isBeingSaved();
  virtual eResult doSave();
  void doSaveEnd();
  
  std::vector<unsigned char> m_writeBuffer;
  
  static std::set<std::filesystem::path> g_openedForWrite;
  static std::mutex g_mutex_openedForWrite;
  std::filesystem::path m_filePath;
  
  int  FlushData();
  void FlushMyBuffer();
  
  
};
namespace StorageStuff {
std::vector< std::string > listFilenames( const DirectoryPath & path );
std::vector< std::string > listFilenames( const std::string & path );
std::vector< std::string > listDirectories( const DirectoryPath & path );
std::vector< std::string > listDirectories( const std::string & path );

const char * FileOperationToString(WritableStorage::FileMode op);
}

}
