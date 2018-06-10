

namespace imajuscule {
    
    bool get_file_contents(const std::string &filename, std::string & str )
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if (!in)
        {
            LG(WARN, "file not found : %s", filename.c_str());
            return false;
        }
        in.seekg(0, std::ios::end);
        int res = in.tellg();
        if(-1 == res) {
            LG(ERR, "tellg error");
            return false;
        }
        Assert(res >= 0);
        str.resize(safe_cast<size_t>(res));
        in.seekg(0, std::ios::beg);
        in.read(&str[0], str.size());
        in.close();
        return true;
    }
}
