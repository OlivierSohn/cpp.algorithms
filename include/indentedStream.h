namespace imajuscule
{
    class indentedStream {
    public:
        indentedStream() : beginline(true) {}
        template <typename T>
        indentedStream & operator << (T t) {
            handleIndent();
            return add(t);
        }
        
        std::string str() const {
            return stream.str();
        }
        void addIndent() {
            indent_level++;
        }
        void removeIndent() {
            indent_level--;
        }
    private:
        std::stringstream stream;
        bool beginline : 1;
        int8_t indent_level = 0;
        void handleIndent() {
            if(beginline) {
                beginline = false;
                add(std::string(indent_level * 4, ' '));
            }
        }
        template <typename T>
        indentedStream & add(T t) {
            stream << t;
            return *this;
        }
    };
    
    template <>
    indentedStream & indentedStream::operator << <const char*>(const char * t);
}
