
namespace imajuscule {
    template <>
    indentedStream & indentedStream::operator << <const char *>(const char * c){
        handleIndent();
        if(strstr(c, "\n")) {
            beginline = true;
        }
        //LG(INFO, "add: %s", c);
        return add(c);
    }
}

