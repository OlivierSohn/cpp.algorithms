namespace imajuscule
{
    // https://stackoverflow.com/questions/9599807/how-to-add-indention-to-the-stream-operator
    class IndentingOStreambuf : public std::streambuf
    {
        std::streambuf*     myDest;
        bool                myIsAtStartOfLine;
        std::string         myIndent;
        std::ostream*       myOwner;
    protected:
        virtual int         overflow( int ch )
        {
            if ( myIsAtStartOfLine && ch != '\n' ) {
                myDest->sputn( myIndent.data(), myIndent.size() );
            }
            myIsAtStartOfLine = ch == '\n';
            return myDest->sputc( ch );
        }
    public:
        explicit            IndentingOStreambuf(
                                std::streambuf* dest, int indent = 4 )
            : myDest( dest )
            , myIsAtStartOfLine( true )
            , myIndent( indent, ' ' )
            , myOwner( NULL )
        {
        }
        explicit            IndentingOStreambuf(
                                std::ostream& dest, int indent = 4 )
            : myDest( dest.rdbuf() )
            , myIsAtStartOfLine( true )
            , myIndent( indent, ' ' )
            , myOwner( &dest )
        {
            myOwner->rdbuf( this );
        }
        virtual             ~IndentingOStreambuf()
        {
            if ( myOwner != NULL ) {
                myOwner->rdbuf( myDest );
            }
        }
    };
}
