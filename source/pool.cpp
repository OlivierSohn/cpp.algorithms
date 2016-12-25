
namespace imajuscule {
    
    Pool * Pool::instance = nullptr;
    Pool::State Pool::state = Pool::Growing;
    
    Pool & Pool::getInstance() {
        if(!instance) {
            instance = new Pool();
        }
        return *instance;
    }
    
    void Pool::allocate_overflow() {
        overflow.reset( new Pool(buffer.size()) );
    }
    
} // ns imajuscule
