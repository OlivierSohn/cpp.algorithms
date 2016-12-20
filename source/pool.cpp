
namespace imajuscule {
    
    Pool * Pool::instance = nullptr;

    Pool & Pool::getInstance() {
        if(!instance) {
            instance = new Pool();
        }
        return *instance;
    }
    
    void Pool::allocate_overflow() {
        overflow.reset( new Pool(elems.size()) );
    }
    
} // ns imajuscule
