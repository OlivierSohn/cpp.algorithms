

namespace imajuscule {

template<typename Left, typename Right>
struct Either {
    Either(Left const & l) : var(l) {}
    Either(Right const & r) : var(r) {}
    
    template<typename FLeft, typename FRight>
    auto either(FLeft leftF, FRight rightF) const {
        if(var.index() == 0) {
            return leftF(std::get<Left>(var));
        }
        else {
            return rightF(std::get<Right>(var));
        }
    }

private:
    std::variant<Left, Right> var;
};

template<typename Left,  typename Right>
std::ostream & operator << (std::ostream &ss, const Either<Left, Right> & r) {
    r.either([&ss](auto & err){
        ss << "error:" << std::endl;
        {
            imajuscule::IndentingOStreambuf i(ss);
            ss << err;
        }
    }, [&ss](auto & suc) {
        ss << "success:" << std::endl;
        {
            imajuscule::IndentingOStreambuf i(ss);
            ss << suc;
        }
    });
    ss << std::endl;
    return ss;
}

/*
 Using these leads to ugly code because we need to explicitely write the Left and Right types, for example:
    return left<std::string, result>("error: bla");
 
 This is why I made the costructors non explicit, so that we can write:
    return {"error: bla"};
 */
template<typename Left, typename Right>
Either<Left, Right> left(Left const & l) {
    return {l};
}
template<typename Left, typename Right>
Either<Left, Right> right(Right const & r) {
    return {r};
}

} // NS imajuscule

