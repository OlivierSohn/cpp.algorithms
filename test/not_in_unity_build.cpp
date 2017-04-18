//
// this is a translation unit separate from unity build
//

#include "public.h"

namespace imajuscule {
    namespace threadlocaltest {
        using Ctxts = fft::Contexts_<imj::Tag, float>;
  
        Ctxts & get_other_context() {
            return Ctxts::getInstance();
        }
    }
}
