/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
    // implementation of custom (imajuscule) fft
    
    namespace imj {
        struct Tag {};
    }
    
    namespace fft {
        
        /*
         * Space complexity, for forward fft of real input of size N:
         *
         * input : 2*N
         */

        template<typename T>
        struct RealInput_<imj::Tag, T> {
            using type = std::vector<complex<T>>;
        };
        
        template<typename T>
        struct RealOutput_<imj::Tag, T> {
            using type = std::vector<complex<T>>;
        };
        
        template<typename T>
        struct ImjContext {
            using vec_roots = std::vector<complex<T>>;
            
            ImjContext(vec_roots * roots) : roots(roots) {}
            
            vec_roots * getRoots() const { return roots; }
            vec_roots * editRoots() { return roots; }
        private:
            vec_roots * roots;
        };
        
        template<typename T>
        struct Context_<imj::Tag, T> {
            using type = ImjContext<T>;
            using InnerCtxt = typename type::vec_roots;
            
            static auto create(int size) {
                auto pv = new InnerCtxt();
                compute_roots_of_unity(size, *pv);
                return type(pv);
            }

            static void destroy(type c) {
                delete c.editRoots();
            }
        };
        
        template<typename T>
        struct Algo_<imj::Tag, T> {
            using RealInput  = typename RealInput_ <imj::Tag, T>::type;
            using RealOutput = typename RealOutput_<imj::Tag, T>::type;
            using Context    = typename Context_   <imj::Tag, T>::type;
            
            Algo_(Context c) : ctxt(c) {}
            
            void run(RealInput const & input,
                     RealOutput & output,
                     unsigned int N) const
            {
                constexpr auto stride = 1;
                tukeyCooley(ctxt.getRoots()->begin(),
                            input.begin(),
                            output.begin(),
                            N/2,
                            stride);
            }
            
            Context ctxt;
        };
    } // NS fft
    
    namespace imj {
        namespace fft {
            using namespace imajuscule::fft;
            
            // this part could be #included to avoid repetitions
            
            template<typename T>
            using RealInput = typename RealInput_<Tag, T>::type;

            template<typename T>
            using RealOutput = typename RealOutput_<Tag, T>::type;
        
            template<typename T>
            using Context = typename Context_<Tag, T>::type;
            
            template<typename T>
            using ScopedContext = ScopedContext_<Tag, T>;
            
            template<typename T>
            using Algo = Algo_<Tag, T>;
        } // NS fft
    } // NS imj
} // NS imajuscule

