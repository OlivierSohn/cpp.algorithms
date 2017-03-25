
namespace imajuscule {
    namespace fft {
        
        void compute_fft(unsigned int fft_length,
                         FFTIter<double> signal_it,
                         FFTIter<double> result_it) {
            auto roots = compute_roots_of_unity<double>(fft_length);
            Algo<double> a(roots);
            a.run(signal_it, result_it, fft_length, 1);
        }
        
    } // NS fft
} // NS imajuscule

