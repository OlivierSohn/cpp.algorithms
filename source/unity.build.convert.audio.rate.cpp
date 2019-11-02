/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

#include "unity.build.cpp"

namespace imajuscule::audio {
    static inline void convert_rate(DirectoryPath const & dir, FileName const & filename) {
        resample_wav(dir,
                     filename,
                     "resampled_" + filename,
                     //44100.,
                     [](double x){
            if(x<200000) {
                return 44100.;
            }
            constexpr auto transition_width = 10000;
            if(x > 200000+transition_width) {
                return 30000.;
            }
            else {
                return 44100. + ((x-200000)/transition_width) * (30000.-44100.);
            }
        },
                     44100);
    }
}


void printUsage() {
    using namespace std;
    cout << "- convert audio rate usage : " << endl;
    cout << "- pass one argument containing the path to the wav file." << endl;
    cout << "- this wav file will be converted to 44.1kHz" << endl;
}

int main(int argc, const char * argv[]) {
    using namespace std;
    using namespace imajuscule;
    using namespace imajuscule::audio;
    
    if(argc != 2) {
        cerr << "1 argument is needed, " << argc-1 << " given" << endl;
        printUsage();
        throw;
    }
    
    DirectoryPath dir;
    FileName filename;
    if(!split_path(argv[1], dir, filename)) {
        cerr << "invalid path : " << argv[1] << endl;
        throw;
    }
   
    convert_rate(dir, filename);

    return 0;
}
