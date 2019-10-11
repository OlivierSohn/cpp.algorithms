

namespace imajuscule::audio {
    void makeDescription(std::string &desc, int16_t num_channels, float length_in_seconds, int32_t sample_rate) {
        desc.clear();
        desc.reserve(14 + 10);
        
        makeDescription(desc, num_channels, length_in_seconds);
        
        // to avoid float rounding error when doing '/ 1000.f' :
        auto s1 = sample_rate/1000;
        auto s2 = sample_rate - 1000 * s1;
        desc += " " + NumTraits<float>::toSignificantString(s1);
        if(s2) {
            auto str2 = NumTraits<float>::toSignificantString(s2);
            while(!str2.empty()) {
                if(str2.back() != '0') {
                    break;
                }
                str2.pop_back();
            }
            desc += "." + str2;
        }
        desc += " kHz";
    }
    
    void makeDescription(std::string &desc, int16_t num_channels, float length_in_seconds) {
        desc.clear();
        desc.reserve(14);
        
        desc.append(num_channels, '.');
        desc += " ";
        constexpr auto n_decimals_total = 2;
        auto str_sec = NumTraits<float>::to_string_with_precision(length_in_seconds, n_decimals_total);
        {
            auto i = str_sec.find('.');
            if(i == std::string::npos) {
                str_sec.push_back('.');
            }
        }
        {
            auto i = str_sec.find('.');
            int n_decimals = str_sec.size() - i - 1;
            if(n_decimals < n_decimals_total) {
                str_sec.append(n_decimals_total - n_decimals, '0');
            }
        }
        desc += str_sec + " s";
    }
    
    void WAVPCMHeader::makeDescription(std::string & desc, bool with_sample_rate) const {
        using imajuscule::audio::makeDescription;
        if(with_sample_rate) {
            makeDescription(desc, num_channels, getLengthInSeconds(), sample_rate);
        }
        else {
            makeDescription(desc, num_channels, getLengthInSeconds());
        }
    }
    
    bool getConvolutionReverbSignature(std::string const & dirname, std::string const & filename, spaceResponse_t & r) {
        WAVReader reader(dirname, filename);
        
        try {
            reader.Initialize();
        }
        catch(std::exception const & e) {
            LG(WARN, "Cannot read '%s' as a '.wav' file. If the file exists in '%s', it might be corrupted.", filename.c_str(), dirname.c_str());
            return false;
        }
        
        r.nChannels = reader.countChannels();
        r.sampleRate = reader.getSampleRate();
        r.nFrames = reader.countFrames();
        r.sampleSize = reader.getSampleSize();
        r.lengthInSeconds = reader.getLengthInSeconds();
        return true;
    }
    
    template<typename WAVRead>
    void readReverb(int nouts, double sample_rate, WAVRead & reader, InterlacedBuffer & ib)
    {
        reader.Initialize();
        
        auto mod = reader.countChannels() % nouts;
        if((reader.countChannels() > nouts) && mod) {
            std::ostringstream msg;
            msg << "Cannot use a '" << reader.countChannels() << "' channels reverb for '" << nouts << "' output channels. The reverb channels count must be a multiple of the output channels count.";
            throw std::runtime_error(msg.str());
        }
        
        double stride = reader.getSampleRate() / sample_rate;
        
        ib.nchannels = reader.countChannels();
        ib.buf.resize(static_cast<int>(reader.countFrames() / stride) * reader.countChannels());
        
        MultiChannelReSampling<decltype(reader)> mci(reader);
        mci.Read(ib.buf.begin(), ib.buf.end(), stride);
    }
    
    void readReverbFromBuffer(int nouts, double sample_rate, void const * buffer, std::size_t const sz, InterlacedBuffer & ib) {
        WAVReaderFromBlock reader(buffer, sz);
        readReverb(nouts, sample_rate, reader, ib);
    }
    
    void readReverbFromFile(int nouts, double sample_rate, std::string const & dirname, std::string const & filename, InterlacedBuffer & ib) {
        WAVReader reader(dirname, filename);
        readReverb(nouts, sample_rate, reader, ib);
    }
}
