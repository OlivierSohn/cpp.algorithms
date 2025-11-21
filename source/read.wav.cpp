

namespace imajuscule::audio {
    std::string to_string(WaveFormat f) {
        switch(f) {
            case WaveFormat::PCM:
                return "PCM";
            case WaveFormat::IEEE_FLOAT:
                return "IEEE_FLOAT";
            case WaveFormat::ALAW:
                return "ALAW";
            case WaveFormat::MULAW:
                return "MULAW";
            case WaveFormat::EXTENSIBLE:
                return "EXTENSIBLE";
            default:
                return "???";
        }
    }

    std::string formatSampleRateKHz(int32_t sample_rate) {
        std::string desc;
        // to avoid float rounding error when doing '/ 1000.f' :
        auto s1 = sample_rate/1000;
        auto s2 = sample_rate - 1000 * s1;
        desc = NumTraits<float>::toSignificantString(s1);
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
        return desc;
    }
    
    void makeDescription(std::string &desc, int16_t num_channels, float length_in_seconds, int32_t sample_rate) {
        desc.clear();
        desc.reserve(14 + 10);
        
        makeDescription(desc, num_channels, length_in_seconds);
        
        desc += " " + formatSampleRateKHz(sample_rate) + " kHz";
    }
        
    void makeDescription(std::string &desc, int16_t num_channels, float length_in_seconds) {
        desc.clear();
        desc.reserve(14);
        
        desc.append(num_channels, '.');
        desc += " ";
        constexpr auto n_decimals_total = 2;
        auto str_sec = formatSecondsDurationWithPrecision(length_in_seconds, n_decimals_total);
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
    
    bool getConvolutionReverbSignature(std::filesystem::path const & filepath, spaceResponse_t & r) {
        WAVReader reader(filepath);
        
        try {
            reader.Initialize();
        }
        catch(std::exception const & e) {
            LG(WARN, "Cannot read '%s' as a '.wav' file. If the file exists, it might be corrupted.", filepath.u8string().c_str());
            return false;
        }
        
        r.nChannels = reader.countChannels();
        r.sampleRate = reader.getSampleRate();
        r.nFrames = reader.countFrames();
        r.sampleSize = reader.getSampleSize();
        r.lengthInSeconds = reader.getLengthInSeconds();
        return true;
    }
    
}
