
namespace imajuscule::audio {

        // http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html

        enum class WaveFormat : uint16_t {
            PCM =          0x0001, // PCM (Linear quantization, represented using signed values)
            IEEE_FLOAT =   0x0003, // IEEE float
            ALAW =         0x0006, // 8-bit ITU-T G.711 A-law
            MULAW =        0x0007, // 8-bit ITU-T G.711 µ-law
            EXTENSIBLE =   0xFFFE  // Determined by SubFormat
        };

    std::string to_string(WaveFormat f);

        void makeDescription(std::string &s, int16_t num_channels, float length_in_seconds, int32_t sample_rate);
        void makeDescription(std::string &s, int16_t num_channels, float length_in_seconds);

    std::string formatSampleRateKHz(int32_t sample_rate);

        struct WAVPCMHeader {
            std::array<int8_t, 4> chunk_id;
            int32_t chunk_size;
            std::array<int8_t, 4> format;
            std::array<int8_t, 4> subchunk1_id;
            int32_t subchunk1_size;
            WaveFormat audio_format;
            int16_t num_channels;
            int32_t sample_rate;
            int32_t byte_rate;
            int16_t block_align;
            int16_t bits_per_sample;
            std::array<int8_t, 4> subchunk2_id;
            int32_t subchunk2_size; // subchunk2_size denotes the number of samples.

            unsigned int countFrames() const {
                if(0 == block_align) {
                    // header has not been read yet
                    return 0;
                }
                Assert(0 == subchunk2_size % block_align);
                Assert(block_align);
                return subchunk2_size / block_align;
            }

            unsigned int countSamples() const {
                return countFrames() * num_channels;
            }

            unsigned int getSampleSize() const {
                Assert(bits_per_sample % 8 == 0);
                return bits_per_sample / 8;
            }

            float getLengthInSeconds() const {
                return countFrames() / static_cast<float>(sample_rate);
            }

            void makeDescription(std::string & desc, bool with_sample_rate) const;
        };

        enum class NChannels {
            ONE = 1,
            TWO = 2,
            MONO = 1,
            STEREO = 2
        };

        template<typename T>
        NChannels numberToNChannels(T n) {
            if (n==1) {
                return NChannels::ONE;
            }
            if (n==2) {
                return NChannels::TWO;
            }
            Assert(0);
            return NChannels::ONE;
        }

        enum class WaveFileFormat {
            BeginSupported,

            RIFF = BeginSupported,
            WaveIR,

            EndSupported,

            Unknown = EndSupported
        };

        constexpr auto getWavChunkId(WaveFileFormat f) {
            switch(f) {
                case WaveFileFormat::RIFF:
                    return "RIFF";
                case WaveFileFormat::WaveIR:
                    return "wvIR";
                default:
                    return "XXXX";
            }
        }
        constexpr auto getWavWAVEId(WaveFileFormat f) {
            switch(f) {
                case WaveFileFormat::RIFF:
                    return "WAVE";
                case WaveFileFormat::WaveIR:
                    return "ver1";
                default:
                    return "XXXX";
            }
        }

        static bool isFormatSupported(WaveFormat format) {
            switch(format) {
                case WaveFormat::PCM:
                case WaveFormat::IEEE_FLOAT:
                    return true;
                default:
                    return false;
            }
        }

        static WAVPCMHeader pcm_(WaveFileFormat const file_format,
                                 WaveFormat const format,
                                 int const num_frames,
                                 int const sample_rate,
                                 NChannels const n_channels,
                                 int const bytes_per_sample) {
            Assert(bytes_per_sample);
            auto const num_channels = to_underlying(n_channels);
            auto const size_data = num_frames * num_channels * bytes_per_sample;
            WAVPCMHeader ret {
                {0,0,0,0},
                static_cast<int32_t>(sizeof(WAVPCMHeader) + size_data),
                {0,0,0,0},
                {'f','m','t',' '},
                16,
                format,
                static_cast<int16_t>(num_channels),
                sample_rate,
                sample_rate * bytes_per_sample * num_channels,
                static_cast<int16_t>(bytes_per_sample * num_channels),
                static_cast<int16_t>(bytes_per_sample * 8),
                {'d','a','t','a'},
                size_data
            };

            memcpy(ret.chunk_id.data(), getWavChunkId(file_format), 4);
            memcpy(ret.format.data(), getWavWAVEId(file_format), 4);

            return ret;
        }

        enum class SampleFormat {
            signed_16,
            signed_24,
            signed_32,
            float_32,
            float_64
        };

        template<typename T>
        struct AudioSample;

        template<>
        struct AudioSample<int16_t> {
            static constexpr auto format = SampleFormat::signed_16;
        };

        template<>
        struct AudioSample<int24_t> {
            static constexpr auto format = SampleFormat::signed_24;
        };

        template<>
        struct AudioSample<int32_t> {
            static constexpr auto format = SampleFormat::signed_32;
        };

        template<>
        struct AudioSample<float> {
            static constexpr auto format = SampleFormat::float_32;
        };

        template<>
        struct AudioSample<double> {
            static constexpr auto format = SampleFormat::float_64;
        };

        template<int SAMPLE_SIZE>
        struct SignedTypeForSampleSize;

        template<int SAMPLE_SIZE>
        struct FloatTypeForSampleSize;

        template<WaveFormat FMT, int SAMPLE_SIZE>
        struct TypeFor;


        template <>
        struct SignedTypeForSampleSize<1> {
            using type = int8_t;
            static_assert(sizeof(type) == 1);
        };

        template <>
        struct SignedTypeForSampleSize<2> {
            using type = int16_t;
            static_assert(sizeof(type) == 2);
        };

        template <>
        struct SignedTypeForSampleSize<3> {
            using type = int24_t;
            static_assert(sizeof(type) == 3);
        };

        template <>
        struct SignedTypeForSampleSize<4> {
            using type = int32_t;
            static_assert(sizeof(type) == 4);
        };


        template <>
        struct FloatTypeForSampleSize<4> {
            using type = float;
            static_assert(sizeof(type) == 4);
        };

        template <>
        struct FloatTypeForSampleSize<8> {
            using type = double;
            static_assert(sizeof(type) == 8);
        };

        template<int SAMPLE_SIZE>
        struct TypeFor<WaveFormat::PCM, SAMPLE_SIZE> {
            using type = typename SignedTypeForSampleSize<SAMPLE_SIZE>::type;
        };

        template<int SAMPLE_SIZE>
        struct TypeFor<WaveFormat::IEEE_FLOAT, SAMPLE_SIZE> {
            using type = typename FloatTypeForSampleSize<SAMPLE_SIZE>::type;
        };

        constexpr auto bytes_per_sample(SampleFormat f) {
            switch(f) {
                case SampleFormat::signed_16:
                    return 2;
                case SampleFormat::signed_24:
                    return 3;
                case SampleFormat::signed_32:
                case SampleFormat::float_32:
                    return 4;
                case SampleFormat::float_64:
                    return 8;
            }
            return 0;
        }

        static
#ifdef NDEBUG
        constexpr
#endif
        bool is_signed(SampleFormat const f) {
            switch(f) {
                case SampleFormat::signed_16:
                case SampleFormat::signed_24:
                case SampleFormat::signed_32:
                    return true;
                case SampleFormat::float_32:
                case SampleFormat::float_64:
                    return false;
            }

            Assert(0);
            return false;
        }

        static
#ifdef NDEBUG
        constexpr
#endif
        bool is_float(SampleFormat const f) {
            switch(f) {
                case SampleFormat::signed_16:
                case SampleFormat::signed_24:
                case SampleFormat::signed_32:
                    return false;
                case SampleFormat::float_32:
                case SampleFormat::float_64:
                    return true;
            }

            Assert(0);
            return false;
        }

        static bool are_compatible(WaveFormat const format, SampleFormat const sf) {
            switch(format) {
                case WaveFormat::PCM:
                    return is_signed(sf);

                case WaveFormat::IEEE_FLOAT:
                    return is_float(sf);

                default:
                    LG(ERR, "unhandled wave format");
                    return false;
            }
        }

        static inline WAVPCMHeader pcm(WaveFormat format, int sample_rate, NChannels n_channels, SampleFormat f) {
            Assert(are_compatible(format, f));
            return pcm_(WaveFileFormat::RIFF,
                        format,
                        0, // number of bytes for samples will be deduced later
                        sample_rate,
                        n_channels,
                        bytes_per_sample(f));
        }

        template<typename Parent>
        struct WAVReader_ : public Parent {
            using Parent::OpenForRead;
            using Parent::ReadData;

            template <typename ...Params>
            WAVReader_(Params&&... params)
            : Parent(std::forward<Params>(params)...)
            , header{}
            {}

            void Initialize()
            {
                auto res = OpenForRead();
                if(res != ILE_SUCCESS) {
                    throw std::runtime_error("open for read failed: " + std::to_string(res));
                }
                readHeader();
            }

            void makeDescription(std::string &s, bool with_sample_rate = false) const {
                header.makeDescription(s, with_sample_rate);
            }

            auto const & getHeader() const { return header; }

            unsigned int getSampleSize() const { return header.getSampleSize(); }
            int countChannels() const { return header.num_channels; }
            int countSamples() const { return header.countSamples(); }
            int countFrames() const { return header.countFrames(); }
            int getSampleRate() const { return header.sample_rate; }
            float getLengthInSeconds() const { return header.getLengthInSeconds(); }
            WaveFormat getFormat() const { return header.audio_format; }
            WaveFileFormat getFileFormat() const {
                auto candidate = WaveFileFormat::BeginSupported;
                for(;
                    candidate != WaveFileFormat::EndSupported;
                    increment(candidate))
                {
                    if(!memcmp(header.chunk_id, getWavChunkId(candidate), 4)) {
                        break;
                    }
                }
                return candidate;
            }

            // 'it' is the begin of the buffer to write to, 'end' is the end.
            // returns the end element

            template<typename ITER>
            ITER Read(ITER it, ITER end) {
                constexpr auto n_reads = sizeof(decltype(*it));
                Assert(n_reads == header.getSampleSize());
                auto dist = std::distance(it, end);
                auto bytes_to_read = dist * n_reads;
                Assert((audio_bytes_read + bytes_to_read) <= header.subchunk2_size);
                audio_bytes_read += bytes_to_read;
                while(it < end) {
                    ReadData(&*it, n_reads, 1);
                    ++it;
                }
                return it;
            }

            template<typename T>
            void ReadOne(T& val) {
                constexpr auto n_reads = sizeof(T);
                Assert(n_reads == header.getSampleSize());
                auto bytes_to_read = n_reads;
                Assert((audio_bytes_read + bytes_to_read) <= header.subchunk2_size);
                audio_bytes_read += bytes_to_read;
                ReadData(&val, n_reads, 1);
            }

            bool HasMore() const {
                Assert(audio_bytes_read <= header.subchunk2_size);
                //LG(INFO, "%d %d", audio_bytes_read, header.subchunk2_size);
                return audio_bytes_read < header.subchunk2_size;
            }

        private:
            WAVPCMHeader header;
            unsigned int audio_bytes_read = 0;

            void readHeader() {
                bool is_wave_ir = false;
                ReadData(&header.chunk_id[0], 4, 1);
                {
                    if(memcmp(header.chunk_id.data(), getWavChunkId(WaveFileFormat::RIFF), 4)) {
                        if(!memcmp(header.chunk_id.data(), getWavChunkId(WaveFileFormat::WaveIR), 4)) {
                            is_wave_ir = true;
                        }
                        else {
                            std::string str(header.chunk_id.data(), header.chunk_id.data()+4);
                            throw std::runtime_error("unhandled wav chunk_id : " + str);
                        }
                    }
                }
                ReadData(&header.chunk_size, sizeof(int32_t), 1);
                ReadData(&header.format[0], 4, 1);
                {
                    if(is_wave_ir) {
                        if(memcmp(header.format.data(), "ver1", 4)) {
                            std::string str(header.format.data(), header.format.data()+4);
                            throw std::runtime_error("unhandled wir format : " + str);
                        }
                    }
                    else if(memcmp(header.format.data(), "WAVE", 4)) {
                        std::string str(header.format.data(), header.format.data()+4);
                        throw std::runtime_error("unhandled wav format : " + str);
                    }
                }
                ReadData(&header.subchunk1_id[0], 4, 1);
                {
                    if(memcmp(header.subchunk1_id.data(), "fmt ", 4)) {
                        std::string str(header.subchunk1_id.data(), header.subchunk1_id.data()+4);
                        throw std::runtime_error("unhandled wav subchunk1_id : " + str);
                    }
                }
                ReadData(&header.subchunk1_size, sizeof(int32_t), 1);
                ReadData(&header.audio_format, sizeof(uint16_t), 1);
                if(!isFormatSupported(header.audio_format)) {
                    throw std::runtime_error("unhandled audio format : " + std::to_string(static_cast<uint16_t>(header.audio_format)));
                }
                ReadData(&header.num_channels, sizeof(int16_t), 1);
                ReadData(&header.sample_rate, sizeof(int32_t), 1);
                ReadData(&header.byte_rate, sizeof(int32_t), 1);
                ReadData(&header.block_align, sizeof(int16_t), 1);
                ReadData(&header.bits_per_sample, sizeof(int16_t), 1);


                if(is_wave_ir) {
                    //
                    // wir format found here : http://freeverb3vst.osdn.jp/tips.shtml says
                    // encoding is IEEE_float 32
                    //
                    if(header.bits_per_sample == 23) {
                        // handle a bug (might be an intentional bug to make it harder to decode?)
                        header.bits_per_sample = 32;
                    }
                    else if(header.bits_per_sample != 32) {
                        throw std::runtime_error("invalid bits_per_sample in wir header : " + std::to_string(header.bits_per_sample));
                    }
                }

                // skip irrelevant subchunks, until "data" subchunk

                while(1) {
                    ReadData(&header.subchunk2_id[0], 4, 1);
                    ReadData(&header.subchunk2_size, sizeof(int32_t), 1);
                    {
                        if(!memcmp(header.subchunk2_id.data(), "data", 4)) {
                            break;
                        }
                        if(!memcmp(header.subchunk2_id.data(), "fact", 4)) {
                            // we skip that, cf. http://www-mmsp.ece.mcgill.ca/Documents/AudioFormats/WAVE/WAVE.html
                            // contains the number of samples per channel
                        }
                        else if(!memcmp(header.subchunk2_id.data(), "PEAK", 4)) {
                        }
                        else {
                            std::string str(header.subchunk2_id.data(),
                                            header.subchunk2_id.data()+4);
                            throw std::runtime_error("unrecognized wav subchunk2_id : " + str);
                        }
                        std::vector<uint8_t> skipped(header.subchunk2_size);
                        ReadData(skipped.data(), skipped.size(), 1);
                    }
                }

                auto s = header.getSampleSize();
                if(s < 2 || s > 4) {
                    throw std::runtime_error("unhandled sample size : " + std::to_string(header.getSampleSize()));
                }
            }

            template<int n_bytes, typename FLT>
            FLT ReadOneSignedAsOneFloat() {
                uint8_t d[n_bytes];
                ReadData(d, n_bytes, 1);

                audio_bytes_read += n_bytes;
                Assert(audio_bytes_read <= header.subchunk2_size);

                int32_t v = uint8array_to_int32<n_bytes>(d);
                constexpr int64_t n_different_values = pow2<8 * n_bytes>();
                constexpr int32_t m = -n_different_values/2;
                constexpr int32_t M = n_different_values/2 - 1;
                Assert(v <= M);
                Assert(v >= m);
                return signed_to_float< FLT, int32_t, M, m >(v);
            }

        public:
            template<typename FLT>
            FLT ReadAsOneFloat() {
                Assert(audio_bytes_read < header.subchunk2_size);
                Assert(HasMore());

                auto const sz = header.getSampleSize();

                auto const fmt = getFormat();
                if(fmt == WaveFormat::IEEE_FLOAT) {
                    switch(sz) {
                        case 4:
                        {
                            float f;
                            ReadOne(f);
                            return f;
                        }
                        case 8:
                        {
                            double d;
                            ReadOne(d);
                            return d;
                        }
                    }
                    throw std::logic_error("unsupported sample size for IEEE_FLOAT");
                }
                else if(likely(fmt == WaveFormat::PCM)) {
                    switch(sz) {
                        case 1:
                            return ReadOneSignedAsOneFloat<1, FLT>();
                        case 2:
                            return ReadOneSignedAsOneFloat<2, FLT>();
                        case 3:
                            return ReadOneSignedAsOneFloat<3, FLT>();
                        case 4:
                            return ReadOneSignedAsOneFloat<4, FLT>();
                    }
                    throw std::logic_error("unsupported sample size for PCM");
                }
                throw std::logic_error("unsupported format");
            }
        };

    using WAVReader = WAVReader_<ReadableStorage>;

    struct ConstMemoryBlock {
        ConstMemoryBlock(void const * buffer, std::size_t sz) : buffer(buffer), totalSize(sz) {};

        eResult OpenForRead() const {
            return ILE_SUCCESS;
        }

        void ReadData(void * p, size_t size, size_t count) {
            auto sz = size*count;
            assert(sz+cur <= totalSize);
            memcpy(p,static_cast<const char*>(buffer)+cur, sz);
            cur += sz;
        }

        void const * buffer;
        std::size_t const totalSize;
        std::size_t cur = 0;
    };

    using WAVReaderFromBlock = WAVReader_<ConstMemoryBlock>;

        bool getConvolutionReverbSignature(std::string const & dirname, std::string const & filename, spaceResponse_t & r);

        struct WAVWriter : public WritableStorage {
            WAVWriter(DirectoryPath d, FileName f, WAVPCMHeader h) : WritableStorage(d,f), header(h) {}

            ~WAVWriter() {
                UpdateFileHeader();
                Finalize();
            }

            eResult Initialize() {
                return doSaveBegin();
            }

            WaveFormat getFormat() const { return header.audio_format; }

            template<typename T>
            void writeSample(T s) {
                // T can be integral, int24_t, float, double
                auto constexpr n = sizeof(T);
                WriteData(&s, n, 1);
                n_sample_bytes_written += n;
            }

            void WriteFromOneFloat(float f) {

                auto const sz = header.getSampleSize();

                auto const fmt = getFormat();
                if(fmt == WaveFormat::IEEE_FLOAT) {
                    switch(sz) {
                        case 4:
                        {
                            writeSample(f);
                            return;
                        }
                        case 8:
                        {
                            writeSample(static_cast<double>(f));
                            return;
                        }
                    }
                    throw std::logic_error("unsupported sample size for IEEE_FLOAT");
                }
                else if(likely(fmt == WaveFormat::PCM)) {
                    switch(sz) {
                        case 1:
                            writeSample(float_to_signed<int8_t>(f));
                            return;
                        case 2:
                            writeSample(float_to_signed<int16_t>(f));
                            return;
                        case 3:
                            writeSample(float_to_signed<int24_t>(f));
                            return;
                        case 4:
                            writeSample(float_to_signed<int32_t>(f));
                            return;
                    }
                    throw std::logic_error("unsupported sample size for PCM");
                }
                throw std::logic_error("unsupported format");
            }
        protected:

            void DoUpdateFileHeader() override;

        private:
            int32_t n_sample_bytes_written = 0;
            WAVPCMHeader header;
        };

        template<typename CONTAINER>
        void write_wav(DirectoryPath dir, std::string const & fname, CONTAINER const & samples, NChannels n_channels, int sample_rate) {
            using namespace imajuscule::audio;
            using SampleType = typename CONTAINER::value_type;

            WAVWriter writer(dir,
                             fname,
                             pcm(WaveFormat::IEEE_FLOAT,
                                 sample_rate,
                                 n_channels,
                                 AudioSample<SampleType>::format));
            auto res = writer.Initialize();
            if(res != ILE_SUCCESS) {
                throw;
            }
            for(auto sample : samples) {
                writer.writeSample(sample);
            }
        }

        template<typename VEC_VEC>
        void write_wav(DirectoryPath dir, std::string const & fname, VEC_VEC const & v_samples, int sample_rate) {

            auto sz = v_samples[0].size();
            for(auto const & v : v_samples) {
                Assert(sz == v.size());
            }

            auto n_channels = v_samples.size();

            std::vector<float> interlaced;
            interlaced.reserve(sz * n_channels);
            for(int k=0; k<sz; ++k) {
                for(int i=0; i<n_channels; ++i) {
                    interlaced.push_back(v_samples[i][k]);
                }
            }

            write_wav(std::move(dir), fname, interlaced, static_cast<NChannels>(n_channels), sample_rate);
        }

    template<int SAMPLE_SIZE, WaveFormat FMT, typename S>
    void resample_wav(DirectoryPath const & dir, FileName const & dest, WAVReader & reader, S sample_rate, int32_t wav_target_sample_rate) {
        using namespace std;

        auto const n_channels = reader.countChannels();

        using T = typename TypeFor<FMT, SAMPLE_SIZE>::type;
        static_assert(sizeof(T) == SAMPLE_SIZE);


        Assert(FMT == reader.getFormat());
        ResampleSincStats stats;
        InterlacedBuffer ib(reader, sample_rate, stats);
        std::cout << stats << std::endl;
        Assert(!reader.HasMore());

        auto header = pcm(FMT,
                          wav_target_sample_rate,
                          NChannels(n_channels),
                          AudioSample<T>::format);

        WAVWriter writer(dir, dest, header);
        {
            auto res = writer.Initialize();

            if(ILE_SUCCESS != res) {
                cerr << "could not create file" << endl;
                throw;
            }
        }

        for(auto const & f : ib.getBuffer()) {
            writer.WriteFromOneFloat(f);
        }
    }

    template<typename S>
    void resample_wav(DirectoryPath const & dir, FileName const & source, FileName const & dest, S target_rate, int32_t wav_target_sample_rate) {
        WAVReader reader(dir, source);
        reader.Initialize();
        auto fmt = reader.getFormat();

        switch(reader.getSampleSize()) {
            case 2:
                if(fmt == WaveFormat::PCM) {
                    resample_wav<2, WaveFormat::PCM>(dir, dest, reader, target_rate, wav_target_sample_rate);
                }
                else if(fmt == WaveFormat::IEEE_FLOAT) {
                    throw std::logic_error("unsupported 16-bit IEEE_FLOAT format");
                }
                else {
                    throw std::logic_error("unsupported 16-bit format");
                }
                break;
            case 3:
                if(fmt == WaveFormat::PCM) {
                    resample_wav<3, WaveFormat::PCM>(dir, dest, reader, target_rate, wav_target_sample_rate);
                }
                else if(fmt == WaveFormat::IEEE_FLOAT) {
                    throw std::logic_error("unsupported 24-bit IEEE_FLOAT format");
                }
                else {
                    throw std::logic_error("unsupported 24-bit format");
                }
                break;
            case 4:
                if(fmt == WaveFormat::PCM) {
                    resample_wav<4, WaveFormat::PCM>(dir, dest, reader, target_rate, wav_target_sample_rate);
                }
                else if(fmt == WaveFormat::IEEE_FLOAT) {
                    resample_wav<4, WaveFormat::IEEE_FLOAT>(dir, dest, reader, target_rate, wav_target_sample_rate);
                }
                else {
                    throw std::logic_error("unsupported 32-bit format");
                }
                break;
            case 8:
                if(fmt == WaveFormat::PCM) {
                    throw std::logic_error("unsupported 64-bit PCM format");
                }
                else if(fmt == WaveFormat::IEEE_FLOAT) {
                    resample_wav<8, WaveFormat::IEEE_FLOAT>(dir, dest, reader, target_rate, wav_target_sample_rate);
                }
                else {
                    throw std::logic_error("unsupported 64-bit format");
                }
                break;
            default:
                throw std::logic_error("unsupported xx-bit format");;
        }
    }

        template<int SAMPLE_SIZE, WaveFormat FMT, typename FILTER>
        void filter_frames(DirectoryPath const & dir, FileName const & filename, FileName const & dest,
                           WAVReader & reader,
                           FILTER filter) {
            using namespace std;

            auto const n_channels = reader.countChannels();
            auto const n_frames = reader.countFrames();

            using T = typename TypeFor<FMT, SAMPLE_SIZE>::type;
            static_assert(sizeof(T) == SAMPLE_SIZE);

            Assert(FMT == reader.getFormat());

            auto header = pcm(FMT,
                              reader.getSampleRate(),
                              NChannels(n_channels),
                              AudioSample<T>::format);

            WAVWriter writer(dir, dest, header);
            {
                auto res = writer.Initialize();

                if(ILE_SUCCESS != res) {
                    cerr << "could not create file" << endl;
                    throw;
                }
            }

            std::vector<T> frame(n_channels);

            for(int i=0; i<n_frames; ++i) {

                for(int j=0; j<n_channels; j++) {
                    reader.ReadOne(frame[j]);
                }

                if(!filter(frame)) {
                    continue;
                }

                for(auto c : frame) {
                    writer.writeSample(c);
                }
            }

            Assert(!reader.HasMore());
        }

        template<typename FILTER>
        void filter_frames(DirectoryPath const & dir, FileName const & source, FileName const & dest,
                           FILTER f) {
            using namespace std;

            WAVReader reader(dir, source);
            reader.Initialize();

            auto fmt = reader.getFormat();

            switch(reader.getSampleSize()) {
                case 2:
                    if(fmt == WaveFormat::PCM) {
                        filter_frames<2, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        throw std::logic_error("unsupported 16-bit IEEE_FLOAT format");
                    }
                    else {
                        throw std::logic_error("unsupported 16-bit format");
                    }
                    break;
                case 3:
                    if(fmt == WaveFormat::PCM) {
                        filter_frames<3, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        throw std::logic_error("unsupported 24-bit IEEE_FLOAT format");
                    }
                    else {
                        throw std::logic_error("unsupported 24-bit format");
                    }
                    break;
                case 4:
                    if(fmt == WaveFormat::PCM) {
                        filter_frames<4, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        filter_frames<4, WaveFormat::IEEE_FLOAT>(dir, source, dest, reader, f);
                    }
                    else {
                        throw std::logic_error("unsupported 32-bit format");
                    }
                    break;
                case 8:
                    if(fmt == WaveFormat::PCM) {
                        throw std::logic_error("unsupported 64-bit PCM format");
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        filter_frames<8, WaveFormat::IEEE_FLOAT>(dir, source, dest, reader, f);
                    }
                    else {
                        throw std::logic_error("unsupported 64-bit format");
                    }
                    break;
                default:
                    throw std::logic_error("unsupported xx-bit format");;
            }
        }

        template<int SAMPLE_SIZE, WaveFormat FMT, typename REWRITE>
        void rewrite_wav(DirectoryPath const & dir, FileName const & filename, FileName const & dest,
                           WAVReader & reader,
                           REWRITE rewrite) {
            using namespace std;

            auto const n_channels = reader.countChannels();
            auto const n_frames = reader.countFrames();

            using T = typename TypeFor<FMT, SAMPLE_SIZE>::type;
            static_assert(sizeof(T) == SAMPLE_SIZE);

            Assert(FMT == reader.getFormat());

            std::vector<std::vector<float>> deinterlaced(n_channels);
            for(auto & v : deinterlaced) {
                v.reserve(n_frames);
            }

            while(reader.HasMore()) {
                for(int i=0; i<n_channels; ++i) {
                    deinterlaced[i].push_back(reader.template ReadAsOneFloat<float>());
                }
            }
            Assert(!reader.HasMore());

            rewrite(deinterlaced);

            auto size = -1;
            for(auto const & v : deinterlaced) {
                if(size < 0) {
                    size = v.size();
                    continue;
                }
                if(size != v.size()) {
                    throw std::logic_error("sizes are not the same");
                }
            }

            auto header = pcm(FMT,
                              reader.getSampleRate(),
                              NChannels(n_channels),
                              AudioSample<T>::format);

            WAVWriter writer(dir, dest, header);
            {
                auto res = writer.Initialize();

                if(ILE_SUCCESS != res) {
                    cerr << "could not create file" << endl;
                    throw;
                }
            }

            for(int i=0; i<size; ++i) {
                for(auto & v : deinterlaced) {
                    writer.WriteFromOneFloat(v[i]);
                }
            }
        }

        template<typename REWRITE>
        void rewrite_wav(DirectoryPath const & dir, FileName const & source, FileName const & dest,
                           REWRITE f) {
            using namespace std;

            WAVReader reader(dir, source);
            reader.Initialize();

            auto fmt = reader.getFormat();

            switch(reader.getSampleSize()) {
                case 2:
                    if(fmt == WaveFormat::PCM) {
                        rewrite_wav<2, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        throw std::logic_error("unsupported 16-bit IEEE_FLOAT format");
                    }
                    else {
                        throw std::logic_error("unsupported 16-bit format");
                    }
                    break;
                case 3:
                    if(fmt == WaveFormat::PCM) {
                        rewrite_wav<3, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        throw std::logic_error("unsupported 24-bit IEEE_FLOAT format");
                    }
                    else {
                        throw std::logic_error("unsupported 24-bit format");
                    }
                    break;
                case 4:
                    if(fmt == WaveFormat::PCM) {
                        rewrite_wav<4, WaveFormat::PCM>(dir, source, dest, reader, f);
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        rewrite_wav<4, WaveFormat::IEEE_FLOAT>(dir, source, dest, reader, f);
                    }
                    else {
                        throw std::logic_error("unsupported 32-bit format");
                    }
                    break;
                case 8:
                    if(fmt == WaveFormat::PCM) {
                        throw std::logic_error("unsupported 64-bit PCM format");
                    }
                    else if(fmt == WaveFormat::IEEE_FLOAT) {
                        rewrite_wav<8, WaveFormat::IEEE_FLOAT>(dir, source, dest, reader, f);
                    }
                    else {
                        throw std::logic_error("unsupported 64-bit format");
                    }
                    break;
                default:
                    throw std::logic_error("unsupported xx-bit format");;
            }
        }

    template<typename WAVRead, typename SR>
    InterlacedBuffer readReverb(SR sample_rate, ResampleSincStats& stats, WAVRead & reader)
    {
        return {reader, sample_rate, stats};
    }

} // namespace imajuscule::audio
