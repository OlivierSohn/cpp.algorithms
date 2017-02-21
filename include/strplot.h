namespace imajuscule {
    
    static constexpr auto delimiter_char = ' ';
    static constexpr auto background_char = ' ';
    static constexpr auto horizontal_line_char = '-';
    static constexpr auto default_curve_char = '+';
    
    template<int Height, int Width>
    struct StringPlot {
        
        // if range is empty, the first curve drawn defines the range
        StringPlot(range<float> r = {}) : range_(std::move(r)) {
            std::fill(p.begin(), p.end(), std::string(Width, background_char));
        }
        
        template<typename T>
        void draw(T const & container, char c = default_curve_char) {
            
            auto r = range_;
            if(r.empty()) {
                auto minmax = std::minmax_element (container.begin(),container.end());
                r.set(*minmax.first, *minmax.second);
                r.extend(0.f);
                range_ = r;
            }
            
            for(auto i=0; i<Width; ++i) {
                auto fi = static_cast<float>(.5f + (container.size()-1)) / (Width-1);
                auto v_index = static_cast<size_t>(fi * i);
                if(container.size() == Width) {
                    assert(v_index == i);
                }
                
                auto height = val_to_height(container[v_index]);
                
                if(height < 0) {
                    p[0][i] = 'E';
                }
                else if(height >= Height) {
                    p[Height-1][i] = 'E';
                }
                else {
                    p[height][i] = c;
                }
            }
        }
        
        auto begin() const { return p.begin(); }
        auto end() const { return p.end(); }
        
        void log() {
            using namespace std;
            auto zero_height = val_to_height(0.f);
            make_line(Height-1-zero_height);
            
            auto bar = std::string(Width, delimiter_char);
            
            cout << bar << endl;
            
            auto height = 0;
            for(auto const & str : *this) {
                cout << str;
                
                if(height == zero_height) {
                    cout << " 0";
                }
                else if(height == 0 ) {
                    cout << " " << range_.getMax();
                }
                else if(height==Height-1) {
                    cout << " " << range_.getMin();
                }
                cout << endl;
                
                height++;
            }

            cout << bar << endl;
        }
        
    private:
        std::array<std::string, Height> p;
        range<float> range_;
        
        template<typename T>
        int val_to_height(T val) const {
            if(range_.delta() == 0.f) {
                return 0;
            }
            float normalized_val = (val - range_.getMin()) / range_.delta();
            
            // each line should contain an equal interval of the range,
            // except the first and last which contain half the equal interval
            // (to have a better aspect for sinus-like curves)
            
            float plot_val = normalized_val * (Height-1) + 0.5f;
            
            return invert_height(static_cast<int>(plot_val));
        }
        
        void make_line(int h) {
            auto index = invert_height(h);
            if(index < 0 || index >= p.size()) {
                return;
            }
            auto &s = p[index];
            std::replace( s.begin(), s.end(), ' ', horizontal_line_char);
        }
        
        constexpr int invert_height(int h) const {
            return Height - 1 - h;
        }
    };
}

