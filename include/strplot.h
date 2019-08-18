/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {

    static constexpr auto delimiter_char = ' ';
    static constexpr auto background_char = ' ';
    static constexpr auto horizontal_line_char = '-';
    static constexpr auto default_curve_char = '+';

    struct StringPlot;
    struct StringPlot {
        template<typename T>
        static constexpr T ignoredValue() {
          static_assert(std::is_floating_point_v<T>);
          return std::numeric_limits<T>::quiet_NaN();
        }

        template<typename T>
        static constexpr bool isIgnoredValue(T val) {
          return std::isnan(val);
        }

        // if range is empty, the first curve drawn defines the range
        StringPlot(int Height, int Width, range<float> r = {}) :
        Height(Height),
        Width(Width),
        p(Height),
        range_(std::move(r))
        {
            std::fill(p.begin(), p.end(), std::string(Width, background_char));
        }

        template<typename T>
        void drawLog(T container, char c = default_curve_char, bool include_zero = false) {
            std::transform(container.begin(),
                           container.end(),
                           container.begin(),
                           [](auto v){
                               if( isIgnoredValue(v) ) {
                                   return v;
                               }
                               if(v < 0.f) {
                                   throw std::logic_error("cannot draw log of negative values");
                               }
                               return std::log(v);
                           });
            draw(container, c, include_zero);
        }

        template<typename T>
        void draw(T const & container, char c = default_curve_char, bool include_zero = true) {
            if(container.empty()) {
                //throw "cannot draw empty container";
                return;
            }

            if(range_.empty()) {
                if(include_zero) {
                    range_.extend(0.f);
                }
                for(auto v : container) {
                    if(isIgnoredValue(v)) {
                        continue;
                    }
                    range_.extend(v);
                }
            }

            for(auto i=0; i<Width; ++i) {
                auto fi = static_cast<float>(.5f + (container.size()-1)) / (Width-1);
                auto v_index = static_cast<size_t>(fi * i);
                if(static_cast<int>(container.size()) == Width) {
                    assert(v_index == i);
                }

                auto val = static_cast<float>(container[v_index]);

                if(std::isnan(val)) {
                    constexpr auto margin = 1;
                    for(int j=margin; j<Height-margin; ++j) {
                        p[j][i] = 'X';
                    }
                    continue;
                }

                auto height = val_to_height(val);

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

        void writeToStream(std::ostream & stream) {

            using namespace std;
            auto zero_height = val_to_height(0.f);
            make_line(Height-1-zero_height);

            auto bar = std::string(Width, delimiter_char);

            stream << bar << endl;

            auto height = 0;
            for(auto const & str : *this) {
                stream << str;

                if(height == zero_height) {
                    stream << " 0";
                }
                else if(height == 0 ) {
                    stream << " " << range_.getMax();
                }
                else if(height==Height-1) {
                    stream << " " << range_.getMin();
                }
                stream << endl;

                height++;
            }

            stream << bar << endl;
        }

        void log() {
            std::cout << *this;
        }

    private:
        int Height;
        int Width;
        std::vector<std::string> p;
        range<float> range_;

        template<typename T>
        int val_to_height(T val) const {
            if(range_.empty() || range_.delta() == 0.f) {
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

        int invert_height(int h) const {
            return Height - 1 - h;
        }

        friend std::ostream& operator<<(std::ostream& stream, StringPlot& obj)
        {
            obj.writeToStream(stream);
            return stream;
        }

    };

}
