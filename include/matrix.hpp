/* Copyright (C) Olivier Sohn - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * Written by Olivier Sohn <olivier.sohn@gmail.com>, 2017
 */

namespace imajuscule {
    
    struct MatrixCoord {
        int row;
        int col;
    };
    
    template<typename T>
    struct Matrix {
        
        Matrix(Matrix const & other) = default;
        Matrix & operator = (Matrix const & other) = delete;
        
        Matrix(Matrix && other) = default;
        Matrix & operator = (Matrix && other) = default;
        
        Matrix() = default;
        
        Matrix(int rows, int columns) :
        n_rows(rows)
        , n_columns(columns) {
            resize(rows, columns);
        }
        
        void resize(int rows, int columns) {
            n_rows = rows;
            n_columns = columns;
            v.resize(rows*columns);
        }
        int countRows() const { return n_rows; }
        int countColumns() const { return n_columns; }
        T * operator [] (int row) { return v.data() + row * n_columns; }
        T const * operator [] (int row) const { return v.data() + row * n_columns; }
        
        T const & getByIndex(int i) const { return v[i]; }
        
        MatrixCoord indexToCoordinates(int i) const {
            auto rowIdx = i/n_columns;
            auto colIdx = i - rowIdx * n_columns;
            return {rowIdx, colIdx};
        }
        
        auto begin() const { return v.begin(); }
        auto end() const { return v.end(); }
        auto begin() { return v.begin(); }
        auto end() { return v.end(); }
        
    private:
        int n_rows, n_columns;
        std::vector<T> v;
    };

} // NS imajuscule

