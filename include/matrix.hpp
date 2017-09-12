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
    
    static inline MatrixCoord indexToCoordinates(int i, int nColumns) {
        auto rowIdx = i/nColumns;
        auto colIdx = i - rowIdx * nColumns;
        return {rowIdx, colIdx};
    }
    
    static inline int coordinatesToIndex(MatrixCoord coord, int nColumns) {
        return nColumns * coord.row + coord.col;
    }
    
    template<typename T>
    struct Matrix {
        
        Matrix(Matrix const & other) = default;
        Matrix & operator = (Matrix const & other) = delete;
        
        Matrix(Matrix && other) = default;
        Matrix & operator = (Matrix && other) = default;
        
        Matrix() = default;
        
        Matrix(int rows, int columns) :
        nRows(rows)
        , nColumns(columns) {
            resize(rows, columns);
        }
        
        void resize(int rows, int columns) {
            nRows = rows;
            nColumns = columns;
            v.resize(rows*columns);
        }
        int countRows() const { return nRows; }
        int countColumns() const { return nColumns; }
        T * operator [] (int row) { return v.data() + row * nColumns; }
        T const * operator [] (int row) const { return v.data() + row * nColumns; }
        
        T const & getByIndex(int i) const { return v[i]; }
        T & editByIndex(int i) { return v[i]; }

        T const & get(MatrixCoord const & c) const { return v[coordinatesToIndex(c, nColumns)]; }
        T & edit(MatrixCoord const & c) { return v[coordinatesToIndex(c, nColumns)]; }
        
        MatrixCoord indexToCoordinates(int i) const {
            return imajuscule::indexToCoordinates(i, nColumns);
        }
        
        auto begin() const { return v.begin(); }
        auto end() const { return v.end(); }
        auto begin() { return v.begin(); }
        auto end() { return v.end(); }
        
    private:
        int nRows, nColumns;
        std::vector<T> v;
    };

} // NS imajuscule

