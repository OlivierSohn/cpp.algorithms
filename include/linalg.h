namespace imajuscule {

  
  /* Returns a square identity matrix of size 'n'.
   */
  template<typename T>
  auto mkIdentity(int n) {
    Matrix<T> res;
    res.resize(n,n);
    for(int i=0; i<n; ++i) {
      res[i][i] = 1;
    }
    return res;
  }
  
  // https://en.wikipedia.org/wiki/Crout_matrix_decomposition
  template<typename T>
  bool crout(const Matrix<T> & A, Matrix<T> & L, Matrix<T> & U) {
    int const n = A.countRows();
    assert(n==A.countColumns());
    assert(n==L.countColumns());
    assert(n==U.countColumns());
    assert(n==L.countRows());
    assert(n==U.countRows());

    T sum = 0;
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        U[i][j] = (i==j)?1:0;
      }
    }
    
    for (int j = 0; j < n; j++) {
      for (int i = j; i < n; i++) {
        sum = 0;
        for (int k = 0; k < j; k++) {
          sum = sum + L[i][k] * U[k][j];
        }
        L[i][j] = A[i][j] - sum;
      }
      
      for (int i = j; i < n; i++) {
        sum = 0;
        for(int k = 0; k < j; k++) {
          sum = sum + L[j][k] * U[k][i];
        }
        if (L[j][j] == 0) {
          // det(L) close to 0, can't divide by 0...
          return false;
        }
        U[j][i] = (A[j][i] - sum) / L[j][j];
      }
    }
    return true;
  }
  
  template<typename T>
  void multiply(Matrix<T> const & A, Matrix<T> const & B, Matrix<T> & C) {
    assert(C.countRows() == A.countRows());
    assert(C.countColumns() == B.countColumns());
    int const nRows = C.countRows();
    int const nCols = C.countColumns();

    int const nMults = A.countColumns();

    for(int a=0; a<nRows; ++a) {
      for(int b=0; b<nCols; ++b) {
        T sum = 0;
        for(int i = 0; i<nMults; ++i) {
          sum += A[a][i] * B[i][b];
        }
        C[a][b] = sum;
      }
    }
  }

  template<typename T>
  bool equals(Matrix<T> const & A, Matrix<T> const & B, T const epsilon) {
    assert(A.countRows() == B.countRows());
    assert(A.countColumns() == B.countColumns());
    int const nRows = A.countRows();
    int const nCols = A.countColumns();
    for(int a=0; a<nRows; ++a) {
      for(int b=0; b<nCols; ++b) {
        if(std::abs(A[a][b] - B[a][b]) > epsilon) {
          return false;
        }
      }
    }
    return true;
  }

  template<typename T>
  bool isUpper(Matrix<T> const & A) {
    int const n = A.countRows();
    assert(n==A.countColumns());

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if(j >= i) {
          continue;
        }
        if(A[i][j] != 0) {
          return false;
        }
      }
    }
    return true;
  }
  
  template<typename T>
  bool isLower(Matrix<T> const & A) {
    int const n = A.countRows();
    assert(n==A.countColumns());
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if(j <= i) {
          continue;
        }
        if(A[i][j] != 0) {
          return false;
        }
      }
    }
    return true;
  }
  
  // solves 'Ax = b' where A is square.
  //
  // TODO optimize by passing allocated matrixes L,U and y in a struct.
  template<typename T>
  bool solve(Matrix<T> const & A, std::vector<T> const & b, std::vector<T> & x) {
    assert(b.size() == A.countRows());

    x.resize(A.countColumns());
    int const sz = A.countColumns();
    
    // LU decomposition
    Matrix<T> L,U;
    L.resize(sz, sz);
    U.resize(sz, sz);
    
    if(!crout(A,L,U)) {
      return false;
    }
    
    // A.x = b
    // L.U.x = b
    //
    // we first solve for y where y = U.x:
    //
    // L.y = b

    std::vector<T> y;
    y.resize(sz);

    for(int i=0; i<sz; ++i) {
      T sum = 0;
      for(int j=0; j<i; ++j) {
        sum += L[i][j] * y[j];
      }
      assert(L[i][i] != 0);
      y[i] = (b[i] - sum) / L[i][i];
    }
    
    // now solve for x:
    // U.x = y
    
    for(int i=sz-1; i>=0; --i) {
      T sum = 0;
      for(int j=i+1; j<sz; ++j) {
        sum += U[i][j] * x[j];
      }
      assert(U[i][i] != 0);
      x[i] = (y[i] - sum) / U[i][i];
    }

    return true;
  }
}
