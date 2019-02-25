
struct cplx {
  float real;
  float imag;
};

inline struct cplx complexFromReal(float r) {
  return (struct cplx) {
    .real = r,
    .imag = 0.f
  };
}

inline struct cplx polar(float const theta) {
  struct cplx c;
  c.imag = sincos(theta,&c.real);
  return c;
}

inline struct cplx cplxMult(struct cplx const a, struct cplx const b) {
  return (struct cplx) {
    .real = b.real * a.real - b.imag * a.imag,
    .imag = b.real * a.imag + b.imag * a.real
  };
}

inline struct cplx cplxSub(struct cplx const a, struct cplx const b) {
  return (struct cplx) {
    .real = a.real - b.real,
    .imag = a.imag - b.imag
  };
}

inline struct cplx cplxAdd(struct cplx const a, struct cplx const b) {
  return (struct cplx) {
    .real = a.real + b.real,
    .imag = a.imag + b.imag
  };
}

inline struct cplx cplxConjugate(struct cplx const a) {
  return (struct cplx) {
    .real = a.real,
    .imag = -a.imag
  };
}

inline struct cplx cplxScalarMult(float const a, struct cplx const b) {
  return (struct cplx) {
    .real = b.real * a,
    .imag = b.imag * a
  };
}

inline struct cplx cplxScalarSub(float const a, struct cplx const b) {
  return (struct cplx) {
    .real = a - b.real,
    .imag = - b.imag
  };
}

inline struct cplx cplxScalarAdd(float const a, struct cplx const b) {
  return (struct cplx) {
    .real = a + b.real,
    .imag = b.imag
  };
}

////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on local memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign(__local struct cplx * a, struct cplx const b) {
  a->real += b.real;
  a->imag += b.imag;
}

inline void butterfly(__local struct cplx *v, int const i, const struct cplx twiddle) {
  struct cplx const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on global memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign_global(__global struct cplx * a, struct cplx const b) {
  a->real += b.real;
  a->imag += b.imag;
}
inline void butterfly_global(__global struct cplx *v, int i,  struct cplx twiddle) {
  struct cplx const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign_global(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly out of place
////////////////////////////////////////////////////////////////////

inline void butterfly_outofplace(int const idx,
                               int const idxD,
                               __local struct cplx const *from,
                               __local struct cplx *to,
                               int const i,
                               int const Ns,
                               const struct cplx twiddle) {
  struct cplx const t = cplxMult(from[idx+i], twiddle);
  struct cplx const fi = from[idx];
  to[idxD+Ns] = cplxSub(fi, t);
  to[idxD]    = cplxAdd(fi, t);
}
