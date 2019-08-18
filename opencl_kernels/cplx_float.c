
struct cplx_float {
  float real;
  float imag;
};

inline struct cplx_float complexFromReal(float r) {
  return (struct cplx_float) {
    .real = r,
    .imag = 0.f
  };
}

inline struct cplx_float polar(float const theta) {
  struct cplx_float c;
  c.imag = sincos(theta,&c.real);
  return c;
}

inline struct cplx_float cplxMult(struct cplx_float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = b.real * a.real - b.imag * a.imag,
    .imag = b.real * a.imag + b.imag * a.real
  };
}

inline struct cplx_float cplxSub(struct cplx_float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = a.real - b.real,
    .imag = a.imag - b.imag
  };
}

inline struct cplx_float cplxAdd(struct cplx_float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = a.real + b.real,
    .imag = a.imag + b.imag
  };
}

inline struct cplx_float cplxConjugate(struct cplx_float const a) {
  return (struct cplx_float) {
    .real = a.real,
    .imag = -a.imag
  };
}

inline struct cplx_float cplxScalarMult(float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = b.real * a,
    .imag = b.imag * a
  };
}

inline struct cplx_float cplxScalarSub(float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = a - b.real,
    .imag = - b.imag
  };
}

inline struct cplx_float cplxScalarAdd(float const a, struct cplx_float const b) {
  return (struct cplx_float) {
    .real = a + b.real,
    .imag = b.imag
  };
}

////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on local memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign(__local struct cplx_float * a, struct cplx_float const b) {
  a->real += b.real;
  a->imag += b.imag;
}

inline void butterfly(__local struct cplx_float *v, int const i, const struct cplx_float twiddle) {
  struct cplx_float const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on global memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign_global(__global struct cplx_float * a, struct cplx_float const b) {
  a->real += b.real;
  a->imag += b.imag;
}
inline void butterfly_global(__global struct cplx_float *v, int i,  struct cplx_float twiddle) {
  struct cplx_float const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign_global(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly out of place
////////////////////////////////////////////////////////////////////

inline void butterfly_outofplace(int const idx,
                               int const idxD,
                               __local struct cplx_float const *from,
                               __local struct cplx_float *to,
                               int const i,
                               int const Ns,
                               const struct cplx_float twiddle) {
  struct cplx_float const t = cplxMult(from[idx+i], twiddle);
  struct cplx_float const fi = from[idx];
  to[idxD+Ns] = cplxSub(fi, t);
  to[idxD]    = cplxAdd(fi, t);
}
