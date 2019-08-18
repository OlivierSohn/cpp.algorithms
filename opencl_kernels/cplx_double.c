
struct cplx_double {
  double real;
  double imag;
};

inline struct cplx_double complexFromReal(double r) {
  return (struct cplx_double) {
    .real = r,
    .imag = 0.
  };
}

inline struct cplx_double polar(double const theta) {
  struct cplx_double c;
  c.imag = sincos(theta,&c.real);
  return c;
}

inline struct cplx_double cplxMult(struct cplx_double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = b.real * a.real - b.imag * a.imag,
    .imag = b.real * a.imag + b.imag * a.real
  };
}

inline struct cplx_double cplxSub(struct cplx_double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = a.real - b.real,
    .imag = a.imag - b.imag
  };
}

inline struct cplx_double cplxAdd(struct cplx_double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = a.real + b.real,
    .imag = a.imag + b.imag
  };
}

inline struct cplx_double cplxConjugate(struct cplx_double const a) {
  return (struct cplx_double) {
    .real = a.real,
    .imag = -a.imag
  };
}

inline struct cplx_double cplxScalarMult(double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = b.real * a,
    .imag = b.imag * a
  };
}

inline struct cplx_double cplxScalarSub(double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = a - b.real,
    .imag = - b.imag
  };
}

inline struct cplx_double cplxScalarAdd(double const a, struct cplx_double const b) {
  return (struct cplx_double) {
    .real = a + b.real,
    .imag = b.imag
  };
}

////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on local memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign(__local struct cplx_double * a, struct cplx_double const b) {
  a->real += b.real;
  a->imag += b.imag;
}

inline void butterfly(__local struct cplx_double *v, int const i, const struct cplx_double twiddle) {
  struct cplx_double const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly on global memory
////////////////////////////////////////////////////////////////////

inline void cplxAddAssign_global(__global struct cplx_double * a, struct cplx_double const b) {
  a->real += b.real;
  a->imag += b.imag;
}
inline void butterfly_global(__global struct cplx_double *v, int i,  struct cplx_double twiddle) {
  struct cplx_double const t = cplxMult(v[i], twiddle);
  v[i] = cplxSub(v[0], t);
  cplxAddAssign_global(&v[0], t);
}


////////////////////////////////////////////////////////////////////
// Functions used when performing the butterfly out of place
////////////////////////////////////////////////////////////////////

inline void butterfly_outofplace(int const idx,
                               int const idxD,
                               __local struct cplx_double const *from,
                               __local struct cplx_double *to,
                               int const i,
                               int const Ns,
                               const struct cplx_double twiddle) {
  struct cplx_double const t = cplxMult(from[idx+i], twiddle);
  struct cplx_double const fi = from[idx];
  to[idxD+Ns] = cplxSub(fi, t);
  to[idxD]    = cplxAdd(fi, t);
}
