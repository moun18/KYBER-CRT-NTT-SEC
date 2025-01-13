#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*coeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int32_t coeffs[KYBER_N];
} poly;

#define poly_compress KYBER_NAMESPACE(poly_compress)
void poly_compress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], const poly *a);
#define poly_decompress KYBER_NAMESPACE(poly_decompress)
void poly_decompress(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES]);

#define poly_tobytes KYBER_NAMESPACE(poly_tobytes)
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], const poly *a);
#define poly_frombytes KYBER_NAMESPACE(poly_frombytes)
void poly_frombytes(poly *r, const uint8_t a[KYBER_POLYBYTES]);

#define poly_frommsg KYBER_NAMESPACE(poly_frommsg)
void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]);
#define poly_tomsg KYBER_NAMESPACE(poly_tomsg)
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *r);

#define poly_getnoise_eta1 KYBER_NAMESPACE(poly_getnoise_eta1)
void poly_getnoise_eta1(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_getnoise_eta2 KYBER_NAMESPACE(poly_getnoise_eta2)
void poly_getnoise_eta2(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_p_ntt KYBER_NAMESPACE(poly_p_ntt)
void poly_p_ntt(poly *r);
#define poly_p_invntt_tomont KYBER_NAMESPACE(poly_p_invntt_tomont)
void poly_p_invntt_tomont(poly *r);

#define poly_pqt_ntt KYBER_NAMESPACE(poly_pqt_ntt)
void poly_pqt_ntt(poly *r, uint8_t t);
#define poly_pqt_invntt_tomont KYBER_NAMESPACE(poly_pqt_invntt_tomont)
void poly_pqt_invntt_tomont(poly *r, uint8_t t);

#define poly_pqt_basemul_montgomery KYBER_NAMESPACE(poly_basemul_montgomery)
void poly_pqt_basemul_montgomery(poly *r, const poly *a, const poly *b, uint8_t t);
#define poly_tomont KYBER_NAMESPACE(poly_tomont)
void poly_pqt_tomont(poly *r, uint8_t t);

#define poly_q_reduce KYBER_NAMESPACE(poly_q_reduce)
void poly_q_reduce(poly *r);

#define poly_p_reduce KYBER_NAMESPACE(poly_p_reduce)
void poly_p_reduce(poly *r);

#define poly_pqt_reduce KYBER_NAMESPACE(poly_pqt_reduce)
void poly_pqt_reduce(poly *r, uint8_t t);

#define poly_add KYBER_NAMESPACE(poly_add)
void poly_add(poly *r, const poly *a, const poly *b);
#define poly_sub KYBER_NAMESPACE(poly_sub)
void poly_sub(poly *r, const poly *a, const poly *b);

#endif
