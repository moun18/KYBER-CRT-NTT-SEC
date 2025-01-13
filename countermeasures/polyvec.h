#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

typedef struct{
  poly vec[KYBER_K];
} polyvec;

#define polyvec_compress KYBER_NAMESPACE(polyvec_compress)
void polyvec_compress(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], const polyvec *a);
#define polyvec_decompress KYBER_NAMESPACE(polyvec_decompress)
void polyvec_decompress(polyvec *r, const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]);

#define polyvec_tobytes KYBER_NAMESPACE(polyvec_tobytes)
void polyvec_tobytes(uint8_t r[KYBER_POLYVECBYTES], const polyvec *a);
#define polyvec_frombytes KYBER_NAMESPACE(polyvec_frombytes)
void polyvec_frombytes(polyvec *r, const uint8_t a[KYBER_POLYVECBYTES]);

#define polyvec_p_ntt KYBER_NAMESPACE(polyvec_p_ntt)
void polyvec_p_ntt(polyvec *r);
#define polyvec_p_invntt_tomont KYBER_NAMESPACE(polyvec_p_invntt_tomont)
void polyvec_p_invntt_tomont(polyvec *r);

#define polyvec_pqt_ntt KYBER_NAMESPACE(polyvec_pqt_ntt)
void polyvec_pqt_ntt(polyvec *r, uint8_t t);
#define polyvec_pqt_invntt_tomont KYBER_NAMESPACE(polyvec_pqt_invntt_tomont)
void polyvec_pqt_invntt_tomont(polyvec *r, uint8_t t);

#define polyvec_pqt_basemul_acc_montgomery KYBER_NAMESPACE(polyvec_pqt_basemul_acc_montgomery)
void polyvec_pqt_basemul_acc_montgomery(poly *r, const polyvec *a, const polyvec *b, uint8_t t);

#define polyvec_q_reduce KYBER_NAMESPACE(polyvec_q_reduce)
void polyvec_q_reduce(polyvec *r);

#define polyvec_p_reduce KYBER_NAMESPACE(polyvec_p_reduce)
void polyvec_p_reduce(polyvec *r);

#define polyvec_pqt_reduce KYBER_NAMESPACE(polyvec_pqt_reduce)
void polyvec_pqt_reduce(polyvec *r, uint8_t t);

#define polyvec_add KYBER_NAMESPACE(polyvec_add)
void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b);

#endif
