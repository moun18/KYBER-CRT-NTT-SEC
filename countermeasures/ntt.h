#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define p_zetas KYBER_NAMESPACE(p_zetas)
extern const int16_t p_zetas[128];
#define pq_zetas KYBER_NAMESPACE(pq_zetas)
extern const int32_t pq_zetas[128];

extern const uint32_t fPQT[8];

#define bit_reverse KYBER_NAMESPACE(bit_reverse)
uint8_t bit_reverse(uint8_t num);

#define p_ntt KYBER_NAMESPACE(p_ntt)
void p_ntt(int32_t poly[256]);

#define p_invntt KYBER_NAMESPACE(p_invntt)
void p_invntt(int32_t poly[256]);

#define pqt_ntt KYBER_NAMESPACE(pqt_ntt)
void pqt_ntt(int32_t poly[256], uint8_t t);

#define pqt_invntt KYBER_NAMESPACE(pqt_invntt)
void pqt_invntt(int32_t poly[256], uint8_t t);

#define pqt_basemul KYBER_NAMESPACE(pqt_basemul)
void pqt_basemul(int32_t r[2], const int32_t a[2], const int32_t b[2], int32_t zeta, uint8_t t);

#define BLOCK_SIZE (3)
#endif
