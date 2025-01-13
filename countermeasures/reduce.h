#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT 1353 // 2^32 mod q 
#define QINV 1806234369 // q^-1 mod 2^32
#define PINV 2340676097 // p^-1 mod 2^32
extern const int32_t PQTINV[8];
extern const int32_t vPQT[8];


#define q_montgomery_reduce KYBER_NAMESPACE(q_montgomery_reduce)
int16_t q_montgomery_reduce(int32_t a);

#define p_montgomery_reduce KYBER_NAMESPACE(p_montgomery_reduce)
int16_t p_montgomery_reduce(int32_t a);

#define pqt_montgomery_reduce KYBER_NAMESPACE(pqt_montgomery_reduce)
int32_t pqt_montgomery_reduce(int64_t a, uint8_t t);

#define q_barrett_reduce KYBER_NAMESPACE(q_barrett_reduce)
int16_t q_barrett_reduce(int32_t a);

#define p_barrett_reduce KYBER_NAMESPACE(p_barrett_reduce)
int16_t p_barrett_reduce(int32_t a);

#define pqt_barrett_reduce KYBER_NAMESPACE(pqt_barrett_reduce)
int32_t pqt_barrett_reduce(int32_t a, uint8_t t);

#endif
