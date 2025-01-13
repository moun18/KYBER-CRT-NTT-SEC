#include <stdint.h>
#include "params.h"
#include "reduce.h"

const int32_t KYBER_PQT[8] = {
  0x1CF53113, 0x19E8DB11, 0x16DC850F, 0x13D02F0D, 0x10C3D90B, 0xDB78309, 0xAAB2D07, 0x79ED705
};

const int32_t PQTINV[8] = {
  0xD05A411B, 0x6FB75F1, 0x7E9C9EF, 0xF57037C5, 0xADB32AA3, 0xB7DAFB39, 0xEC62B0B7, 0x17BD5DCD 
};

const int32_t vPQT[8] = {
  9, 10, 11, 13, 15, 19, 24, 34
};

/*************************************************
* Name:        q_montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q, where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-q2^15,...,q2^15-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int16_t q_montgomery_reduce(int32_t a)
{
  int32_t r;

  r = (int32_t)a*QINV;
  r = (a - (int64_t)r*KYBER_Q) >> 32;
  return r - (KYBER_Q & ((1 << 31) - (((KYBER_Q - r - 1) >> 31) & 1)));
}

/*************************************************
* Name:        p_montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q, where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-p2^15,...,p2^15-1}
*
* Returns:     integer in {-p+1,...,p-1} congruent to a * R^-1 modulo q.
**************************************************/
int16_t p_montgomery_reduce(int32_t a)
{
  int32_t r;

  r = (int32_t)a*PINV;
  r = (a - (int64_t)r*KYBER_P) >> 32;
  return r - (KYBER_P & ((1 << 31) - (((KYBER_P - r - 1) >> 31) & 1)));
}

/*************************************************
* Name:        pqt_montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              32-bit integer congruent to a * R^-1 mod q, where R=2^32
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-pqt2^31,...,pqt2^31-1}
               - uint8_t t: index of the moduli chosen.
*
* Returns:     integer in {-pqt+1,...,pqt-1} congruent to a * R^-1 modulo pqt.
**************************************************/
int32_t pqt_montgomery_reduce(int64_t a, uint8_t t)
{
  int32_t r;

  r = (int32_t)a*PQTINV[t];
  r = (a - (int64_t)r*KYBER_PQT[t]) >> 32;
  return r - (KYBER_PQT[t] & ((1 << 31) - (((KYBER_PQT[t] - r - 1) >> 31) & 1)));
}

/*************************************************
* Name:        q_barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
**************************************************/
int16_t q_barrett_reduce(int32_t a) {
  int32_t r;
  const int32_t v = ((1ULL<<32) + KYBER_Q/2)/KYBER_Q;

  r  = ((int64_t)v*a) >> 32;
  r *= KYBER_Q;
  r = a - r;
  return r - (KYBER_Q & ((1 << 31) - (((KYBER_Q - r - 1) >> 31) & 1)));
}

/*************************************************
* Name:        p_barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              centered representative congruent to a mod p in {-(p-1)/2,...,(p-1)/2}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {-(p-1)/2,...,(p-1)/2} congruent to a modulo p.
**************************************************/
int16_t p_barrett_reduce(int32_t a) {
  int32_t r;
  const int32_t v = ((1ULL<<32) + KYBER_P/2)/KYBER_P;

  r  = ((int64_t)v*a) >> 32;
  r *= KYBER_P;
  r = a - r;
  return r - (KYBER_P & ((1 << 31) - (((KYBER_P - r - 1) >> 31) & 1)));
}

/*************************************************
* Name:        pqt_barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              centered representative congruent to a mod pqt in {-(pqt-1)/2,...,(pqt-1)/2}
*
* Arguments:   - int16_t a: input integer to be reduced
               - uint8_t t: index of the moduli chosen.
*
* Returns:     integer in {-(pqt-1)/2,...,(pqt-1)/2} congruent to a modulo pqt.
**************************************************/
int32_t pqt_barrett_reduce(int32_t a, uint8_t t) {
  int32_t r;
  r  = ((((uint64_t)vPQT[t])*(uint64_t)a) >> 32);
  r *= KYBER_PQT[t];
  r = a - r;
  return r - (KYBER_PQT[t] & ((1 << 31) - (((KYBER_PQT[t] - r - 1) >> 31) & 1)));
}