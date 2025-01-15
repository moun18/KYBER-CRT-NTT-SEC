#include <stdint.h>
#include "params.h"
#include "ntt.h"
#include "reduce.h"
#include "randombytes.h"

/* Code to generate zetas and zetas_inv used in the number-theoretic transform:

#define KYBER_ROOT_OF_UNITY 17

static const uint8_t tree[128] = {
  0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120,
  4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
  2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
  6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
  1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121,
  5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
  3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123,
  7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127
};

void init_ntt() {
  unsigned int i;
  int16_t tmp[128];

  tmp[0] = MONT;
  for(i=1;i<128;i++)
    tmp[i] = fqmul(tmp[i-1],MONT*KYBER_ROOT_OF_UNITY % KYBER_Q);

  for(i=0;i<128;i++) {
    zetas[i] = tmp[tree[i]];
    if(zetas[i] > KYBER_Q/2)
      zetas[i] -= KYBER_Q;
    if(zetas[i] < -KYBER_Q/2)
      zetas[i] += KYBER_Q;
  }
}
*/

const int16_t p_zetas[128] = {
  5569, 4279, 2332, 876, 4466, 953, 4350, 1028, 3838, 7186, 1843, 3907, 5486, 3207, 5144, 4620, 721, 4500, 4, 792, 3196, 2966, 3512, 4086, 2523, 289, 
  3455, 481, 3066, 269, 7176, 7544, 3598, 5752, 2108, 2610, 2153, 3839, 7384, 2642, 808, 6364, 388, 14, 2772, 3505, 2700, 4611, 6620, 4990, 4852, 571,
  5524, 3050, 4782, 2073, 3361, 4912, 4770, 7378, 1454, 3695, 1915, 2801, 1566, 2828, 6912, 1358, 49, 2021, 746, 1769, 4617, 127, 2103, 1620, 5839, 3972,
  2994, 1375, 3415, 242, 1830, 1333, 2780, 5089, 1411, 2862, 5963, 5481, 2217, 1149, 4753, 4012, 3233, 2611, 2351, 4638, 4285, 3520, 5670, 1234, 6221, 2798,
  972, 431, 847, 6405, 825, 2049, 6290, 1098, 2336, 1668, 7662, 3919, 181, 5114, 6361, 7475, 5298, 4388, 871, 3476, 4639, 4483, 4319, 2571
};

const int32_t pq_zetas[128] = {
  24769113, 19391123, 12268889, 13911167, 23323982, 22268172, 23154884, 8857221, 10288697, 14316889, 25387548, 20796374, 22649074, 21440878, 20282984, 2616160,
  25010057, 23124310, 3886590, 16176978, 10126754, 16294367, 23983594, 17808644, 14788448, 1129396, 19659134, 16046090, 594503, 18342497, 24786082, 21260871, 24736418,
  6235043, 4188253, 11378171, 8589511, 10711153, 9024878, 22623187, 5154759, 6719558, 23926703, 16122433, 22738532, 22009570, 20426479, 11433939, 5905628, 13062690, 10927234,
  8380542, 573918, 15971849, 19076705, 22853048, 23729970, 4851623, 480992, 9831377, 8896052, 15695978, 23375198, 8997252, 14779810, 22730907, 13740540, 7014111, 4961975,
  10110217, 21046686, 1883614, 7316929, 14494174, 17445654, 17852264, 15060599, 15097137, 9919165, 1030629, 23599447, 19440853, 22269049, 2090565, 2975327, 20935814, 9379912,
  23322378, 8808389, 24116140, 6469619, 20163774, 22817323, 13491848, 20649761, 14926794, 5302241, 16457340, 10212334, 11517339, 12932793, 15816413, 13724487, 9004930, 
  1191527, 15961549, 9164280, 23095491, 21023722, 13297860, 7241792, 9771330, 22046806, 2390459, 22850956, 11786573, 1067840, 10566489, 7088243, 16475539, 11403902, 1886233,
  17098777, 12031922, 7985198, 7132451, 20789105, 22914994
};

const uint32_t fPQT[8] = {
  0x861939E, 0x861939E, 0xCF414A1, 0x3CF129B, 0x5553D9C, 0xCF414A1, 0x5553D9C, 0x5553D9C
};

/*************************************************
* Name:        fpmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int16_t a: first factor
*              - int16_t b: second factor
*
* Returns 16-bit integer congruent to a*b*R^{-1} mod p
**************************************************/
static int16_t fpmul(int16_t a, int16_t b) {
  return p_montgomery_reduce((int32_t)a*b);
}

/*************************************************
* Name:        fpqtmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int32_t a: first factor
*              - int32_t b: second factor
*              - uint8_t t: index of the moduli chosen.
*
* Returns 32-bit integer congruent to a*b*R^{-1} mod pqt
**************************************************/
static int32_t fpqtmul(int32_t a, int32_t b, uint8_t t) {
  return pqt_montgomery_reduce((int64_t)a*b, t);
}

static int16_t p_zetas_power(uint8_t index);

static int32_t pq_zetas_power(uint8_t index);

static void p_BlindedNTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list);

static void pqt_BlindedNTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list, const uint8_t t);

static void p_BlindedINTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list);

static void pqt_BlindedINTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list, const uint8_t t);

uint8_t bit_reverse(uint8_t num) {
  uint8_t i;
  uint8_t revnum = 0u;

  for(i=0; i<7; ++i) {
    revnum |= ((num & (1 << i)) >> i) << (6u - i);
  }

  return revnum;
}

/*************************************************
* Name:        p_ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rp.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zp
**************************************************/
void p_ntt(int32_t r[256]) {

  uint8_t mask_list[128]; 

  for(uint16_t i = 0; i < (1 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }
  
  for(uint16_t i = (7 << (7 - BLOCK_SIZE)); i < (8 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }

  randombytes(&mask_list[1 << (7 - BLOCK_SIZE)], 6 << (7 - BLOCK_SIZE));

  p_BlindedNTT_internal(r, BLOCK_SIZE, mask_list);
}

/*************************************************
* Name:        p_invntt_tomont
*
* Description: Inplace inverse number-theoretic transform in Rp and
*              multiplication by Montgomery factor 2^16.
*              Input is in bitreversed order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zp
**************************************************/
void p_invntt(int32_t r[256]) {

  uint16_t j;

  const int16_t f = 4124; // mont^2/128

  uint8_t mask_list[128]; 

  for(uint16_t i = 0; i < (1 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }
  
  for(uint16_t i = (7 << (7 - BLOCK_SIZE)); i < (8 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }

  randombytes(&mask_list[1 << (7 - BLOCK_SIZE)], 6 << (7 - BLOCK_SIZE));

  p_BlindedINTT_internal(r, BLOCK_SIZE, mask_list);

  for(j = 0; j < 256; j++)
    r[j] = fpmul(r[j], f);
}

/*************************************************
* Name:        pqt_ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rpqt.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zpqt
*              - uint8_t t: index of the moduli chosen.
**************************************************/
void pqt_ntt(int32_t r[256], uint8_t t) {

  uint8_t mask_list[128]; 

  for(uint16_t i = 0; i < (1 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }
  
  for(uint16_t i = (7 << (7 - BLOCK_SIZE)); i < (8 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }

  randombytes(&mask_list[1 << (7 - BLOCK_SIZE)], 6 << (7 - BLOCK_SIZE));

  pqt_BlindedNTT_internal(r, BLOCK_SIZE, mask_list, t);
}

/*************************************************
* Name:        pqt_invntt_tomont
*
* Description: Inplace inverse number-theoretic transform in Rpqt and
*              multiplication by Montgomery factor 2^16.
*              Input is in bitreversed order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zpqt
*              - uint8_t t: index of the moduli chosen.
**************************************************/
void pqt_invntt(int32_t r[256], uint8_t t) {

  uint16_t j;

  const int32_t f = fPQT[t]; // mont^2/128

  uint8_t mask_list[128]; 

  for(uint16_t i = 0; i < (1 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }
  
  for(uint16_t i = (7 << (7 - BLOCK_SIZE)); i < (8 << (7 - BLOCK_SIZE)); i++){
    mask_list[i] = 0;
  }

  randombytes(&mask_list[1 << (7 - BLOCK_SIZE)], 6 << (7 - BLOCK_SIZE));

  pqt_BlindedINTT_internal(r, BLOCK_SIZE, mask_list, t);

  for(j = 0; j < 256; j++)
    r[j] = fpqtmul(r[j], f, t);
}

/*************************************************
* Name:        pqt_basemul
*
* Description: Multiplication of polynomials in Zpqt[X]/(X^2-zeta)
*              used for multiplication of elements in Rpqt in NTT domain
*
* Arguments:   - int32_t r[2]: pointer to the output polynomial
*              - const int32_t a[2]: pointer to the first factor
*              - const int32_t b[2]: pointer to the second factor
*              - int32_t zeta: integer defining the reduction polynomial
*              - uint8_t t: index of the moduli chosen.
**************************************************/
void pqt_basemul(int32_t r[2], const int32_t a[2], const int32_t b[2], int32_t zeta, uint8_t t)
{
  int32_t mem;

  mem = fpqtmul(a[0], b[0], t);
  r[0]  = fpqtmul(a[1], b[1], t);
  r[1]  = fpqtmul(a[0] + a[1], b[0] + b[1], t);
  r[1] -= r[0];
  r[1] -= mem;
  r[0]  = fpqtmul(r[0], zeta, t);
  r[0] += mem;
}

int16_t p_zetas_power(uint8_t index)
{
  return -(p_zetas[index & 127] & ((1 << 15) - (index >> 7))) + (p_zetas[index & 127] & (~((1 << 15) - (index >> 7))));
}

int32_t pq_zetas_power(uint8_t index)
{
  return -(pq_zetas[index & 127] & ((1 << 31) - (index >> 7))) + (pq_zetas[index & 127] & (~((1 << 31) - (index >> 7))));
}

void p_BlindedNTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list)
{
  uint16_t len, start, i, j, block_num, block_pad, pos1, pos2;
  int16_t omega1, omega2, omega3, omega4, data;
  uint8_t k1, k2, block_counter;
  k1 = 1;
  block_num = 1 << (7 - block_size);
  block_pad = 0;

  for(len = 7; len > block_size; len--){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << len){
      k2 = bit_reverse(k1);
      for(j = start; j < (start + (1 << len)); j += 1 << block_size){
        pos1 = block_counter >> block_size;
        pos2 = (pos1 >> (len + 1 - block_size)) << (len + 1 - block_size);
        pos2 += pos1 & ((1 << (len - block_size)) - 1);
        pos2 |= block_pad;
        omega1 = p_zetas_power((256 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        pos2 += 1 << (len - block_size);
        pos2 &= block_num - 1;
        pos2 |= block_pad;
        omega2 = p_zetas_power((256 + k2 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        for(i = j; i < (j + (1 << block_size)); i++){
          data = fpmul(omega2, r[i + (1 << len)]);
          r[i] = fpmul(omega1, r[i]) + data;
          r[i + (1 << len)] = r[i] - 2 * data;
          block_counter++;
        }
      }
      k1++;
    }
    block_pad += block_num;
  }

  for(len = block_size; len >= 1; len--){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << block_size){
      pos1 = (block_counter >> (block_size + 1)) << 1;
      pos1 |= block_pad;
      pos2 = pos1 + block_num + (((block_counter & ((1 << block_size) - 1)) >> (len - 1)) & 1);
      omega1 = p_zetas_power((256 + mask_list[pos2] - mask_list[pos1]) & 255);
      omega3 = p_zetas_power((256 + mask_list[pos2 + 1] - mask_list[pos1]) & 255);
      for(i = start; i < (start + (2 << block_size)); i += 2 << len){
        k2 = bit_reverse(k1);
        omega2 = p_zetas_power((k2 + 256 + mask_list[pos2] - mask_list[pos1 + 1]) & 255);
        omega4 = p_zetas_power((k2 + 256 + mask_list[pos2 + 1] - mask_list[pos1 + 1]) & 255);
        for(j = i; j < (i + (1 << (len - 1))); j++){
          data = fpmul(omega2, r[j + (1 << len)]);
          r[j] = fpmul(omega1, r[j]) + data;
          r[j + (1 << len)] = r[j] - 2 * data;
          ++block_counter;

          data = fpmul(omega4, r[j + (3 << (len - 1))]);
          r[j + (1 << (len - 1))] = fpmul(omega3, r[j + (1 << (len - 1))]) + data;
          r[j + (3 << (len - 1))] = r[j + (1 << (len - 1))] - 2 * data;
          ++block_counter;
        }
        k1++;
      }
    }
    block_pad += block_num;
  }
}

void p_BlindedINTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list)
{
  uint16_t len, start, i, j, block_num, block_pad, pos1, pos2;
  int16_t omega1, omega2, omega3, omega4, data;
  uint8_t k1, k2, block_counter;
  k1 = 127;
  block_num = 1 << (7-block_size);
  block_pad = 6 << (7-block_size);
 
  for(len = 1; len <= block_size; len++){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << block_size){
      pos1 = (block_counter >> (block_size + 1)) << 1;
      pos1 |= block_pad;
      pos2 = pos1 + block_num + (((block_counter & ((1 << block_size) - 1)) >> (len - 1)) & 1);
      omega1 = p_zetas_power((256 + mask_list[pos2] - mask_list[pos1]) & 255);
      omega3 = p_zetas_power((256 + mask_list[pos2 + 1] - mask_list[pos1]) & 255);
      for(i = start; i < (start + (2 << block_size)); i += 2 << len){
        k2 = bit_reverse(k1);
        omega2 = p_zetas_power((k2 + 128 + 256 + mask_list[pos2] - mask_list[pos1 + 1]) & 255);
        omega4 = p_zetas_power((k2 + 128 + 256 + mask_list[pos2 + 1] - mask_list[pos1 + 1]) & 255);
        for(j = i; j < (i + (1 << (len - 1))); j++){
          data = r[j];
          r[j] = fpmul(omega1, data + r[j + (1 << len)]);
          r[j + (1 << len)] = fpmul(omega2, data - r[j + (1 << len)]);
          block_counter++;

          data = r[j + (1 << (len - 1))];
          r[j + (1 << (len - 1))] = fpmul(omega3, r[j + (3 << (len - 1))] + data);
          r[j + (3 << (len - 1))] = fpmul(omega4, data - r[j + (3 << (len - 1))]);
          block_counter++;
        }
        k1--;
      }
    }
    block_pad -= block_num;
  }

  for(len = block_size + 1; len <= 7; len++){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << len){
      k2 = bit_reverse(k1);
      for(j = start; j < (start + (1 << len)); j += 1 << block_size){
        pos1 = block_counter >> block_size;
        pos2 = (pos1 >> (len + 1 - block_size)) << (len + 1 - block_size);
        pos2 += pos1 & ((1 << (len - block_size)) - 1);
        pos2 |= block_pad;
        omega1 = p_zetas_power((256 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        pos2 += 1 << (len - block_size);
        pos2 &= block_num - 1;
        pos2 |= block_pad;
        omega2 = p_zetas_power((256 + 128 + k2 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        for(i = j; i < (j + (1 << block_size)); i++){
          data = r[i];
          r[i] = fpmul(omega1, r[i + (1 << len)] + data);
          r[i + (1 << len)] = fpmul(omega2, data - r[i + (1 << len)]);
          ++block_counter;
        }
      }
      k1--;
    }
    block_pad -= block_num;
  }
}


void pqt_BlindedNTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list, uint8_t t)
{
  uint16_t len, start, i, j, block_num, block_pad, pos1, pos2;
  int32_t omega1, omega2, omega3, omega4, data;
  uint8_t k1, k2, block_counter;
  k1 = 1;
  block_num = 1 << (7 - block_size);
  block_pad = 0;

  for(len = 7; len > block_size; len--){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << len){
      k2 = bit_reverse(k1);
      for(j = start; j < (start + (1 << len)); j += 1 << block_size){
        pos1 = block_counter >> block_size;
        pos2 = (pos1 >> (len + 1 - block_size)) << (len + 1 - block_size);
        pos2 += pos1 & ((1 << (len - block_size)) - 1);
        pos2 |= block_pad;
        omega1 = pq_zetas_power((256 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        pos2 += 1 << (len - block_size);
        pos2 &= block_num - 1;
        pos2 |= block_pad;
        omega2 = pq_zetas_power((256 + k2 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        for(i = j; i < (j + (1 << block_size)); i++){
          data = fpqtmul(omega2, r[i + (1 << len)], t);
          r[i] = fpqtmul(omega1, r[i], t) + data;
          r[i + (1 << len)] = r[i] - 2 * data;
          block_counter++;
        }
      }
      k1++;
    }
    block_pad += block_num;
  }

  for(len = block_size; len >= 1; len--){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << block_size){
      pos1 = (block_counter >> (block_size + 1)) << 1;
      pos1 |= block_pad;
      pos2 = pos1 + block_num + (((block_counter & ((1 << block_size) - 1)) >> (len - 1)) & 1);
      omega1 = pq_zetas_power((256 + mask_list[pos2] - mask_list[pos1]) & 255);
      omega3 = pq_zetas_power((256 + mask_list[pos2 + 1] - mask_list[pos1]) & 255);
      for(i = start; i < (start + (2 << block_size)); i += 2 << len){
        k2 = bit_reverse(k1);
        omega2 = pq_zetas_power((k2 + 256 + mask_list[pos2] - mask_list[pos1 + 1]) & 255);
        omega4 = pq_zetas_power((k2 + 256 + mask_list[pos2 + 1] - mask_list[pos1 + 1]) & 255);
        for(j = i; j < (i + (1 << (len - 1))); j++){
          data = fpqtmul(omega2, r[j + (1 << len)], t);
          r[j] = fpqtmul(omega1, r[j], t) + data;
          r[j + (1 << len)] = r[j] - 2 * data;
          ++block_counter;

          data = fpqtmul(omega4, r[j + (3 << (len - 1))], t);
          r[j + (1 << (len - 1))] = fpqtmul(omega3, r[j + (1 << (len - 1))], t) + data;
          r[j + (3 << (len - 1))] = r[j + (1 << (len - 1))] - 2 * data;
          ++block_counter;
        }
        k1++;
      }
    }
    block_pad += block_num;
  }
}

void pqt_BlindedINTT_internal(int32_t r[256], const uint16_t block_size, const uint8_t* mask_list, uint8_t t)
{
  uint16_t len, start, i, j, block_num, block_pad, pos1, pos2;
  int32_t omega1, omega2, omega3, omega4, data;
  uint8_t k1, k2, block_counter;
  k1 = 127;
  block_num = 1 << (7-block_size);
  block_pad = 6 << (7-block_size);
 
  for(len = 1; len <= block_size; len++){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << block_size){
      pos1 = (block_counter >> (block_size + 1)) << 1;
      pos1 |= block_pad;
      pos2 = pos1 + block_num + (((block_counter & ((1 << block_size) - 1)) >> (len - 1)) & 1);
      omega1 = pq_zetas_power((256 + mask_list[pos2] - mask_list[pos1]) & 255);
      omega3 = pq_zetas_power((256 + mask_list[pos2 + 1] - mask_list[pos1]) & 255);
      for(i = start; i < (start + (2 << block_size)); i += 2 << len){
        k2 = bit_reverse(k1);
        omega2 = pq_zetas_power((k2 + 128 + 256 + mask_list[pos2] - mask_list[pos1 + 1]) & 255);
        omega4 = pq_zetas_power((k2 + 128 + 256 + mask_list[pos2 + 1] - mask_list[pos1 + 1]) & 255);
        for(j = i; j < (i + (1 << (len - 1))); j++){
          data = r[j];
          r[j] = fpqtmul(omega1, data + r[j + (1 << len)], t);
          r[j + (1 << len)] = fpqtmul(omega2, data - r[j + (1 << len)], t);
          block_counter++;

          data = r[j + (1 << (len - 1))];
          r[j + (1 << (len - 1))] = fpqtmul(omega3, r[j + (3 << (len - 1))] + data, t);
          r[j + (3 << (len - 1))] = fpqtmul(omega4, data - r[j + (3 << (len - 1))], t);
          block_counter++;
        }
        k1--;
      }
    }
    block_pad -= block_num;
  }

  for(len = block_size + 1; len <= 7; len++){
    block_counter = 0;
    for(start = 0; start < 256; start += 2 << len){
      k2 = bit_reverse(k1);
      for(j = start; j < (start + (1 << len)); j += 1 << block_size){
        pos1 = block_counter >> block_size;
        pos2 = (pos1 >> (len + 1 - block_size)) << (len + 1 - block_size);
        pos2 += pos1 & ((1 << (len - block_size)) - 1);
        pos2 |= block_pad;
        omega1 = pq_zetas_power((256 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        pos2 += 1 << (len - block_size);
        pos2 &= block_num - 1;
        pos2 |= block_pad;
        omega2 = pq_zetas_power((256 + 128 + k2 + mask_list[pos1 | (block_pad + block_num)] - mask_list[pos2]) & 255);
        for(i = j; i < (j + (1 << block_size)); i++){
          data = r[i];
          r[i] = fpqtmul(omega1, r[i + (1 << len)] + data, t);
          r[i + (1 << len)] = fpqtmul(omega2, data - r[i + (1 << len)], t);
          ++block_counter;
        }
      }
      k1--;
    }
    block_pad -= block_num;
  }
}
