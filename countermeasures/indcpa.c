#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "params.h"
#include "indcpa.h"
#include "polyvec.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "symmetric.h"
#include "randombytes.h"

int16_t CFP_dec_S;
int16_t CFP_dec_U;
int16_t CFP_dec_V;
int16_t CFP_enc_R;
int16_t CFP_enc_A;
int16_t CFP_enc_E1;
int16_t CFP_enc_E2;
int16_t CFP_enc_M;
int16_t CFP_enc_T;
int16_t CFP_key_A;
int16_t CFP_key_S;
int16_t CFP_key_E;
const int16_t CFPDec = 0;

/*************************************************
* Name:        pack_pk
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r: pointer to the output serialized public key
*              polyvec *pk: pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  polyvec_tobytes(r, pk);
  memcpy(r+KYBER_POLYVECBYTES, seed, KYBER_SYMBYTES);
}

/*************************************************
* Name:        unpack_pk
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk: pointer to output public-key polynomial vector
*              - uint8_t *seed: pointer to output seed to generate matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  polyvec_frombytes(pk, packedpk);
  memcpy(seed, packedpk+KYBER_POLYVECBYTES, KYBER_SYMBYTES);
}

/*************************************************
* Name:        pack_sk
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r: pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void pack_sk(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes(r, sk);
}

/*************************************************
* Name:        unpack_sk
*
* Description: De-serialize the secret key; inverse of pack_sk
*
* Arguments:   - polyvec *sk: pointer to output vector of polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk(polyvec *sk, const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk: pointer to the input vector of polynomials b
*              poly *v: pointer to the input polynomial v
**************************************************/
static void pack_ciphertext(uint8_t r[KYBER_INDCPA_BYTES], polyvec *b, poly *v)
{
  polyvec_compress(r, b);
  poly_compress(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b: pointer to the output vector of polynomials b
*              - poly *v: pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext(polyvec *b, poly *v, const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress(b, c);
  poly_decompress(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int32_t *r: pointer to output buffer
*              - unsigned int len: requested number of 16-bit integers (uniform mod q)
*              - const uint8_t *buf: pointer to input buffer (assumed to be uniformly random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
static unsigned int rej_uniform(int32_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;

    if(val0 < KYBER_Q)
      r[ctr++] = val0;
    if(ctr < len && val1 < KYBER_Q)
      r[ctr++] = val1;
  }

  return ctr;
}

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)

/*************************************************
* Name:        gen_matrix
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a: pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed: boolean deciding whether A or A^T is generated
**************************************************/
#if(XOF_BLOCKBYTES % 3)
#error "Implementation of gen_matrix assumes that XOF_BLOCKBYTES is a multiple of 3"
#endif

#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES];
  xof_state state;

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed)
        xof_absorb(&state, seed, i, j);
      else
        xof_absorb(&state, seed, j, i);

      xof_squeezeblocks(buf, GEN_MATRIX_NBLOCKS, &state);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        xof_squeezeblocks(buf, 1, &state);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}


static void hide(poly* v, int16_t t, uint8_t rand){
  int32_t hide = (int32_t)t * (int32_t)KYBER_Q;
  for (int i = 0; i < KYBER_N; i++){
    v->coeffs[i] += hide;
  }

  poly_pqt_reduce(v, rand);
}

static int32_t lift(int32_t p, int32_t q, uint8_t t){
    return pqt_montgomery_reduce((int64_t)p * 4137947LL + (int64_t)q * 20631166LL, t);
}

static int32_t csubp(int32_t r){
  return r - (KYBER_P & ((1 << 31) - (((KYBER_P - r - 1) >> 31) & 1)));
}
static int32_t caddp(int32_t r){
  return r + (KYBER_P & ((1 << 31) - ((r >> 31) & 1)));
}
/*************************************************
* Name:        indcpa_keypair_derand
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
*                             (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
*              - const uint8_t *coins: pointer to input randomness
*                             (of length KYBER_SYMBYTES bytes)
**************************************************/
void indcpa_keypair_derand(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                           uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES],
                           const uint8_t coins[KYBER_SYMBYTES])
{  
  int32_t rand;
  randombytes((uint8_t *)&rand, 4);
  CFP_key_A = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_key_S = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_key_E = p_barrett_reduce(rand);
  uint8_t hide_mod;
  randombytes((uint8_t *)&hide_mod, 1);
  hide_mod = hide_mod & 7;
  uint16_t hide_val;
  randombytes((uint8_t *)&hide_val, 2);
  unsigned int i, j, k;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;
  poly polyA;
  
  memcpy(buf, coins, KYBER_SYMBYTES);
  buf[KYBER_SYMBYTES] = KYBER_K;
  hash_g(buf, buf, KYBER_SYMBYTES+1);

  gen_a(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&skpv.vec[i], noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(&e.vec[i], noiseseed, nonce++);


  for (j = 0u; j < KYBER_N; j++){
    polyA.coeffs[j] = CFP_key_A;
  }
  
  poly_p_ntt(&polyA);

  
  for (i = 0u; i < KYBER_K; i++)
  {
    hide(&e.vec[i], hide_val, hide_mod);
    hide(&skpv.vec[i], hide_val, hide_mod);
    for (k = 0u; k < KYBER_K; k++){
      hide(&a[k].vec[i], hide_val, hide_mod);
    }

    for (j = 0u; j < KYBER_N; j++){
      skpv.vec[i].coeffs[j] = lift(CFP_key_S, skpv.vec[i].coeffs[j], hide_mod);
      e.vec[i].coeffs[j] = lift(CFP_key_E, e.vec[i].coeffs[j], hide_mod);
      
      for (k = 0u; k < KYBER_K; k++){
        a[k].vec[i].coeffs[j] = lift(polyA.coeffs[j], a[k].vec[i].coeffs[j], hide_mod);
      }
    }
  }

  polyvec_pqt_ntt(&skpv, hide_mod);
  polyvec_pqt_ntt(&e, hide_mod);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_pqt_basemul_acc_montgomery(&pkpv.vec[i], &a[i], &skpv, hide_mod);
    poly_pqt_tomont(&pkpv.vec[i], hide_mod);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_pqt_reduce(&pkpv, hide_mod);
  for (i = 0; i < KYBER_K; i++){
    for (j = 0; j < KYBER_N; j++){
      e.vec[i].coeffs[j] = skpv.vec[i].coeffs[j];
      a[0].vec[i].coeffs[j] = pkpv.vec[i].coeffs[j];
    }
  }

  polyvec_p_reduce(&e);
  polyvec_p_invntt_tomont(&e);
  polyvec_p_reduce(&a[0]);
  polyvec_p_invntt_tomont(&a[0]);
  int32_t valpk = p_barrett_reduce((int32_t)KYBER_K * (int32_t)CFP_key_A * (int32_t)CFP_key_S);
  for (i = 0; i < KYBER_N; i++){
    for (j = 0; j < KYBER_K; j++){
      if (p_montgomery_reduce(a[0].vec[j].coeffs[i]) != csubp(CFP_key_E + p_barrett_reduce( valpk * ((i << 1) - 254))) ||
          p_montgomery_reduce(e.vec[j].coeffs[i]) != CFP_key_S) {
          return;
      }
    }
  }
  
  polyvec_q_reduce(&pkpv);
  polyvec_q_reduce(&skpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}


/*************************************************
* Name:        indcpa_enc
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c: pointer to output ciphertext
*                            (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m: pointer to input message
*                                  (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk: pointer to input public key
*                                   (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins used as seed
*                                      (of length KYBER_SYMBYTES) to deterministically
*                                      generate all randomness
**************************************************/
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  int32_t rand;
  randombytes((uint8_t *)&rand, 4);
  CFP_enc_R = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_key_A = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_enc_E1 = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_enc_E2 = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_enc_M = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_enc_T = p_barrett_reduce(rand);
  uint8_t hide_mod;
  randombytes((uint8_t *)&hide_mod, 1);
  hide_mod = hide_mod & 7;
  uint16_t hide_val;
  randombytes((uint8_t *)&hide_val, 2);
  unsigned int i, j, y;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], b;
  poly v, k, epp, polyA;

  unpack_pk(&pkpv, seed, pk);
  poly_frommsg(&k, m);
  gen_at(at, seed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2(&epp, coins, nonce++);

  
  hide(&epp, hide_val, hide_mod);
  hide(&k, hide_val, hide_mod);
  for (j = 0u; j < KYBER_N; j++){
    polyA.coeffs[j] = CFP_enc_A;
    v.coeffs[j] = CFP_enc_T;
    epp.coeffs[j] = lift(CFP_enc_E2, epp.coeffs[j], hide_mod);
    k.coeffs[j] = lift(CFP_enc_M, k.coeffs[j], hide_mod);
  }
  
  poly_p_ntt(&polyA);
  poly_p_ntt(&v);

  for (i = 0u; i < KYBER_K; i++)
  {
    hide(&pkpv.vec[i], hide_val, hide_mod);
    hide(&sp.vec[i], hide_val, hide_mod);
    hide(&ep.vec[i], hide_val, hide_mod);
    for (j = 0u; j < KYBER_N; j++){
      pkpv.vec[i].coeffs[j] = lift(v.coeffs[j], pkpv.vec[i].coeffs[j], hide_mod);
      sp.vec[i].coeffs[j] = lift(CFP_enc_R, sp.vec[i].coeffs[j], hide_mod);
      ep.vec[i].coeffs[j] = lift(CFP_enc_E1, ep.vec[i].coeffs[j], hide_mod);
      
      for (y = 0u; y < KYBER_K; y++){
        at[y].vec[i].coeffs[j] = lift(polyA.coeffs[j], at[y].vec[i].coeffs[j], hide_mod);
      }
    }
  }

  polyvec_pqt_ntt(&sp, hide_mod);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pqt_basemul_acc_montgomery(&b.vec[i], &at[i], &sp, hide_mod);

  polyvec_pqt_basemul_acc_montgomery(&v, &pkpv, &sp, hide_mod);

  polyvec_pqt_invntt_tomont(&b, hide_mod);
  poly_pqt_invntt_tomont(&v, hide_mod);

  polyvec_add(&b, &b, &ep);
  poly_add(&v, &v, &epp);
  poly_add(&v, &v, &k);
  polyvec_pqt_reduce(&b, hide_mod);
  poly_pqt_reduce(&v, hide_mod);

  int32_t valu = p_barrett_reduce((int32_t)KYBER_K * (int32_t)CFP_enc_A * (int32_t)CFP_enc_R);
  int32_t valv = p_barrett_reduce((int32_t)KYBER_K * (int32_t)CFP_enc_T * (int32_t)CFP_enc_R);

  for (i = 0; i < KYBER_N; i++){
    if (p_barrett_reduce(v.coeffs[i]) != csubp(csubp(CFP_enc_M + CFP_enc_E2 + p_barrett_reduce( valv * ((i << 1) - 254))))){
            
            return;
    }
    for (j = 0; j < KYBER_K; j++){
      if (p_barrett_reduce(b.vec[j].coeffs[i]) != csubp(CFP_enc_E1 + p_barrett_reduce( valu * ((i << 1) - 254)))){
          return;
      }
    }
  }
  
  polyvec_q_reduce(&b);
  poly_q_reduce(&v);

  pack_ciphertext(c, &b, &v);
}

/*************************************************
* Name:        indcpa_dec
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m: pointer to output decrypted message
*                            (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c: pointer to input ciphertext
*                                  (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  int32_t rand;
  randombytes((uint8_t *)&rand, 4);
  CFP_dec_S = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_dec_U = p_barrett_reduce(rand);
  randombytes((uint8_t *)&rand, 4);
  CFP_dec_V = p_barrett_reduce(rand);
  uint8_t hide_mod;
  randombytes((uint8_t *)&hide_mod, 1);
  hide_mod = hide_mod & 7;
  uint16_t hide_val;
  randombytes((uint8_t *)&hide_val, 2);
  polyvec b, skpv;
  poly v, mp;

  unpack_ciphertext(&b, &v, c);
  unpack_sk(&skpv, sk);

  
  hide(&v, hide_val, hide_mod);

  for (size_t j = 0u; j < KYBER_N; j++){
    mp.coeffs[j] = CFP_dec_S;
    v.coeffs[j] = lift(CFP_dec_V, v.coeffs[j], hide_mod);
  }
  
  poly_p_ntt(&mp);

  for (size_t i = 0u; i < KYBER_K; i++) {
    hide(&skpv.vec[i], hide_val, hide_mod);
    hide(&b.vec[i], hide_val, hide_mod);
    
    for (size_t j = 0u; j < KYBER_N; j++){
      skpv.vec[i].coeffs[j] = lift(mp.coeffs[j], skpv.vec[i].coeffs[j], hide_mod);
      b.vec[i].coeffs[j] = lift(CFP_dec_U, b.vec[i].coeffs[j], hide_mod);
    }
  }

  
  polyvec_pqt_ntt(&b, hide_mod);
  polyvec_pqt_basemul_acc_montgomery(&mp, &skpv, &b, hide_mod);
  poly_pqt_invntt_tomont(&mp, hide_mod);

  poly_sub(&mp, &v, &mp);
  poly_pqt_reduce(&mp, hide_mod);

  int32_t val = p_barrett_reduce((int32_t)KYBER_K * (int32_t)CFP_dec_S * (int32_t)CFP_dec_U);
  for (int i = 0; i < KYBER_N; i++){
      if (p_barrett_reduce(mp.coeffs[i]) != caddp(CFP_dec_V - p_barrett_reduce( val * ((i << 1) - 254)))){
              
            return;
      }
  }

  poly_q_reduce(&mp);
  poly_tomsg(m, &mp);
}
