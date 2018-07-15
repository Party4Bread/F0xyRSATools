//
// Created by solsa on 2018-07-16.
//


#include "stdafx.h"

bool shanks(mpz_t p, mpz_t q, const mpz_t n, long long int iterlimit){
    mpz_t a, f, h1, h2, k, pp, ppp, qq, qqq, qqqq, r, te;
    mpz_inits(a, f, h1, h2, k, pp, ppp, qq, qqq, qqqq, r, te, NULL);
    int i = 0;

    mpz_sqrt(k,n);

    mpz_set(a,k);
    mpz_set(h1,k);
    mpz_set_ui(h1,1);
    mpz_set_ui(ppp,0);
    mpz_set_ui(qqq,1);
    mpz_set(qqqq,n);
    mpz_set_ui(r,0);

    while(true) {
        mpz_sub(pp,k,r);
        mpz_sub(qq,ppp,pp);
        mpz_mul(qq,qq,a);
        mpz_add(qq,qq,qqqq);
        mpz_add(a,pp,k);
        mpz_fdiv_qr(a,r,a,qq);
        mpz_set(te,h2);
        mpz_addmul(te,a,h1);
        mpz_set(h2,h1);
        mpz_set(h1,te);
        mpz_set(ppp,pp);
        mpz_set(qqqq,qqq);
        mpz_set(qqq,qq);
        mpz_sqrt(te,qq);
        if ( (++i%2) != 0 || !mpz_perfect_square_p(qq)) continue;
        mpz_sub(te,h2,te);
        mpz_gcd(f,te,n);
        if (mpz_cmp_ui(f,1)>0 && mpz_cmp(f,n)<0) {
            mpz_set(p,f);
            mpz_div(q,n,f);
            return true;
        }
    }
}