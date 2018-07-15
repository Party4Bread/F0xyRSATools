//
// Created by solsa on 2018-07-16.
//

#include "stdafx.h"

//TODO:Need speed optimization on kraitchik
bool kraitchik(mpz_t p, mpz_t q, const mpz_t n, long long int iterlimit){
    mpz_t x,y,k,t1,t2,t3,t4;

    mpz_inits(x,y,k,t1,t2,t3,t4,NULL);

    mpz_sqrt(x,n);
    mpz_add_ui(x,x,1);

    while(true){
        mpz_set_ui(k,1);
        while(true){
            mpz_mul(t1,x,x);
            mpz_mul(t2,k,n);
            mpz_sub(t1,t1,t2);
            if(mpz_cmp_ui(t1,0)<0)break;
            if(mpz_perfect_square_p(t1)){
                mpz_sqrt(y,t1);
                mpz_add(t1,x,y);
                mpz_sub(t2,x,y);
                mpz_mod(t3,t1,n);
                mpz_mod(t4,t2,n);
                if(mpz_cmp_ui(t3,0)!=0&&mpz_cmp_ui(t4,0)!=0){
                    mpz_gcd(p,t1,n);
                    mpz_gcd(q,t2,n);
                    mpz_clears(x,y,k,t1,t2,t3,t4,NULL);
                    return true;
                }
                if(iterlimit!=INF_ITER){
                    iterlimit--;
                    if(iterlimit<0)return false;
                }
            }
            mpz_add_ui(k,k,1);
        }
        mpz_add_ui(x,x,1);
    }
}
