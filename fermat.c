//
// Created by solsa on 2018-07-16.
//


#include "stdafx.h"

bool fermat(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit){
    mpz_t a, b;

    mpz_init(a);
    mpz_init(b);

    mpz_sub_ui(a, n, 1);
    mpz_sqrt(a, a);
    mpz_add_ui(a, a, 1);

    while (true)
    {
        mpz_mul(b, a, a);
        mpz_sub(b, b, n);
        if (mpz_perfect_square_p(b))
            break;

        mpz_add_ui(a, a, 1);
        if(iterlimit!=INF_ITER){
            iterlimit--;
            if(iterlimit<0)return false;
        }
    }

    mpz_sqrt(p, b);
    mpz_sqrt(q, b);
    mpz_add(p, a, p);
    mpz_sub(q, a, q);

    mpz_clear(a);
    mpz_clear(b);
    return true;
}

//bool fermat2(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit){
//    mpz_t x,y,w;
//    mpz_inits(x,y,w,NULL);
//    mpz_sub_ui(x, n, 1);
//    mpz_sqrt(x,x);
//    mpz_add_ui(x,x,1);
//
//    mpz_pow_ui(y,x,2);
//    mpz_sub(y,y,n);
//    mpz_sqrt(y,y);
//    while(iterlimit==INF_ITER||iterlimit>0){
//        mpz_mul(w,x,x);
//        mpz_submul(w,y,y);
//        int r = mpz_cmp(w,n);
//        if(r>0){
//            mpz_add_ui(y,y,1);
//        }else if(r<0){
//            mpz_add_ui(x,x,1);
//        }else{
//            mpz_add(p, x, y);
//            mpz_sub(q, x, y);
//            mpz_clears(x,y,w,NULL);
//            return true;
//        }
//        if(iterlimit!=INF_ITER)iterlimit--;
//    }
//    return false;
//}

