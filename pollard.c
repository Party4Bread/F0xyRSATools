//
// Created by solsa on 2018-07-16.
//

#include "stdafx.h"


bool pollardrho(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit){
    void g(mpz_t src){
        mpz_mul(src,src,src);
        mpz_add_ui(src,src,1);
        mpz_mod(src,src,n);
    }
    mpz_t x, y, d, a;
    unsigned long c=1;
    mpz_init_set_ui(x,2);
    mpz_init_set_ui(y,2);
    mpz_init_set_ui(d,1);
    mpz_init(a);
    while(mpz_cmp_ui(d,1)==0) {
        g(x);
        g(y);
        g(y);
        mpz_sub(a,x,y);
        mpz_abs(a,a);
        mpz_gcd(d,a,n);
        if(iterlimit!=INF_ITER){
            iterlimit--;
            if(iterlimit<0)break;
        }
    }
    if(mpz_cmp_ui(d,1)==0||mpz_cmp(d,n)==0){
        return false;
    }
    else{
        mpz_set(q,d);
        mpz_div(p,n,d);
        mpz_clears(x, y, d, a,NULL);
        return true;
    }
}

//TODO:Increase Memory Efficiency, Increase Speed, iterlimit
bool pollardrho_brent(mpz_t p, mpz_t q, const mpz_t n,unsigned long long iterlimit){
    mpz_t y,c,m,g,r,v,x,k,ys,upper , i,t1;
    mpz_init_set_ui(y,1);
    mpz_init_set_ui(c,1);
    mpz_init_set_ui(m,1);
    mpz_init_set_ui(g,1);
    mpz_init_set_ui(r,1);
    mpz_init_set_ui(v,1);
    mpz_inits(x,k,ys,upper,i,t1,NULL);

    void f(mpz_t src){
        mpz_mul(src,src,src);
        mpz_add_ui(src,src,c);
        mpz_mod(src,src,n);
    }

    while(mpz_cmp_ui(g,1)==0&&iterlimit>0){
        iterlimit--;
        mpz_set(x,y);
        for(mpz_set_ui(i,0); mpz_cmp(i,r)<0;mpz_add_ui(i,i,1)){
            f(y);
        }
        mpz_set_ui(k,0);
        while(mpz_cmp(k,r)<0&&mpz_cmp_ui(g,1)==0){
            mpz_set(ys,y);
            mpz_sub(t1,r,k);
            if(mpz_cmp(m,t1)<0){
                mpz_set(upper,m);
            }else{
                mpz_set(upper,t1);
            }
            for(mpz_set_ui(i,0); mpz_cmp(i,upper)<0;mpz_add_ui(i,i,1)){
                f(y);
                mpz_sub(t1,x,y);
                mpz_abs(t1,t1);
                mpz_mul(v,t1,v);
                mpz_mod(v,v,n);
            }
            mpz_gcd(g,v,n);
            mpz_add(k,k,m);
        }
        mpz_mul_ui(r,r,2);
    }
    if(mpz_cmp(g,n)==0){

        mpz_set_ui(g, 1);
        while (mpz_cmp_ui(g,1)==0)
        {
            f(ys);
            mpz_sub(t1,x,ys);
            mpz_abs(t1,t1);
            mpz_gcd(g,t1,n);
        }
    }
    mpz_set(q,g);
    mpz_div(p,n,g);
    mpz_clears(y,c,m,g,r,v,x,k,ys,upper , i,t1,NULL);
    return true;
}




//bool pollardrho(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit,long long c){
//    mpz_t xF, cS, x ,f, cnt, t1;
//    mpz_init_set_ui(xF,2);
//    mpz_init_set_ui(cS,2);
//    mpz_init_set_ui(x,2);
//    mpz_init_set_ui(f,1);
//    mpz_init(cnt);
//    mpz_init(t1);
//
//    while(mpz_cmp_ui(f,1)==0&&(iterlimit==INF_ITER||mpz_cmp_ui(cS,iterlimit)==-1)){
//        mpz_set_ui(cnt,1);
//        while(mpz_cmp(cnt,cS)<1&&mpz_cmp_ui(f,1)<1){
//            mpz_mul(x,x,x);
//            mpz_add_ui(x,x,c);
//            mpz_mod(x,x,n);
//
//            mpz_sub(t1,x,xF);
//            mpz_gcd(f,t1,n);
//            mpz_add_ui(cnt,cnt,1);
//        }
//        mpz_mul_ui(cS,cS,2);
//        mpz_set(xF,x);
//    }
//    if(mpz_cmp_ui(f,1)!=0){
//        mpz_div(p,n,f);
//        mpz_set(q,f);
//        return true;
//    }
//    else{
//        return false;
//    }
//}
