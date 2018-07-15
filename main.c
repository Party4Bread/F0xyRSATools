#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>
#define INF_ITER -1337


bool fermat(mpz_t p, mpz_t q, mpz_t n, long long iterlimit);
bool kraitchik(mpz_t p, mpz_t q, mpz_t n, long long iterlimit);
bool wiener(mpz_t p, mpz_t q, mpz_t n, mpz_t e);
bool smallprime(mpz_t p, mpz_t q, mpz_t n, mpz_t e);
bool pollardrho(mpz_t p, mpz_t q, mpz_t n, long long iterlimit,long long c);
bool pollardrho_brent(mpz_t p, mpz_t q, mpz_t n, long long iterlimit,long long c);
void test_factorization();

int main()
{
    test_factorization();
    return 0;
}

void test_factorization(){
    mpz_t N;
    mpz_t p, q;

    mpz_init(N);
    mpz_init(p);
    mpz_init(q);

    mpz_set_str(N, "4E733FEBB94DB17CA3E6AA26EC33B4960C150C52300E06C60B3318F0744FEF2D687A8F5BF598894A22EEC4ABDAE01B197E4CC5603DE67EB670E261EB4E4CC5E26241EDCDE494CCE415BBC5A410ABCEFDFF6199BBCDF62E9D434FAA88A1D16012520F80D126208206FF80191E20ED7423CDCE5B8A555B4161534E789A74F0A701", 16);

    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    puts("-----------Fermat Method-----------");
    fermat(p, q, N, INF_ITER);
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);

    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    puts("-----------Kraitchik Method-----------");
    kraitchik(p, q, N, INF_ITER);
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);

    mpz_set_str(N, "21", 10);

    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    puts("-----------Pollard Rho Method-----------");
    pollardrho(p, q, N, INF_ITER, 1);
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
}
bool fermat(mpz_t p, mpz_t q, mpz_t n, long long iterlimit){
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

bool kraitchik(mpz_t p, mpz_t q, mpz_t n, long long int iterlimit){
    mpz_t x,y,k,t1,t2,t3,t4;

    mpz_init(x);
    mpz_init(y);
    mpz_init(k);
    mpz_init(t1);
    mpz_init(t2);
    mpz_init(t3);
    mpz_init(t4);

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
                    mpz_clear(x);
                    mpz_clear(y);
                    mpz_clear(k);
                    mpz_clear(t1);
                    mpz_clear(t2);
                    mpz_clear(t3);
                    mpz_clear(t4);
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

//Not Working WTF
bool pollardrho(mpz_t p, mpz_t q, mpz_t n, long long iterlimit,long long c){
    mpz_t xF, cS, x ,f, cnt, t1;
    mpz_init_set_ui(xF,2);
    mpz_init_set_ui(cS,2);
    mpz_init_set_ui(x,2);
    mpz_init_set_ui(f,1);
    mpz_init(cnt);
    mpz_init(t1);

    while(mpz_cmp_ui(f,1)==0&&(iterlimit==INF_ITER||mpz_cmp_ui(cS,iterlimit)==-1)){
        mpz_set_ui(cnt,1);
        while(mpz_cmp(cnt,cS)<1&&mpz_cmp_ui(f,1)<1){
            mpz_mul(x,x,x);
            mpz_add_ui(x,x,c);
            mpz_mod(x,x,n);

            mpz_sub(t1,x,xF);
            mpz_gcd(f,t1,n);
            mpz_add_ui(cnt,cnt,1);
        }
        mpz_mul_ui(cS,cS,2);
        mpz_set(xF,x);
    }
    if(mpz_cmp_ui(f,1)!=0){
        mpz_div(q,n,f);
        mpz_set(p,f);
        return true;
    }
    else{
        return false;
    }
}
