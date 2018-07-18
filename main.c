#include "stdafx.h"
#include <time.h>

bool fermat(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit);
bool kraitchik(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit);
bool wiener(mpz_t p, mpz_t q, const mpz_t n, mpz_t e);
bool smallprime(mpz_t p, mpz_t q, const mpz_t n, mpz_t e);
bool pollardrho(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit);
bool pollardrho_brent(mpz_t p, mpz_t q, const mpz_t n, long long iterlimit);
void test_factorization();

void speedtest();

int main(int argc, char **argv)
{
    if(argc!=2){
        printf("Usage : %s [number]",argv[0]);
        return 0;
    }
    mpz_t n,p,q;
    mpz_inits(p,q,NULL);
    mpz_init_set_str(n,argv[1],10);
    if(!mpz_probab_prime_p(n,36)){
        printf("It's prime....");
        return 0;
    }
    size_t s = mpz_sizeinbase(n,2);

    if(s<=24){
        shanks(p,q,n,INF_ITER);
    }else if(s<=128){
        if(!pollardrho_brent(p,q,n,INF_ITER)){
            if(!fermat(p,q,n,100000000)){
                shanks(p,q,n,INF_ITER);
            }
        }
    }else{
        printf("From 128bits using alternative approach is recommended");
        if(!pollardrho_brent(p,q,n,1000000000)){
            if(!fermat(p,q,n,100000000)){
                shanks(p,q,n,INF_ITER);
            }
        }
    }
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    //speedtest();
    //test_factorization();
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

    //mpz_set_str(N, "21", 10);
    mpz_set_str(N, "1790271774967384279", 10);
    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    puts("-----------Pollard Rho Method-----------");
    pollardrho(p, q, N, INF_ITER);
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
}


void speedtest(){
    mpz_t a,p,q;
    mpz_inits(p,q,a,NULL);
    mpz_set_str(a, "4E733FEBB94DB17CA3E6AA26EC33B4960C150C52300E06C60B3318F0744FEF2D687A8F5BF598894A22EEC4ABDAE01B197E4CC5603DE67EB670E261EB4E4CC5E26241EDCDE494CCE415BBC5A410ABCEFDFF6199BBCDF62E9D434FAA88A1D16012520F80D126208206FF80191E20ED7423CDCE5B8A555B4161534E789A74F0A701", 16);


    clock_t before;
    double  result;
    before  = clock();
    fermat(p,q,a,INF_ITER);
    result = (double)(clock() - before) / CLOCKS_PER_SEC;

    puts("-----------Fermat Method-----------");
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    printf("걸린시간은 %g 입니다.\n", result);


    mpz_set_str(a, "4126541", 10);
    before  = clock();
    shanks(p,q,a,INF_ITER);
    result = (double)(clock() - before) / CLOCKS_PER_SEC;

    puts("-----------Shanks Method-----------");
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    printf("걸린시간은 %g 입니다.\n", result);

    mpz_set_str(a, "9214228665501282294291731833", 10);
    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    before  = clock();
    pollardrho(p,q,a,INF_ITER);
    result = (double)(clock() - before) / CLOCKS_PER_SEC;

    puts("-----------pollard rho Method-----------");
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    printf("걸린시간은 %g 입니다.\n", result);

    mpz_set_ui(p,1);
    mpz_set_ui(q,1);
    before  = clock();
    pollardrho_brent(p,q,a,INF_ITER);
    result = (double)(clock() - before) / CLOCKS_PER_SEC;

    puts("-----------pollard rho brent Method-----------");
    gmp_printf("p = %Zd\n", p);
    gmp_printf("q = %Zd\n", q);
    printf("걸린시간은 %g 입니다.\n", result);

}