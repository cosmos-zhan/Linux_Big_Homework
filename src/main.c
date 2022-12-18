#include "stdio.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "float.h"
#include "method.h"
#include "method_utils.h"
//#define N 3

int main(int argc, char* argv[])
{

    
    //long double a = 1.5;
    //long double b = 1;
    //long double c[N] = { 1, 1.2, -0.9 };
    //long double d[N] = { 0.1, -1.4, -1.2 };
    
    
    
    srand(time(NULL));
    long double a = fabs(gaussrand());
    long double b = gaussrand();
    while (a == 0 || !isnormal(a)) {
        a = fabs(gaussrand());
    }

    while (b == 0 || !isnormal(b)) {
        b = gaussrand();
    }
    
    //printf("a: %Lf b: %Lf\n", a, b);
    
    int N = 0;
    sscanf(argv[1], "%d", &N);
    printf("N = %d\n\n", N);
    
    long double *c = (long double*)malloc(sizeof(long double) * N);
    //long double c[N];
    for (int i = 0; i < N; ++i) {
        c[i] = gaussrand();
       while (c[i] == 0 || !isnormal(c[i])) {
            c[i] = gaussrand();
       }
    }
    long double *d = (long double*)malloc(sizeof(long double) * N);
    //long double d[N];
    for (int i = 0; i < N; ++i) {
        d[i] = gaussrand();
        while (d[i] == 0 || !isnormal(d[i])) {
            d[i] = gaussrand();
       }
    }
    

    printf("最小值为：%Lf \n\n", get_min_qsort(a, b, c, d, N));
    printf("最小值为：%Lf \n\n", get_min_ERF(a, b, c, d, N));
    printf("最小值为：%Lf \n\n", get_min_Tree(a, b, c, d, N));
    printf("最小值为：%Lf \n\n", get_min_genetic(a, b, c, d, N));
    printf("最小值为：%Lf \n\n", get_min_SA(a, b, c, d, N));
    
    free(c);
    free(d);
	
}

