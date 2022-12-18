#include "stdio.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "float.h"

#define max(a,b) a>b?a:b
#define N 3
#define PI 3.141592654

//计算正太分布
long double gaussrand()
{
    static long double U, V;
    static int phase = 0;
    long double Z;
    if (phase == 0)
    {
        U = rand() / (RAND_MAX + 1.0);
        V = rand() / (RAND_MAX + 1.0);
        Z = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);
    }
    else {
        Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
    }
    phase = 1 - phase;
    return Z;
}

//交换
void swap_long_double(long double* a, long double* b) {
    long double tmp = *a;
    *a = *b;
    *b = tmp;
}

//交换
void swap_int(int* a, int* b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

//快速排序
void quick_sort(long double q[], int index[], int l, int r)
{
    if (l >= r) return;

    int i = l - 1, j = r + 1;
    long double k = q[(l + r) >> 1];
    while (i < j)
    {
        do i++; while (q[i] < k);
        do j--; while (q[j] > k);
        if (i < j) {
            swap_long_double(&q[i], &q[j]);
            swap_int(&index[i], &index[j]);
        }
    }
    quick_sort(q, index, l, j), quick_sort(q, index, j + 1, r);
}

//计算所有断点
void get_breakpoints(long double d_c[], const long double c[], const long double d[], int len) {
    //printf("breakpoints： \n");
    for (int i = 0; i < len; ++i) {
        d_c[i] = -(d[i] / c[i]);
        //printf("c%d = %Lf, d%d = %Lf, c%d/d%d = %Lf \n", i, c[i], i, d[i], i, i, d_c[i]);
    }
}

//计算函数值
long double f(const long double x, const long double a, const long double b, const long double c[], const long double d[], const int len) {
    long double sum = 0;
    for (int i = 0; i < len; ++i) {
        sum += max(0, c[i] * x + d[i]);
    }
    sum += (a / 2) * x * x + b * x;
    return sum;
}

//计算一阶导数的值
long double g(const long double x, const long double a, const long double b) {
    return a * x + b;
}

//使用快速排序寻找最小值
long double get_min_qsort(const long double a, const long double b, const long double c[], const long double d[], int len)
{
    printf("快速排序:\n");
    long double *d_c = (long double*)malloc(sizeof(long double) * len);
    int *index = (int*)malloc(sizeof(int) * len);
	//long double d_c[N];
    //int index[N];
    //设置索引数组
    for (int i = 0; i < len; ++i) {
        index[i] = i;
    }

    //计算初始截距
    long double intercept = b;
    for (int i = 0; i < len; ++i) {
        if (c[i] < 0)
            intercept += c[i];
    }

    //printf("初始截距: %Lf \n", intercept);

    //获取所有断点
    get_breakpoints(d_c, c, d, len);

    clock_t start, end;
    long double cpu_time_used;

    start = clock();

    //断点排序
    quick_sort(d_c, index, 0, len - 1);
    /*printf("排序后断点：\n");
    for (int i = 0; i < len; ++i) {
        printf("%Lf  ", d_c[i]);
    }
    printf("\n");
    printf("排序后断点索引：\n");
    for (int i = 0; i < len; ++i) {
        printf("%d  ", index[i]);
    }
    printf("\n");*/

    long double min_x = 0;

    //循环求解
    for (int i = 0; i < len; ++i) {
        long double gk_down = g(d_c[i], a, intercept);
        long double gk_up = g(d_c[i], a, intercept + fabs(c[index[i]]));

        //printf("gk_up: %Lf \n", gk_up);
        //printf("gk_down: %Lf \n", gk_down);


        //if (g(d_c[i], a, intercept) > 0) {
        //    //printf("intercept: %Lf \n", intercept);
        //    min_x = -intercept / a;
        //    break;
        //}
        //else if (g(d_c[i], a, intercept) == 0) {
        //    min_x = d_c[i];
        //    break;
        //}
        //else if (g(d_c[i], a, intercept + fabs(c[index[i]])) >= 0) {
        //    min_x = d_c[i];
        //    break;
        //}
        //else if (i == N - 1) {
        //    min_x = -(intercept + fabs(c[index[i]])) / a;
        //    break;
        //}

        if (gk_down > 0) {
            //printf("intercept: %Lf \n", intercept);
            min_x = -intercept / a;
            break;
        }
        else if (gk_down == 0) {
            min_x = d_c[i];
            break;
        }
        else if (gk_up >= 0) {
            min_x = d_c[i];
            break;
        }
        else if (i == len - 1) {
            min_x = -(intercept + fabs(c[index[i]])) / a;
            break;
        }

        intercept += fabs(c[index[i]]);
    }

    end = clock();
    cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
    printf("花费时间：%Lf s\n", cpu_time_used);
    printf("min_x: %Lf \n", min_x);
	
	free(d_c);
	free(index);

    //返回最小值
    return f(min_x, a, b, c, d, len);
	
}

//查找最小值
long double find_min(const long double a[], const int len) {
    long double min = DBL_MAX;
    for (int i = 0; i < len; ++i) {
        if (a[i] < min) min = a[i];
    }
    //printf("find_min_c: %Lf \n", min);
    return min;
}

//分割索引序列
int division_U(int U[], int l, int r, int k, const long double d_c[]) {
	int i = l, j = r - 1;
    long double v = d_c[U[k]];

    swap_int(&U[l], &U[k]);

    int x = U[l];

    while (i < j){
        while (i < j && d_c[U[j]] >= v) {
            j--;
        }
        if (i < j) {
            U[i] = U[j];
            i++;
        }
        while (i < j && d_c[U[i]] < v) {
            i++;
        }
        if (i < j) {
            U[j] = U[i];
            j--;
        }
    }

    U[i] = x;

    return i;
}

//计算区间内所有斜率的和
long double sum_fabs_c(const long double c[], const int U[], const int l, const int r, const long double d_c[]) {
    long double sum = 0;
    for (int i = l; i < r; ++i) {
    //printf("%Lf %Lf\n", d_c[U[i]], fabs(c[U[i]]));
        sum += fabs(c[U[i]]);
    }
    return sum;
}

//使用Efficient root finding寻找最小值
long double get_min_ERF(const long double a, const long double b, const long double c[], const long double d[], int len)
{
    printf("高效查找:\n");
	long double *d_c = (long double*)malloc(sizeof(long double) * len);
	int *U = (int*)malloc(sizeof(int) * len);
	//long double d_c[N];
    //int U[N];

    //设置索引数组
    for (int i = 0; i < len; ++i) {
        U[i] = i;
    }

    clock_t start, end;
    long double cpu_time_used;


    //获取所有断点
    get_breakpoints(d_c, c, d, len);

    int l = 0, r = len;

    /*for (int i = l; i < r; ++i) {
        printf("U[%d]: %d ", i, U[i]);
    }
    printf("\n");
    for (int i = l; i < r; ++i) {
        printf("xk[%d]: %Lf ", i, d_c[U[i]]);
    }
    printf("\n");*/

    start = clock();
    int k = 0;
    long double xk = find_min(d_c, len);
    long double gk = a * xk + b;
    for (int i = 0; i < len; ++i) {
        if ((c[i] * xk + d[i]) >= 0) {
            gk += c[i];
        }
    }
    //printf("初始gk = %Lf\n", gk);

    if (gk >= 0) {
        xk = xk - gk / a;
    }
    else {

        int cnt = 0;
        //int flag = 0;

        //求解x
        while (l < r) {

            int k_old = k;
            long double xk_old = xk;
            long double gk_old = gk;

            int p = l + rand() % (r - l);
            k = U[p];
            xk = d_c[k];


            //flag = 1;
            int tc = cnt++;
            //printf("\n\np%d: %d  k%d: %d  xk_old%d: %Lf  xk%d: %Lf\n", tc, p, tc, k, tc, xk_old, tc, xk);

            int dv = division_U(U, l, r, p, d_c);
            //k = U[dv];
            //printf("p%d: %d  k%d: %d  xk_old%d: %Lf  xk%d: %Lf\n", tc, p, tc, k, tc, xk_old, tc, xk);

            //printf("l%d: %d  dv%d: %d r%d: %d\n", tc, l, tc, dv, tc, r);

            if (xk >= xk_old) {
                gk = gk_old + a * (xk - xk_old) + sum_fabs_c(c, U, l, dv, d_c);
            }
            else {
                gk = gk - fabs(c[k_old]) - a * (xk_old - xk) - sum_fabs_c(c, U, dv, r, d_c);
            }

            //printf("gk_old%d: %Lf  gk%d: %Lf\n", tc, gk_old, tc, gk);

            if (gk < 0) {
                gk = gk + fabs(c[k]);
                if (gk >= 0) {
                    break;
                }
                else {
                    l = dv + 1;
                }
            }
            else {
                if (l < dv) {
                    gk = gk + fabs(c[k]);
                }
                r = dv;
            }

            if (l >= r) {
            //   printf("gk/a: %Lf\n", gk / a);
                xk = xk - gk / a;
                break;
            }
        }
    }
    end = clock();
    cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
    printf("花费时间：%Lf s\n", cpu_time_used);
    printf("min_x: %Lf \n", xk);
    
    free(d_c);
    free(U);
    //返回最小值
    return f(xk, a, b, c, d, len);
}


int main(int argc, char* argv[])
{

    long double a = 1.5;
    long double b = 1;
    long double c[N] = { 1, 1.2, -0.9 };
    long double d[N] = { 0.1, -1.4, -1.2 };

    //long double a = fabs(gaussrand());
    //long double b = gaussrand();
	
    //printf("a: %Lf b: %Lf\n", a, b);
    
    //int N = 0;
    //sscanf(argv[1], "%d", &N);
    
    printf("N = %d\n\n", N);
    
    //long double *c = (long double*)malloc(sizeof(long double) * N);
    //long double c[N];
    //for (int i = 0; i < N; ++i) {
    //    c[i] = gaussrand();
    //}
    
    //long double *d = (long double*)malloc(sizeof(long double) * N);
    //long double d[N];
    //for (int i = 0; i < N; ++i) {
    //    d[i] = gaussrand();
    //}

    printf("最小值为：%Lf \n\n", get_min_qsort(a, b, c, d, N));
    printf("最小值为：%Lf \n\n", get_min_ERF(a, b, c, d, N));
    
    //free(c);
    //free(d);
	
}
