#include "method.h"
#include "method_utils.h"
#include "stdio.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "time.h"
#include "float.h"

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

        //int cnt = 0;
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
            //int tc = cnt++;
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

//利用二叉树进行求解
long double get_min_Tree(const long double a, const long double b, long double c[], const long double d[], int len)
{
    printf("二叉树:\n");
    long double* d_c = (long double*)malloc(sizeof(long double) * len);

    //计算初始截距
    long double intercept = b;
    for (int i = 0; i < len; ++i) {
        if (c[i] < 0)
            intercept += c[i];
    }

    //获取所有断点
    get_breakpoints(d_c, c, d, len);
    
    clock_t start, end;
    long double cpu_time_used;

    start = clock();

    // TreeNode* root = createTree(d_c, len);
    //初始化根结点；
	TreeNode* root=(TreeNode*)malloc(sizeof(TreeNode));
	root->data=d_c[0];
    root->index=0;
	root->left=root->right=NULL;
	//将其他数组数据进行排序；
	for(int i=1;i<len;i++){
		Sort_Tree(root,d_c[i],i);
    }

    //mid_Order(root);
    preorder_add(root->left, &intercept,c);
    long double x_min=0;
    long double gk;
    while (root) {
        gk = g(root->data, a, intercept);
        if (gk == 0) {
            x_min = root->data;
            break;
        }
        if (gk < 0) {
            if (g(root->data, a, (intercept + fabs(c[root->index]))) >= 0) {
                x_min = root->data;
                break;
            }
            TreeNode* last = (TreeNode*)malloc(sizeof(TreeNode));
            last=get_suc(root);

            if (last != NULL) {
                int n=root->index;
                if (g(last->data, a, (intercept + fabs(c[n]))) >= 0) {
                    x_min = -(intercept + fabs(c[n])) / a;
                    break;
                }
            }
            if (root->right != NULL) {
                intercept+=fabs(c[root->index]);
                root = root->right;
                preorder_add(root->left, &intercept,c);
            }
            else {
                break;
            }
        }
        if (gk > 0) {
            if (root->left != NULL) {
                TreeNode* pre = NULL;
                pre=get_pre(root);
                if (pre != NULL) {
                    if (g(pre->data, a, intercept) <= 0) {
                        x_min = -intercept / a;
                        break;
                    }
                }
                root = root->left;
                preorder_sub(root->right, &intercept,c);
                intercept -= fabs(c[root->index]);
            }
            else {
                break;
            }

        }
    }
    end = clock();
    cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
    printf("花费时间：%Lf s\n", cpu_time_used);
    printf("min_x: %Lf \n", x_min);
    //返回最小值
    return f(x_min, a, b, c, d, len);
}

//使用Simulated Annealing寻找最小值
long double get_min_SA(const long double a, const long double b, const long double c[], const long double d[], int m) {
	int L = 100;
	double alpha = 0.99;
	int T0 = 100;
	int Tmin = 0.01;
	int T = T0;
	long double result = DBL_MAX;
	long double x[100] = { 0 };
	for (int i = 0; i <= 99; i++) {
		x[i] = (rand() % 10000 / (float)10000 * 2 - 1);
	}
	printf("模拟退火算法查找:\n");
	clock_t start, end;
    	long double cpu_time_used;
	start = clock();
	
	while (T > Tmin) {
		for (int i = 0; i < L; i++) {
			long double funTmp = f(x[i], a, b, c, d, m);
			long double x_new = x[i] + (rand() % 10000 / (float)10000 * 2 - 1) * T;
			if (1) {
				long double funTmp_new = f(x_new, a, b, c, d, m);
				if (funTmp_new - funTmp < 0) {
					x[i] = x_new;
				}
				else {
					long double p = exp(-(funTmp_new - funTmp) / T);
					if ((rand() % 10000 / (float)10000) < p) {
						x[i] = x_new;
					}
				}
			}
		}
		T = T * alpha;
	}
	for (int i = 0; i < L; i++) {
		long double temp = f(x[i], a, b, c, d, m);
		if (result > temp)
			result = temp;
		//printf("result=%.10f ", result);
	}
	end = clock();
    	cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
	//printf("%f\n", result);
	printf("花费时间：%Lf s\n", cpu_time_used);
    	printf("min_x: %Lf \n", x[99]);
    	
	return result;
}

//使用遗传算法得到的最小值
long double get_min_genetic(const long double a, const long double b, const long double c[], const long double d[], int len)
{
//    long double a = gaussrand();
//    if (a < 0) { a = -a; }
//    long double b = gaussrand();

//	int m = 10;//m的取值
//    long double numa[200000];//c,d取值
//    for (int i = 0; i < 2*m; i++ )
//    {
//        numa[i] = gaussrand();
//    }
    printf("遗传算法\n");
    int m=len;
    long double numx[20];//随机生成20个x作为初始种群
    for (int i = 0; i < 20; i++)
    {
        //srand(time(0));
        numx[i] = (long double)(rand()%10000)/1000 ;//生成20个100以内的随机数
    }

    clock_t start, end;
    long double cpu_time_used;
    start = clock();  //开始计时

	for (int i = 0 ; i < 1000; i++)//迭代1000次
	{
        
        long double count[20];//20个x结果
        for (int j = 0; j < 20; j++)
        {
            count[j] = a * 0.5 * numx[j] * numx[j] + b * numx[j];//外部
            for (int k = 0; k < m; k++)//内部循环
            {
                long double num = c[k] * (long double)(numx[j]) + d[k];
                if (num > 0) { count[j] += num; }
            }//cout  << numx[j] << " " << count[j] << endl;
        }
        
        int numer[20][16];//20个16位2进制
        float pp[3];
        long double b1=count[0];
        pp[0] = numx[0]; pp[1] = numx[1]; pp[2] = numx[2];
        for( int q=0;q<20;q++)//最优的3个解保留
        {
            if (count[q] < b1)
            {
                b1 = count[q];
                pp[2] = pp[1];
                pp[1] = pp[0];
                pp[0] = numx[q];
             }//cout << pp[0] << endl;
        }
        
        //得到3个最小值的解的二进制
        for (int q = 0; q < 3; q++)
        {
            int tep = pp[q] * 1000;
            for (int p = 0; p < 16; p++)
            {
                numer[q][p] = (int)(tep) % 2;
               // cout <<  "  " << numer[q][p] ;
                tep = tep / 2;
               
            }
            //cout << endl;
        }

        //交叉变异
        float pc = 0.9;//0.9的概率从父辈遗传，0.1的概率发生随机变异
        for (int q = 3; q < 20;q++)
        { //cout << endl;
            for (int p = 0; p < 16; p++)
            {
               
                float tep = rand() / (RAND_MAX + 1.0);
                if (pc > tep)//模拟遗传
                {
                    if (tep < 0.3) { numer[q][p] = numer[0][p]; }
                    else if (tep < 0.6) { numer[q][p] = numer[1][p]; }
                    else { numer[q][p] = numer[2][p]; }
                }
                else//模拟变异
                {
                    numer[q][p] = rand() % 2;
                }
                //cout << "  " << numer[q][p];
            }//cout << endl;
        }

        for (int q = 0; q < 20; q++)
        {
            long long temp = 0;
            for (int p = 0; p < 16; p++)
            {
                temp += numer[q][p] * pow(2 , p);
            }
            numx[q] = (long double)(temp) / 1000;
        }
        
	}

    end = clock();
    cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
    printf("花费时间：%Lf s\n", cpu_time_used);
    printf("min_x: %Lf \n", numx[0]);
    long double x_min = numx[0];
    return f(x_min, a, b, c, d, len);
}

