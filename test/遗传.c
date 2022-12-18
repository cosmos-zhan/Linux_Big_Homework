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
        numx[i] = long double(rand()%10000)/1000 ;//生成20个100以内的随机数
    }

    clock_t start, end;
    long double cpu_time_used;
    start = clock();  //开始计时

	for (int i = 0 ; i < 5; i++)//迭代5次
	{
        
        long double count[20];//20个x结果
        for (int j = 0; j < 20; j++)
        {
            count[j] = a * 0.5 * numx[j] * numx[j] + b * numx[j];//外部
            for (int k = 0; k < m; k++)//内部循环
            {
                long double num = c[k] * long double(numx[j]) + d[k];
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
                numer[q][p] = int(tep) % 2;
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
            int temp = 0;
            for (int p = 0; p < 16; p++)
            {
                temp += numer[q][p] * pow(2 , p);
            }
            numx[q] = float(temp) / 1000;
        }
        
	}

    end = clock();
    cpu_time_used = ((long double)(end - start)) / CLOCKS_PER_SEC;
    printf("花费时间：%Lf s\n", cpu_time_used);
    printf("min_x: %Lf \n", numx[0]);
    long double x_min = numx[0];
    return f(x_min, a, b, c, d, len);
}
