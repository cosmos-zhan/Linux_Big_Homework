#ifndef METHOD_H_
#define METHOD_H_

//快速排序求解最小值
long double get_min_qsort(const long double a, const long double b, const long double c[], const long double d[], int len);

//论文中的高效方法求解最小值
long double get_min_ERF(const long double a, const long double b, const long double c[], const long double d[], int len);

//排序二叉树方法求最小值
long double get_min_Tree(const long double a, const long double b, long double c[], const long double d[], int len);

//模拟退火方法求解最小值
long double get_min_SA(const long double a, const long double b, const long double c[], const long double d[], int m);

//遗传算法求解最小值
long double get_min_genetic(const long double a, const long double b, const long double c[], const long double d[], int len);

#endif /* METHOD_H_ */
