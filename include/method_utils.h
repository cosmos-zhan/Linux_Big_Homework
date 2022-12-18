#ifndef METHOD_UTILS_H_
#define METHOD_UTILS_H_

//计算正太分布
long double gaussrand();

//交换长浮点型
void swap_long_double(long double* a, long double* b);

//交换整形
void swap_int(int* a, int* b);

//快速排序
void quick_sort(long double q[], int index[], int l, int r);

//计算所有断点
void get_breakpoints(long double d_c[], const long double c[], const long double d[], int len);

//计算函数值
long double f(const long double x, const long double a, const long double b, const long double c[], const long double d[], const int len);

//计算一阶导数的值
long double g(const long double x, const long double a, const long double b);

//查找断点最小值
long double find_min(const long double a[], const int len);

//分割索引序列
int division_U(int U[], int l, int r, int k, const long double d_c[]);

//计算区间内所有斜率的和
long double sum_fabs_c(const long double c[], const int U[], const int l, const int r, const long double d_c[]);

typedef struct TreeNodes
{
    double data;
    int index;
    struct TreeNodes* father;
    struct TreeNodes* left;
    struct TreeNodes* right;
}TreeNode;

//遍历二叉树
void Sort_Tree(TreeNode* bt, double key,int dex) ;
void preorder_add(TreeNode* root,long double* sum,long double c[]);
void preorder_sub(TreeNode* root,long double* sum,long double c[]);

//查找一个中序遍历二叉树节点的前驱节点
TreeNode* get_pre(TreeNode* root);

//查找一个中序遍历二叉树节点的后继节点
TreeNode* get_suc(TreeNode* root);

#endif /* METHOD_UTILS_H_ */
