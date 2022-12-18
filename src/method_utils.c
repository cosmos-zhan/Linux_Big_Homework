#include "method_utils.h"
#include "math.h"
#include "stdio.h"
#include "float.h"
#include "stdlib.h"
#include "string.h"

#define max(a,b) a>b?a:b
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

//交换长浮点型
void swap_long_double(long double* a, long double* b) {
    long double tmp = *a;
    *a = *b;
    *b = tmp;
}


//交换整形
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


//查找断点最小值
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

//排序二叉树
void Sort_Tree(TreeNode* bt, double key,int dex)  //在二叉排序树中插入查找关键字key；
{
    TreeNode* parent;
    TreeNode* p = (TreeNode*)malloc(sizeof(TreeNode));
    p->data = key;   
    p->index=dex;            //保存结点数据；
    p->left = p->right = NULL;  //左右子树置空；
    TreeNode* head = bt;
    while (head)                //查找关键字所在的位置；
    {
        parent = head;
        if (key < head->data)     //如果关键字小于结点的数据；
            head = head->left; //在左子树上查找；
        else                   //若关键字大于结点的数据
            head = head->right; //再右子树上查找；
    }
    //判断添加到左子树还是右子树；
    if (key < parent->data){  //小于父结点；
        parent->left = p; //添加到左子树；
        p->father=parent;
    }
    else{                 //大于父结点；
        parent->right = p; //添加到右子树；
        p->father=parent;
    }   
}

//遍历二叉树
void mid_Order(TreeNode* T){
    if(T != NULL){
        mid_Order(T->left);//先访问左结点
        printf("i=%d,data=%f     ",T->index,T->data);
        mid_Order(T->right);//最后访问右结点
    }
    return;
}

void preorder_add(TreeNode* root,long double* sum,long double c[]) {
    if (root == NULL) {
        return;
    }
    *sum += fabs(c[root->index]);
    preorder_add(root->left, sum,c);
    preorder_add(root->right, sum,c);
}
void preorder_sub(TreeNode* root,long double* sum,long double c[]) {
    if (root == NULL) {
        return;
    }
    *sum -= fabs(c[root->index]);
    preorder_sub(root->left, sum,c);
    preorder_sub(root->right, sum,c);
}

//查找二叉树的中序前驱、后继节点
TreeNode* get_pre(TreeNode* root)/*查找一个中序遍历二叉树节点的前驱节点*/
{
    if (root == NULL)
    {
        return NULL;
    }
    if (root->left != NULL)
    {
        TreeNode* start = root->left;
        while (start->right != NULL)
        {
            start = start->right;
        }
        return start;
    }
    else if (root->father != NULL)
    {
        TreeNode* fa = root->father;
        TreeNode* current = root;
        while (fa && fa->left == current)
        {
            current = fa;
            fa = fa->father;
        }
        if (fa != NULL)
        {
            return fa;
        }
        return NULL;
    }
    return NULL;
}

TreeNode* get_suc(TreeNode* root)/*查找一个中序遍历二叉树节点的后继节点*/
{
    if (root == NULL)
    {
        return NULL;
    }
    if (root->right != NULL)
    {
        TreeNode* rig = root->right;
        while (rig->left != NULL)
        {
            rig = rig->left;
        }
        return rig;
    }
    else if (root->father != NULL)
    {
        TreeNode* fa = root->father;
        TreeNode* current = root;
        while (fa && fa->right == current)
        {
            current = fa;
            fa = fa->father;
        }
        if (fa != NULL)
        {
            return fa;
        }
        return NULL;
    }
    return NULL;
}

