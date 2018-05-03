#include "LKH.h"

/*
 * The Minimum1TreeCost function returns the cost of a minimum 1-tree.
 *
 * The minimum 1-tre is found by determining the minimum spanning tree and 
 * then adding an edge corresponding to the second nearest neighbor of one 
 * of the leaves of the tree (any node which has degree 1). The leaf chosen
 * is the one that has the longest second nearest neighbor distance.
 * 这里为何要选择所有叶子节点中具有最长的最近邻居的点为N1,是个疑问
 * The V-value of a node is its degree minus 2. Therefore, Norm being the 
 * sum of squares of all V-values, is a measure of a minimum 1-tree/s 
 * discrepancy from a tour. If Norm is zero, then the 1-tree constitutes a 
 * tour, and an optimal tour has been found.
 */

GainType Minimum1TreeCost(int Sparse)
{
    Node *N, *N1 = 0;
    GainType Sum = 0;
    int Max = INT_MIN;

    MinimumSpanningTree(Sparse);
    N = FirstNode;
    do {
        N->V = -2;
        Sum += N->Pi;
    }
    while ((N = N->Suc) != FirstNode);
    Sum *= -2; // 这里为何要乘以-2？是因为下界为2pi + cij吗
    // 由于已经生成了最小生成树，则该树中的节点的度通过下面的循环计算
    // 首先每个进入最小生成树的节点至少度为1，所以有v++
    // 剩下的度由其指向该节点的子节点的Dad字段决定，所以有Dad->V++
    // 循环完成后即可得到最小生成树中所有节点的正确V值
    while ((N = N->Suc) != FirstNode) {
        N->V++;
        N->Dad->V++;
        Sum += N->Cost;
        N->Next = 0;
    }
    // 为何要将生成树的根节点（FirstNode）的Dad指向其直接后继？
    FirstNode->Dad = FirstNode->Suc;
    FirstNode->Cost = FirstNode->Suc->Cost;
    // 生成树中，v为-1的点是叶子节点, N1为最终选出的连成环的叶子节点，作为mst向1-tree转换的特殊节点
    do {
        if (N->V == -1) {
            // Connect 求度为1的某个节点到其他节点且形成的边不在MST中的最短边
            Connect(N, Max, Sparse);
            // 只有当前节点N的最近边大于上一个节点的最近边（被Max记录）时，N1节点才会被替换，即头部注释中的longest
            // 最大的问题在于，Connect函数求出的是每个度为1的节点最近的邻居，而不是第二近的邻居 ？？？？
            if (N->NextCost > Max && N->Next) {
                N1 = N;
                Max = N->NextCost;
            }
        }
    }
    while ((N = N->Suc) != FirstNode);
    assert(N1);
    // 特殊节点选择成功后，两个节点的degree都加一
    N1->Next->V++;
    N1->V++;
    // Sum 为所有MST中的边的cost再加上N1节点新增成环边的cost
    Sum += N1->NextCost;
    // Norm 为每个节点度减2的平方和，即v的平方和
    Norm = 0;
    do
        Norm += N->V * N->V;
    while ((N = N->Suc) != FirstNode);
    // 这里的if else的作用应该是保证特殊节点一定为FirstNode
    if (N1 == FirstNode) // 如果N1节点就是First节点，则将其后继节点的Dad设置为0，即将其后继节点作为1-tree的Root节点
        N1->Suc->Dad = 0;
    else { // 否则的话仍将FirstNode作为1-tree的根，并将N1移动到FirstNode的直接前继节点，并替换FirstNode
        FirstNode->Dad = 0;
        Precede(N1, FirstNode);
        FirstNode = N1;
    }
    if (Norm == 0) {
        for (N = FirstNode->Dad; N; N1 = N, N = N->Dad)
            Follow(N, N1);
        for (N = FirstNode->Suc; N != FirstNode; N = N->Suc) {
            N->Dad = N->Pred;
            N->Cost = D(N, N->Dad);
        }
        FirstNode->Suc->Dad = 0;
    }
    return Sum;
}
