#include "LKH.h"

/*
 * The GenerateCandidates function associates to each node a set of incident 
 * candidate edges. The candidate edges of each node are sorted in increasing
 * order of their Alpha-values.
 *
 * The parameter MaxCandidates specifies the maximum number of candidate edges 
 * allowed for each node, and MaxAlpha puts an upper limit on their 
 * Alpha-values.
 *
 * A non-zero value of Symmetric specifies that the candidate set is to be
 * complemented such that every candidate edge is associated with both its 
 * two end nodes (in this way MaxCandidates may be exceeded). 
 *
 * The candidate edges of each node is kept in an array (CandidatSet) of
 * structures. Each structure (Candidate) holds the following information:
 *
 *      Node *To    : pointer to the other end node of the edge
 *      int Cost    : the cost (length) of the edge
 *      int Alpha   : the alpha-value of the edge
 *
 * The algorithm for computing Alpha-values in time O(n^2) and space O(n) 
 * follows the description in
 *
 *      Keld Helsgaun,
 *      An Effective Implementation of the Lin-Kernighan Traveling 
 *      Salesman Heuristic,
 *      Report, RUC, 1998. 
 */

static int Max(const int a, const int b)
{
    return a > b ? a : b;
}
// 注意，这里的MaxCandidates为Node数-1（AscentCandidates)，即为每个节点生成除本身之外的所有candidate结构
void GenerateCandidates(int MaxCandidates, GainType MaxAlpha,
                        int Symmetric)
{
    Node *From, *To;
    Candidate *NFrom, *NN;
    int a, d, Count;

    if (TraceLevel >= 2)
        printff("Generating candidates ... ");
    if (MaxAlpha < 0 || MaxAlpha > INT_MAX)
        MaxAlpha = INT_MAX;
    /* Initialize CandidateSet for each node */
    FreeCandidateSets();
    From = FirstNode;
    do
        From->Mark = 0; // mark用于标记是否已经访问
    while ((From = From->Suc) != FirstNode);

    if (MaxCandidates > 0) {
        do { // 为每个Node的CandidateSet指针分配空间
            assert(From->CandidateSet =
                   (Candidate *) malloc((MaxCandidates + 1) *
                                        sizeof(Candidate)));
            From->CandidateSet[0].To = 0;
        }
        while ((From = From->Suc) != FirstNode);
    } else { // 否则的话通过读入文件将候选集写入
        AddTourCandidates();
        do {
            if (!From->CandidateSet)
                eprintf("MAX_CANDIDATES = 0: No candidates");
        } while ((From = From->Suc) != FirstNode);
        return;
    }

    /* Loop for each node, From */
    do { // 这个do while 相当于论文中的外层for循环，用于以From指针从FirstNode按拓补顺序开始遍历最小生成树
        NFrom = From->CandidateSet;
        if (From != FirstNode) { // 不是FirstNode 则从该点一直向Dad节点回溯，直到根节点（Dad = 0)
            // 这里相当于论文中的第一个内层for循环，用于从From节点通过Dad域向根节点回溯来求（From，To）的beta值
            // 不同的是，这里需要判断条件即论文中的i!=1
            From->Beta = INT_MIN;
            for (To = From; To->Dad != 0; To = To->Dad) {
                To->Dad->Beta =
                    !FixedOrCommon(To, To->Dad) ? // 若To的父节点不包含在To的固定边中，且没有与给定的MergeTourFile相同的边，则按照max计算，若相同则直接就是beta值
                    Max(To->Beta, To->Cost) : To->Beta;
                To->Dad->Mark = From;
            }
        }
        Count = 0; // Count 用于统计当前已经插入了几个Candidate
        /* Loop for each node, To */
        To = FirstNode;
        // 这段do while循环相当于论文中的第二个内层for循环，用于从To节点向树的底部前向求第一个循环中未被mark标记的（From，To）的beta值
        // 不同的是，这里需要处理论文中说到的两个节点相等的问题(From == FirstNode),以及任何一个节点等于1的问题（紧接着的两个else)
        // 当From和To的节点中有一个等于1，则表示论文中求alpha值的第二种情况，需要用d - N1节点新增的那条成环的最长边
        do {
            if (To == From)
                continue; // 对应论文中的 i 不等于 j
            d = c && !FixedOrCommon(From, To) ? c(From, To) : D(From, To);
            if (From == FirstNode) // 特殊情况一：外部节点指针From为FirstNode
                a = To == From->Dad ? 0 : d - From->NextCost; // 当To刚好为FirstNode的Dad时，表示该节点为1-tree中的一条边，所以a为0，否则的话为d - nextCost
            else if (To == FirstNode)
                a = From == To->Dad ? 0 : d - To->NextCost;
            else {
                if (To->Mark != From)
                    To->Beta =
                        !FixedOrCommon(To, To->Dad) ?
                        Max(To->Dad->Beta, To->Cost) : To->Dad->Beta;
                a = d - To->Beta;
            }
            if (FixedOrCommon(From, To))
                a = INT_MIN;
            else {
                if (From->FixedTo2 || To->FixedTo2 || Forbidden(From, To))
                    continue;
                if (InInputTour(From, To)) {
                    a = 0;
                    if (c)
                        d = D(From, To);
                } else if (c) {
                    if (a > MaxAlpha ||
                        (Count == MaxCandidates &&
                         (a > (NFrom - 1)->Alpha ||
                          (a == (NFrom - 1)->Alpha
                           && d >= (NFrom - 1)->Cost))))
                        continue;
                    if (To == From->Dad) {
                        d = From->Cost;
                        a = 0;
                    } else if (From == To->Dad) {
                        d = To->Cost;
                        a = 0;
                    } else {
                        a -= d;
                        a += (d = D(From, To));
                    }
                }
            }
            if (a <= MaxAlpha && IsPossibleCandidate(From, To)) {
                /* Insert new candidate edge in From->CandidateSet */
                NN = NFrom;
                while (--NN >= From->CandidateSet) {
                    if (a > NN->Alpha || (a == NN->Alpha && d >= NN->Cost))
                        break;
                    *(NN + 1) = *NN;
                }
                NN++;
                NN->To = To;
                NN->Cost = d;
                NN->Alpha = a;
                if (Count < MaxCandidates) {
                    Count++;
                    NFrom++;
                }
                NFrom->To = 0;
            }
        }
        while ((To = To->Suc) != FirstNode);
    }
    while ((From = From->Suc) != FirstNode);

    AddTourCandidates();
    if (Symmetric)
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2)
        printff("done\n");
}
