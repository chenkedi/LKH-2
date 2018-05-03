#include "LKH.h"
#include "Heap.h"

/*
 * dad为最小生成树结构的索引，该树以id等于1的点为root
 * The MinimumSpanningTree function determines a minimum spanning tree using 
 * Prim's algorithm.
 *
 * At return the Dad field of each node contains the father of the node, and 
 * the Cost field contains cost of the corresponding edge. The nodes are 
 * placed in a topological ordered list, i.e., for any node its father precedes 
 * the node in the list. The fields Pred and Suc of a node are pointers to the 
 * predecessor and successor node in this list.
 *
 * The function can be used to determine a minimum spanning tree in a dense 
 * graph, or in a sparse graph (a graph determined by a candidate set).
 *
 * When the graph is sparse a priority queue, implemented as a binary heap, 
 * is used  to speed up the determination of which edge to include next into 
 * the tree. The Rank field of a node is used to contain its priority (usually 
 * equal to the shortest distance (Cost) to nodes of the tree).        
 */

void MinimumSpanningTree(int Sparse)
{
    Node *Blue;         /* Points to the last node included in the tree */
    Node *NextBlue = 0; /* Points to the provisional next node to be included */
    Node *N;
    Candidate *NBlue;
    int d;

    Blue = N = FirstNode;
    Blue->Dad = 0;              /* The root of the tree has no father */
    if (Sparse && Blue->CandidateSet) {
        /* The graph is sparse */
        /* Insert all nodes in the heap */
        Blue->Loc = 0;          /* A blue node is not in the heap */
        while ((N = N->Suc) != FirstNode) {
            N->Dad = Blue;
            N->Cost = N->Rank = INT_MAX;
            HeapLazyInsert(N);
        }
        /* Update all neighbors to the blue node */
        for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) { // 遍历blue的candidateSet
            if (FixedOrCommon(Blue, N)) {
                N->Dad = Blue;
                N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                N->Rank = INT_MIN;
                HeapSiftUp(N);
            } else if (!Blue->FixedTo2 && !N->FixedTo2) {
                N->Dad = Blue;
                N->Cost = N->Rank = NBlue->Cost + Blue->Pi + N->Pi;
                HeapSiftUp(N);
            }
        }
        /* Loop as long as there are more nodes to include in the tree */
        while ((NextBlue = HeapDeleteMin())) {
            Follow(NextBlue, Blue);
            Blue = NextBlue;
            /* Update all neighbors to the blue node */
            for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
                if (!N->Loc)
                    continue;
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                    N->Rank = INT_MIN;
                    HeapSiftUp(N);
                } else if (!Blue->FixedTo2 && !N->FixedTo2 &&
                           (d = NBlue->Cost + Blue->Pi + N->Pi) < N->Cost) {
                    N->Dad = Blue;
                    N->Cost = N->Rank = d;
                    HeapSiftUp(N);
                }
            }
        }
    } else { // 这种实现的时间复杂度应该是N平方，应该存在改进空间
        /* The graph is dense */
        while ((N = N->Suc) != FirstNode)
            N->Cost = INT_MAX;
        /* Loop as long as there a more nodes to include in the tree */
        while ((N = Blue->Suc) != FirstNode) {
            int Min = INT_MAX;
            /* Update all non-blue nodes (the successors of Blue in the list) */
            do {
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = D(Blue, N);
                    NextBlue = N;
                    Min = INT_MIN;
                } else {
                    if (!Blue->FixedTo2 && !N->FixedTo2 &&
                        !Forbidden(Blue, N) &&
                        (!c || c(Blue, N) < N->Cost) &&
                        (d = D(Blue, N)) < N->Cost) {
                        N->Cost = d;
                        N->Dad = Blue;
                    }
                    if (N->Cost < Min) {
                        Min = N->Cost;
                        NextBlue = N;
                    }
                }
            }
            while ((N = N->Suc) != FirstNode);
            Follow(NextBlue, Blue); // nextBlue的前一个节点与后一个节点进行双向互联，nextBlue自己与自己互连，nextBlue与Blue的前一个节点互联，这两个节点互联
            Blue = NextBlue;
        }
    }
}
