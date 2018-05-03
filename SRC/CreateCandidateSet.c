#include "LKH.h"

/*
 * 该函数主要做以下几点：
 * 1、是否有给定的PiFile和CandidateFile，如果又则读取，若存在PiFile，则不需要进行Ascend过程，可以认为次梯度优化已经结束，直接可以求cost
 *    若存在CandidateSetFile，则需要经CandidateSet加入当当前执行的上下文中
 * 2、根据CandidateSetType的特殊类型，调用对应特殊的createCandidate函数
 * 3、没有PiFile时，执行Ascent进行次梯度优化
 * The CreateCandidateSet function determines for each node its set of incident
 * candidate edges.
 * 当点权pi值在文件中没有给出时，该函数才会使用次梯度优化迭代计算lower bound；当给出了该
 * 文件，则直接读取，并在构建最小1-tree后直接计算出lower bound
 * The Ascent function is called to determine a lower bound on the optimal tour 
 * using subgradient optimization. But only if the penalties (the Pi-values) is
 * not available on file. In the latter case, the penalties is read from the 
 * file, and the lower bound is computed from a minimum 1-tree.      
 *
 * The function GenerateCandidates is called to compute the Alpha-values and to 
 * associate to each node a set of incident candidate edges.  
 *
 * The CreateCandidateSet function itself is called from LKHmain.
 *
 * 在构建候选集时，需要将DELAUNAY这种candidateType单独考虑
 */

void CreateCandidateSet()
{
    GainType Cost, MaxAlpha, A;
    Node *Na;
    int CandidatesRead = 0, i;
    double EntryTime = GetTime();

    Norm = 9999;
    // C 表示距离计算函数指针，可以指向多种不同计算函数
    if (C == C_EXPLICIT) {// 当距离有距离矩阵表示时

        Na = FirstNode;
        // 从首节点一直往下遍历，为每个节点的C(距离）乘以precision，以放大到整数级别进行计算
        // 由于距离计算出来后是以下三角矩阵拉直成一个数组储存，每个节点的C都指向数组中，该节点这一行起始位置的前一个位置（因为id从1开始，而数组从0开始）
        // 所以通过如下方式可以为每个距离乘以精度
        do {
            for (i = 1; i < Na->Id; i++)
                Na->C[i] *= Precision;
        }
        while ((Na = Na->Suc) != FirstNode);
    }
    /**
     * 第一、二个大if用于读取可能存在的给定的Pi文件和初始解和初始candidate文件
     */
    // 当距离计算函数为直接返回1，或者maxTrials为0且问题已经给定初始tour（等同于读给定了初始候选集）
    // 或者初始化算法为给定的两种时，直接从文件中读取
    if (Distance == Distance_1 ||
        (MaxTrials == 0 &&
         (FirstNode->InitialSuc || InitialTourAlgorithm == SIERPINSKI ||
          InitialTourAlgorithm == MOORE))) {
        // 从文件中读取候选集
        CandidatesRead = ReadCandidates(MaxCandidates);
        // AddTourCandidates方法用于向候选集中添加problem给定的需要存在的路线
        AddTourCandidates();
        if (ProblemType == HCP || ProblemType == HPP)// 如果问题为汉密尔顿环，则开始次梯度优化
            Ascent();
        goto End_CreateCandidateSet;
    }


    if (TraceLevel >= 2)
        printff("Creating candidates ...\n");
    if (MaxCandidates > 0 && // candidateSetType表示使用何种度量来生成邻居，如alpha表示使用Alpha度量
        (CandidateSetType == QUADRANT || CandidateSetType == NN)) {
        ReadPenalties(); // 如果maxCandidate大于0，且为给定的candidateSetType中的一种，若已经给定每个点的pi值，则读取
        // 如果ReadCandidates方法没有读取到候选集文件，且MaxCandidates要求大于0，则根据候选集
        // 度量类型调用对应的createCandidate函数重新构建候选集
        if (!(CandidatesRead = ReadCandidates(MaxCandidates)) &&
            MaxCandidates > 0) {
            if (CandidateSetType == QUADRANT)
                CreateQuadrantCandidateSet(MaxCandidates);
            else if (CandidateSetType == NN)
                CreateNearestNeighborCandidateSet(MaxCandidates);
        }
        // 然后添加问题中指定的路线
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        goto End_CreateCandidateSet;
    }

    /**
     * Minimum1TreeCost会在Ascent中被调用，用于计算一个最小生成树的cost
     */

    if (!ReadPenalties()) {
        /* No PiFile specified or available */
        Na = FirstNode; // 若没有pi File，则初始化各个节点pi值
        do
            Na->Pi = 0;
        while ((Na = Na->Suc) != FirstNode);
        CandidatesRead = ReadCandidates(MaxCandidates); // 若有candidate文件则读取
        Cost = Ascent(); // 进行次梯度优化，求得cost下界，并得到alpha度量矩阵，生成候选集合
        // Subgradient 表示每个节点的pi值是否应该由次梯度优化决定
        if (Subgradient && SubproblemSize == 0) { // 只有当子问题中的节点数为0时才写pi文件
            WritePenalties();
            PiFile = 0;
        }
        // 存在PiFile，并且（MaxCandidate为0，或者已有candidate文件）则表示次梯度优化已经完成，可以结束createCandidate函数
    } else if ((CandidatesRead = ReadCandidates(MaxCandidates)) ||
               MaxCandidates == 0) {
        AddTourCandidates();
        if (CandidateSetSymmetric)
            SymmetrizeCandidateSet();
        goto End_CreateCandidateSet;
    } else {
        // 无已有的Candidate file，但存在已有的PiFile，并且maxCandidates不为0,
        // 则直接格局Pi值和Minimum1TreeCost来生成cost，而不用进行Ascent
        if (CandidateSetType != DELAUNAY && MaxCandidates > 0) {
            if (TraceLevel >= 2)
                printff("Computing lower bound ... ");
            Cost = Minimum1TreeCost(0);
            if (TraceLevel >= 2)
                printff("done\n");
        } else {
            // 只有当候选集类型为DELAUNAY类型时，才需要调用对应的createCandidate函数后再求cost
            CreateDelaunayCandidateSet();
            Na = FirstNode;
            do { // 这里为何要将每个节点的pi值赋值为0？
                Na->BestPi = Na->Pi;
                Na->Pi = 0;
            }
            while ((Na = Na->Suc) != FirstNode);
            if (TraceLevel >= 2)
                printff("Computing lower bound ... ");
            Cost = Minimum1TreeCost(1);
            if (TraceLevel >= 2)
                printff("done\n");
            Na = FirstNode;
            do {
                Na->Pi = Na->BestPi;
                Cost -= 2 * Na->Pi;
            }
            while ((Na = Na->Suc) != FirstNode);
        }
    }


    LowerBound = (double) Cost / Precision;
    if (TraceLevel >= 1) { // 根据已知的optimum求与最优解之间的gap
        printff("Lower bound = %0.1f", LowerBound);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.2f%%",
                    100.0 * (Optimum - LowerBound) / Optimum);
        if (!PiFile)
            printff(", Ascent time = %0.2f sec.",
                    fabs(GetTime() - EntryTime));
        printff("\n");
    }
    MaxAlpha = (GainType) fabs(Excess * Cost); // 得到最大的alpha阈值，用于与maxCandidate共同限制候选集大小
    if ((A = Optimum * Precision - Cost) > 0 && A < MaxAlpha)
        MaxAlpha = A;
    if (CandidateSetType == DELAUNAY || MaxCandidates == 0)
        OrderCandidateSet(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
    else
        GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);

  End_CreateCandidateSet:
    if (ExtraCandidates > 0) {
        AddExtraCandidates(ExtraCandidates,
                           ExtraCandidateSetType,
                           ExtraCandidateSetSymmetric);
        AddTourCandidates();
    }
    ResetCandidateSet(); // 在每个节点的候选集内部进行排序，并去除不可能的候选点
    if (MaxTrials > 0 ||
        (InitialTourAlgorithm != SIERPINSKI &&
         InitialTourAlgorithm != MOORE)) {
        Na = FirstNode;
        do {
            if (!Na->CandidateSet || !Na->CandidateSet[0].To) {
                if (MaxCandidates == 0)
                    eprintf("MAX_CANDIDATES = 0: Node %d has no candidates",
                            Na->Id);
                else
                    eprintf("Node %d has no candidates", Na->Id);
            }
        }
        while ((Na = Na->Suc) != FirstNode);
        if (!CandidatesRead && SubproblemSize == 0)
            WriteCandidates();
    }
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do
            for (i = 1; i < Na->Id; i++)
                Na->C[i] += Na->Pi + NodeSet[i].Pi;
        while ((Na = Na->Suc) != FirstNode);
    }
    if (TraceLevel >= 1) {
        CandidateReport();
        printff("Preprocessing time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
}
