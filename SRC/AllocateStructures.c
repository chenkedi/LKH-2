#include "Segment.h"
#include "LKH.h"
#include "Heap.h"
#include "Sequence.h"

/*      
 * The AllocateStructures function allocates all necessary 
 * structures except nodes and candidates.
 */

#define Free(s) { free(s); s = 0; }

void AllocateStructures()
{
    int i, K;

    // 先清空所有数据结构的指针
    Free(Heap);
    Free(BestTour);
    Free(BetterTour);
    Free(HTable);
    Free(Rand);
    Free(CacheSig);
    Free(CacheVal);
    Free(T);
    Free(G);
    Free(t);
    Free(p);
    Free(q);
    Free(SwapStack);
    Free(tSaved);

    // 调用Heap.h生成符合dimension规定的二叉堆
    MakeHeap(Dimension);
    assert(BestTour = (int *) calloc(1 + Dimension, sizeof(int)));
    assert(BetterTour = (int *) calloc(1 + Dimension, sizeof(int)));
    assert(HTable = (HashTable *) malloc(sizeof(HashTable))); // HashTable是C中的哈希表
    HashInitialize((HashTable *) HTable); // 初始化HashTable，将每个EntrySet中的hash赋值为整型最大值，cost赋值为长整型最小值
    // SRandom为随机数指定seed，可用Java的Random代替。生成的随机数序列存储在Random.c定义的
    SRandom(Seed);
    // 申请随机数序列空间
    assert(Rand = (unsigned *)
           malloc((Dimension + 1) * sizeof(unsigned)));
    // 用初始化好的Random()函数填充Rand随机数组，数组维数为Dimension + 1
    for (i = 1; i <= Dimension; i++)
        Rand[i] = Random();
    SRandom(Seed); // 使用Seed再次初始化Random序列
    if (WeightType != EXPLICIT) {
        for (i = 0; (1 << i) < (Dimension << 1); i++); // 通过左移的指数递增，找到第一个比2 * Dimension 大的2的幂数的位移位数作为值
        i = 1 << i; // i 值设置为上部for循环中生成的2的幂数
        assert(CacheSig = (int *) calloc(i, sizeof(int))); // WeightType非Explicit则意味着要在算法运行的过程中计算，所以使用Cache缓存距离以免重复计算
        assert(CacheVal = (int *) calloc(i, sizeof(int))); // 这里为何用最小的大于Dimension * 2的2的幂数作为缓存的大小？？？
        CacheMask = i - 1;
    }
    AllocateSegments();
    // K 表示 k-opt中，k要>=2
    K = MoveType;
    // 第一次move的后续move的k值
    if (SubsequentMoveType > K)
        K = SubsequentMoveType;
    assert(T = (Node **) malloc((1 + 2 * K) * sizeof(Node *)));
    assert(G = (GainType *) malloc(2 * K * sizeof(GainType)));
    assert(t = (Node **) malloc(6 * K * sizeof(Node *)));
    assert(tSaved = (Node **) malloc((1 + 2 * K) * sizeof(Node *)));
    assert(p = (int *) malloc(6 * K * sizeof(int)));
    assert(q = (int *) malloc(6 * K * sizeof(int)));
    assert(incl = (int *) malloc(6 * K * sizeof(int)));
    assert(cycle = (int *) malloc(6 * K * sizeof(int)));
    assert(SwapStack =
           (SwapRecord *) malloc((MaxSwaps + 6 * K) * sizeof(SwapRecord)));
}

/*      
 * The AllocateSegments function allocates the segments of the two-level tree.
 */

void AllocateSegments()
{
    Segment *S = 0, *SPrev; // 初始化分段的指针
    SSegment *SS = 0, *SSPrev;
    int i;

    FreeSegments();
#ifdef THREE_LEVEL_TREE
    GroupSize = (int) pow((double) Dimension, 1.0 / 3.0); // 当segment的数据结构是三层数或者是双向链表时，GroupSize都不同
#elif defined TWO_LEVEL_TREE // 在Segment.h中默认定义了两层数，所以程序中默认使用该数据结构
    GroupSize = (int) sqrt((double) Dimension); // 每个segment的初始大小为根号Dimension，可以参考Two-level-tree的论文
#else
    GroupSize = Dimension;
#endif
    Groups = 0; // 当前segments的个数
    for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S) {
        assert(S = (Segment *) malloc(sizeof(Segment)));
        S->Rank = ++Groups; // 该segement在同级分段中的序号
        if (!SPrev) // 首次循环将S赋值给FirstSegment
            FirstSegment = S;
        else // 后续连接S与Sprev
            SLink(SPrev, S);
    }
    SLink(S, FirstSegment); // 最后首尾相连
#ifdef THREE_LEVEL_TREE
    SGroupSize = sqrt((double) Groups);
#else
    SGroupSize = Dimension;
#endif
    SGroups = 0; // 父segments的个数,如上同理构建第二层树的父层
    for (i = Groups, SSPrev = 0; i > 0; i -= SGroupSize, SSPrev = SS) {
        SS = (SSegment *) malloc(sizeof(SSegment));
        SS->Rank = ++SGroups;
        if (!SSPrev)
            FirstSSegment = SS;
        else
            SLink(SSPrev, SS);
    }
    SLink(SS, FirstSSegment);
}
