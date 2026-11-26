/*
 * solver.c - Cache Lab Bonus: reverse-engineer an unknown cache.
 *
 * Known ranges: b in 3..7, s in 1..9, E in {1,2,3,4,5,6,8,12,16}, three
 * replacement policies, write-allocate on or off.  Every probe is a hit or
 * a miss, so each parameter is pinned down with a run of loads whose
 * expected outcome is the same for every cache that shares that parameter.
 *
 *   b             the first power-of-two offset that misses block 0
 *   E             capacity of one set, measured with E+1 blocks that all
 *                 map to set 0 (stride 512 blocks is a multiple of every
 *                 possible 2^s because s <= 9)
 *   s             the smallest stride 2^k whose E+1 blocks collide
 *   write_allocate  store then load the same fresh block
 *   policy        fill one set, touch its oldest block, insert one more
 *                 block, then look for the touched block: it survives
 *                 under LRU, is gone under FIFO, and RANDOM mixes the two
 *
 * A set's capacity is tested by filling it and then reading every block
 * back once.  Under any policy a batch of at most E blocks stays resident,
 * so the whole read-back hits; a batch of more than E blocks leaves at
 * least one block out, so the read-back must miss.  Hits never evict
 * anything, so the read-back itself does not disturb the answer.
 */

#include "solver.h"

#define MAX_S 9
#define STRIDE_E 512u      /* 2^9 blocks: same set for every legal s */

/* Byte address of block number blk. */
static uint32_t block_addr(uint32_t blk, int b)
{
    return blk << b;
}

/* Return 1 when e blocks that all map to set 0 stay resident together. */
static int capacity_ok(int e, int b)
{
    int i;

    cache_flush();
    for (i = 0; i < e; i++)
        cache_load(block_addr((uint32_t)i * STRIDE_E, b));
    for (i = 0; i < e; i++)
        if (!cache_load(block_addr((uint32_t)i * STRIDE_E, b)))
            return 0;
    return 1;
}

/* Return 1 when e + 1 blocks spaced 2^k apart cannot all stay resident,
 * i.e. when 2^k is a multiple of the number of sets. */
static int stride_conflict(int k, int e, int b)
{
    int i;
    uint32_t stride = 1u << k;

    cache_flush();
    for (i = 0; i <= e; i++)
        cache_load(block_addr((uint32_t)i * stride, b));
    for (i = 0; i <= e; i++)
        if (!cache_load(block_addr((uint32_t)i * stride, b)))
            return 1;
    return 0;
}

/* Does the block last touched in a filled set survive one extra insert? */
static int victim_test(int e, int b)
{
    int i;

    cache_flush();
    for (i = 0; i < e; i++)
        cache_load(block_addr((uint32_t)i * STRIDE_E, b));
    cache_load(block_addr(0, b));                       /* worst case for LRU */
    cache_load(block_addr((uint32_t)e * STRIDE_E, b));  /* forces an eviction */
    return cache_load(block_addr(0, b));
}

void solve(void)
{
    static const int cand[] = { 1, 2, 3, 4, 5, 6, 8, 12, 16 };
    int b = 0;
    int s = 1;
    int e = 1;
    int wa = 0;
    int policy = POLICY_LRU;
    int i, k, hits, misses, reps;

    /* Block size: 1 << b is the first offset past block 0. */
    cache_flush();
    cache_load(0);
    while (b < 31 && cache_load(1u << b))
        b++;

    /* Associativity: the largest candidate that still fits in one set. */
    for (i = 0; i < (int)(sizeof(cand) / sizeof(cand[0])); i++)
        if (capacity_ok(cand[i], b))
            e = cand[i];

    /* Set index bits: first stride that overflows a single set. */
    for (k = 1; k <= MAX_S; k++) {
        if (stride_conflict(k, e, b)) {
            s = k;
            break;
        }
    }

    /* Write allocate: a store miss either brings the block in or not. */
    cache_flush();
    cache_store(block_addr(0, b));
    wa = cache_load(block_addr(0, b)) ? 1 : 0;

    /* Replacement policy: ask the same question many times. */
    if (e > 1) {
        reps = 40 * e;
        if (reps < 200)
            reps = 200;
        hits = 0;
        misses = 0;
        for (i = 0; i < reps; i++) {
            if (victim_test(e, b))
                hits++;
            else
                misses++;
        }
        if (misses == 0)
            policy = POLICY_LRU;
        else if (hits == 0)
            policy = POLICY_FIFO;
        else
            policy = POLICY_RANDOM;
    }

    report(b, s, e, policy, wa);
}
