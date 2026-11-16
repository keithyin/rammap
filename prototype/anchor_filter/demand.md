
在使用 rammap 的过程中，发现了一个问题，在低复杂区域，高频重复seed未被过滤导致 anchor 组合爆炸，这个问题在基因组比对上可能没有很好的解决方法，但是在一些特殊场景中，这个问题可以基于数据先验进行缓解。
比如：pacbio ccs中，subreads到draft的比对。该比对存在两个特点：
    1. 极低的概率出现 long-join 的情况，
    2. subreads 和 draft 的比对基本是在一个对角线上进行的（两者比对不一定在 0 附近， ）。

存在低复杂区域的数据：
    /root/projects/rammap/test_data/low_complexity_target.fa
    /root/projects/rammap/test_data/low_complexity_queries.fa

低复杂区域没那么多的数据：
    /root/projects/rammap/test_data/s_aureus_target.fa
    /root/projects/rammap/test_data/s_aureus_queries.fa

