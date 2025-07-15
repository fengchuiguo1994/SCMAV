```
~/Tools/ChIA-PIPE-master/util/cpu-dir/cpu cluster -s 7970 -M -B 1000 -5 5,-20 -3 5,480 -t 5 -O THG0099V_bulk.dis8k.e500 -j -x -v 1 -g THG0099V_bulk.dis8000.bam 1>nr.dis8k.log 2>nr.dis8k.log.err
```
```
bamToBed -bedpe -i THG0099V_bulk.dis8000.bam | awk '{print ($6+$5-$3-$2)/2}' 
7911
7964.5
7903
7933.5
7969
7944
7941
7943
7977
7949.5
```
加了-M参数，所以从中点间的长度开始取值。cpu在做cluster时，不会去读取insert size那一列，而是去每个reads重新计算中间点之间的距离，大于-s才会被保留下来，延申，如果是正：取中点，中点-20得到起点（-5 5,-20），起点+500得到终点（-5 5,-20 -3 5,480）；如果是负：取中点，中点+20得到起点（-5 5,-20），起点-500得到终点（-5 5,-20 -3 5,480）。<br/>

另：cpu的insert size数字不是标准定义而来，而是两个片段的各自最远端的差（不考虑比对方向）。标准定义的是两个片段的末端的差（+的起点和-的终点）
```
samtools fixmate THG0099V_bulk.dis8000.bam THG0099V_bulk.dis8000.fix.bam
```