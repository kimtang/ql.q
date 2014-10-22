
\l quant.q
\l ql.q
num:100000;steps:256;
paths:.ql.paths `type_`spot`drift`diffu`matur`steps`repl! (`euler;100f;{[t;s] 0.05 *s };{[t;s] 0.2 * s};1.0;steps;num)

r:exp[-0.05] * avgs {max 0,avg[x]-110} each paths
