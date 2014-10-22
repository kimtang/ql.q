\l quant.q
\l ql.q
/ testing black scholes formula
tab:flip `type_`direct`spot`strike`rate`vola`matur!(`bls`vega`delta`theta`vega`rho;`call;100;100;0.01;0.25;1.0);
update price:.ql.bls tab from tab
num:100000;
tab1:([] type_:num?`bls`vega`delta`theta`vega`rho;direct:num?`call`put;spot:num?100.0;strike:num?100.0;rate:num?0.10;vola:num?0.50;matur:num?10.0 );tab1
update price:.ql.bls tab1 from tab1
ftab:{[x]([] type_:x?`bls`vega`delta`theta`vega`rho;direct:x?`call`put;spot:x?100.0;strike:x?100.0;rate:x?0.10;vola:x?0.50;matur:x?10.0 )};
/ measure time
num:2;
scal:1000000;
perf:flip `num`time!(scal*1+til num;value each "\\t .ql.bls ftab ",/: string scal*1+til num);perf


/ testing implied vola
\l ql.q
num:1000;
tab2:flip `type_`direct`spot`strike`rate`price`matur!(`impl;num?`call`put;100;100;0.01;(03.00 + 0.005*til num);1.0)
t:{[x]flip `type_`direct`spot`strike`rate`price`matur!(`impl;`call;100;100;0.01;(03.00 + 0.005*til x);1.0)}
num:300;
flip `line`time!(1+til num; value each "\\t .ql.bls t ",/: string 1+til num)
update blsprice: .ql.bls res2 from res2:update type_:`bls, vola: .ql.bls tab2 from tab2

/ testing binomial tree
\l ql.q
t:([] spot:100;rate:0.01;vola:0.25;matur:1.0;num:2500;payoff:({x};{max 0,x-100};{max 0,100-x};{abs 100-x}));t
/ update price: .ql.binbaum t from t
update price: .ql.binbaum each t from t
/ .ql.binbaum t

/ testing gaussian random number generator
\l ql.q
num:`int$1e6;
select count i by 0.01 xbar r  from ([] n:til num ;r:.ql.randn num)
\t value ".ql.randn num"
select n, sums r from ([] n:til 1+num ;r:0,.ql.randn num)

/ testing geometric brownian motion
\l ql.q
num:1000000;steps:4;
paths:.ql.paths `type_`spot`drift`diffu`matur`steps`repl!
    (`euler;100f;{[t;s] 0.01 *s };{[t;s] 0.25 * s};1.0;steps;num)
select count i by 1 xbar paths from t:([] paths:last each paths )
exp[-0.01] * avg {max 0,last[x]-100} each paths
exp[-0.01] * avg {max 0,neg last[x]-100} each paths
exp[-0.01] * avg {last[x]} each paths

/ testing simulate multidimensional sde
\l ql.q
num:1000000
L::(0.25 0f; 0 0.25)
steps:5;
arg:`type_`spot`drift`diffu`matur`steps`repl!(`euler;(100 100f);{[t;s] 0.01 *s };{[t;s] s*L};1.0;steps;num)
p:.ql.paths arg
select count i by 0.5 xbar p from ([] p:last each raze p) 
select count i by 0.5 xbar p from ([] p:last flip raze p) / this version is slightly faster
1_select count i by 1 xbar p from ([]p:raze {[x] {[x]max 0,x-100}each last x} each raze p)
1_select count i by 1 xbar p from ([]p:raze {[x] {[x]max 0,100-x}each last x} each raze p)
