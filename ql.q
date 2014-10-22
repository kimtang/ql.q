/ the file quant.q is needed from gordon baker
/ gbkr.com
/\l quant.q

\d .ql_impl

/------- Finance 
phi:{[x] exp[neg (x xexp 2)%2f] % sqrt 2 * .quant.pi}; / black scholes formula  

bls:()!(); 

bls[`d]:{[x] d1:(log[x[`spot] % x[`strike]] + x[`matur]*x[`rate]+0.5*x[`vola] xexp 2) %x[`vola] * sqrt x[`matur];
	d2:d1-x[`vola] *sqrt x[`matur];
	sig:neg 1f-2f*`float$(x[`direct]=`call );
	(`d1`d2`sig)!(d1;d2;sig) };

bls[`bls]:{[x] neg x[`sig]*(  x[`strike]*exp[neg x[`rate]*x[`matur]] *.quant.cdf[`gauss] x[`sig]*x[`d2] )-x[`spot]*.quant.cdf[`gauss] x[`sig]*x[`d1]};

bls[`delta]:{[x] x[`sig]*.quant.cdf[`gauss] x[`sig]*x[`d1] };
bls[`theta]:{[x] nd2:x[`strike]*x[`rate]*.quant.cdf[`gauss] x[`sig]*x[`d2];
	nd1:%[neg x[`vola]*x[`spot]*phi x[`d1];2f*sqrt x[`matur]];
	:nd1+x[`sig]*nd2*exp neg x[`matur]};
bls[`gamma]:{[x] :phi[x[`d1]] % x[`spot]*x[`vola]*sqrt x[`matur]   };

bls[`vega]: {[x] :phi[x[`d1]]*x[`spot]*sqrt x[`matur] };

bls[`rho]:{[x] :x[`strike]*x[`matur]*exp[neg x[`rate] *x[`matur] ] *.quant.cdf[`gauss] x[`sig]*x[`d2]  };

bls[`98]:{[x] a:update num :til count x from x ; impl:select from a where type_=`impl; nimpl:select from a where not type_=`impl;
              t:$[0=count nimpl;nimpl;flip flip[nimpl],'bls[`d]nimpl] ;
              g:$[0=count t;nimpl;t group t`type_];
              nimpl:$[0=count t;nimpl;raze key [g] {[x;y] p:bls[x] y ; update p from y}' value g];
              impl:update p:.ql.bls each impl from impl;
              :exec p from `num xasc nimpl,impl
    };
/ dictionary
bls[`99]:{[x] 
    :$[0<sum 0<type each value x;.ql.bls flip x;
	   $[`impl=x`type_; .ql_impl.bls[`impl] x; .ql_impl.bls[x`type_] x,.ql_impl.bls[`d] x ]]
    };

bls[`impl]:{ [x] p:`type_ _ x;
    www:{[p;x] 1e-10< abs p[`price]-.ql.bls p,(`type_`vola)!(`bls;x)  } p ;
    /ww:{[x] abs p[`price]-.ql.bls p,(`type_`vola)!(`bls;x)};
    {[p;x] x+(p[`price]-.ql.bls p,(`type_`vola)!(`bls;x))  % .ql.bls p,(`type_`vola)!(`vega;x) }[p]/[www;0.65]
    };

/ binomial tree
binbaum:()!();
binbaum[99h]:{[x]
    s:x[`spot];r:x[`rate];v:x[`vola];t:x[`matur];n:x[`num];f:x[`payoff];
    dt: t % n;
    beta: avg exp dt*(0f;v xexp 2)+(neg r;  r);
    u: beta + sqrt neg 1-beta xexp 2;
    d: reciprocal u;
    p: (neg d-exp r*dt) % u-d;
    q:1-p;
    S: s */{[u;d;x](u xexp x;d xexp reverse x)}[u;d] til n;
    /V:{max 0,k-x} each S;
    V:f each S;
    exp[neg r*t]*first {[p;q;x](p*1_x )+ q*-1_x}[p;q]/[{not 1=count x};V]
    };
binbaum[98h]:{
    break "dont use has bug";;
    s:x[`spot];r:x[`rate];v:x[`vola];t:x[`matur];n:x[`num];f:x[`payoff];
    dt: t % n;
    beta: avg exp dt*/:(0f;v xexp 2)+(neg r;  r);
    u:: beta + sqrt neg 1-beta xexp 2;
    d:: reciprocal u;
    p:: (neg d-exp r*dt) % u-d;
    q::1-p;
    S: s{[x;y] nn:y[0];uu::y[1];dd::y[2]; x */ {[x](uu xexp x;dd xexp reverse x)} til nn }' flip(n;u;d);
    V:f{[x;y] x y}'S;    
    exp[neg r*t] * first flip{ (p*flip[1_flip[x]] )+ q*flip[-1_flip[x]] }/[{not 1=max count each x};V];
    };
binbaum[0h]:{.ql.binbaum[99h;`spot`rate`vola`matur`num`payoff!x til 6 ]};

/ TODO: change to the implementation from stat.q
randn:()!();
randn[1]: {hn:`int$ x % 2;
		   x#raze exec  s*f,s*g from select s:sqrt -2.0*log u1,f:cos .quant.pi2 * u2,g:sin .quant.pi2 * u2 from ([]u1:hn?1.0;u2:hn?1.0)
		   };
/ old implementation, this is also very bad randn[1]: {:x#1_raze{ u*sqrt -2f*log[s] % s:{u$u::-1+2?2.0}/[1<;2]}\[`int$ x % 2 ;0] };
/ old implementation, which is very bad: randn[2]: {a1::x[0];a2:x[1]; flip {randn[1] a1} each til a2 };  
randn[2]: {a1:x[0];a2:x[1]; (a1,a2) # randn[1] a1*a2 };


ostyle:()!();
ostyle[`native]:{x};
ostyle[`wtime]:{[x] t: flip[x] 0; p:flip[x] 1;
                $[ 1 = count first[p] 0 ;(t; flip p) ;(t;flip each flip p)] };
ostyle[`opaths]:{[x] ostyle[`wtime;x] 1};

pathsde:()!();
pathsde[`euler]:{[x] 
    spot:x[`spot];drift:x[`drift];diffu:x[`diffu];matur:x[`matur];steps:x[`steps]-1;repl:x[`repl];
    h:matur % steps;
    delta:sqrt h;
    / get rid of the time vector
    :{ [drift;diffu;delta;repl;h;x] x + (h;(drift[x 0; x 1] * h) +diffu[x 0; x 1]*delta * .ql.randn repl) }[drift;diffu;delta;repl;h]\[steps;(0f;repl#spot)]
    };

/pathsde[`milstein]:{[x] 
/    spot:x[`spot];drift::x[`drift];diffu::x[`diffu];matur:x[`matur];steps:x[`steps]-1;repl::x[`repl];
/    h::matur % steps;
/    delta::sqrt h;
    / get rid of the time vector
/    {[x] diffu[x 0; x 1] * delta * .ql.randn repl }\[steps;(of;repl#spot)]
/    flip first 1_flip { [x] x + (h;(drift[x 0; x 1] * h) +diffu[x 0; x 1]*delta * .ql.randn repl) }\[steps;(0f;repl#spot)]
/    };

pathsde[`runge]:{[x] 
    spot:x[`spot];drift:x[`drift];diffu:x[`diffu];matur:x[`matur];steps:x[`steps]-1;repl:x[`repl];
    h:matur % steps;
    delta:sqrt h;
    / get rid of the time vector
    :{ [drift;diffu;delta;repl;h;x]
        deltaw:delta * .ql.randn repl;
        yhat:x[1]+ (drift[x 0; x 1] * h )+diffu[x 0; x 1]*delta;
        x + (h;(drift[x 0; x 1] * h )+ (diffu[x 0; x 1]*deltaw) + reciprocal[2f*delta] * ( diffu[x 0;yhat]-diffu[x 0; x 1] ) * neg[h - deltaw*deltaw])
        }[drift;diffu;delta;repl;h]\[steps;(0f;repl#spot)]
    };

/ paths:{ [x] spot::x[`spot];r:x[`rate];v:x[`vola];t:x[`matur];st:x[`steps];repl:x[`repl];
/    dt:t % st;
/    nudt: dt* r- 0.5* v xexp 2;
/    sidt: v * sqrt dt;
/    inc:nudt+sidt*.ql.randn (repl;st);
/    logp: sums each {log[spot],x} each inc;
/    : exp each logp
/    };

/ here we will implement the multidimensional sde
pathmde:()!();
pathmde[`euler]:{[x] 
    spot:x[`spot];drift:x[`drift];diffu:x[`diffu];matur:x[`matur];steps:x[`steps]-1;repl:x[`repl];
    d:count spot;h:matur % steps;delta:sqrt h;    
    :{[spot;drift;diffu;repl;d;h;delta;x] 
        current:x 1;t:x 0;wnt:delta*(repl,d) # .ql.randn d*repl;
        x+(h;current{[drift;diffu;t;h;x;y] (drift[t;x]*h) +diffu[t;x] wsum\: y }[drift;diffu;t;h]'wnt)
    }[spot;drift;diffu;repl;d;h;delta]\[steps;(0f;repl#enlist spot)]
    };

pathmde[`neuler]:{[x] 
    spot::x[`spot];drift::x[`drift];diffu::x[`diffu];matur:x[`matur];steps:x[`steps]-1;repl::x[`repl];
    d::count spot;h::matur % steps;delta::sqrt h; wnt: delta*flip steps cut d cut .ql.randn d*repl*steps;
    str:(0f;repl#enlist spot);
    :enlist[str],{[x;y] cur:x 1;t::x 0;x+(h;(drift[t;cur]*h) +diffu[t;cur;y])}\[str;wnt]};


mc:{[x]:((),x`payoffs)@\:/:.ql.paths x};
points:()!();
points[`gkps7]:`w`p!{[x] w:(0.129484966168870;0.279705391489277;0.381830050505119;0.417959183673469);
                   p:(0.949107912342759;0.741531185599394;0.405845151377397);
                   w:w,1_reverse w;p: p,0f,neg reverse p;
                   :(w;p)
              }[];
points[`gkps15]:`w`p!{[x] w:(0.022935322010529;0.063092092629979;0.104790010322250;0.140653259715525;0.169004726639267;0.190350578064785;0.204432940075298;0.209482141084728);
                    p:(0.991455371120813;0.949107912342759;0.864864423359769;0.741531185599394;0.586087235467691;0.405845151377397;0.207784955007898);
                    w:w,1_reverse w;p: p,0f,neg reverse p;
                    :(w;p)
              }[];

integrate:{ [x] f:x[`f];a:x[`a];b:x[`b];tol:x[`tol]; hl:%[b-a; 2f];cer:%[b+a; 2f];
    ff:f cer + points[`gkps15;`p]* hl;
    large: hl* sum points[`gkps15;`w] * ff;
    small: hl* sum points[`gkps7;`w] * ff 1+ 2*til 7 ;
    :$[ tol> abs large - small ;:large; :sum integrate each flip `f`a`b`tol!(f;(a;cer);(cer;b);tol % 2f) ]
    };

\d .
\d .ql

bls:{ [x]
    :.ql_impl.bls[`$ string type x] x
    /xx:?[`impl =x`type_;{x};.ql_impl.bls[`dm] ][x] ;.ql_impl.bls[xx`type_]xx
    };
binbaum:{[x] .ql_impl.binbaum[type x;x]};

randn:{[x] :.ql_impl.randn[count x] x };

/ mvnrnd[0 0 0f; (1 -0.5 0; -0.5 1 -0; 0 0 0.2); 10]
mvnrnd:{[m;v;n] m + cholcov[v] mmu n cut randn n* count m}

paths:{[x] 
    res:$[1=count x`spot;.ql_impl.pathsde[x`type_] x;.ql_impl.pathmde[x`type_] x];
    style:$[` = x`output;`opaths;x`output];
    :.ql_impl.ostyle[style] res};

/ cholesky decomposition
cholcov:{ f:{[x;y] A:x;i:y 0; j:y 1; s:(A . y) - sum ((*) . A y) til j; A[i;j]:$[i>j;s % A[j;j];(i=j) & s>0f;sqrt s;'`not_pos_def];:A};
    r:f over enlist[x] , raze (tn,/:\:tn)@'where each wn:tn>=\:tn:til n:count x;
    :flip r* wn;}

/
 l:(1 2 3 4f;0 2 3 4f; 0 0 3 4f; 0 0 0 4f);
 a:l mmu flip l;
 cholcov a;
\

mc:{[x] .ql_impl.mc x}
integrate: .ql_impl.integrate
solver:{ raze lsq[enlist y;flip x] }
diag:{x @'til count x}
zeros:{[x;y] (x,y)#*[x;y]#0f }
ones:{[x] til[x]{[x;y] y[x]:1f; :y}'zeros[x;x] }

/ covarinace matrix
cvm:{(x+flip(not n=\:n)*x:(n#'0.0),'(x$/:'(n:til count x)_\:x)%count first x)-a*\:a:avg each x}

/ correlation matrix
crm:{cvm[x]%u*/:u:dev each x}
\d .
