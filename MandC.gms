




Sets
i index of generating(NG) units /1*8/
t index of hour /1*24/
d inex of days /1*4/
y index of year /y1*y10/
e index of ESS systems /1*3/
pv set of solar pv panels /1*3/
N  all nodes /n1*n24/
N_ex(N) load-points /n5*n12/
N_can(N)   future load-points /n3,n4,n13*n24/
S(N) all substation nodes /n1*n4/
S_ex(S)  existing substation nodes /n1, n2 /
S_can(S) future substation nodes /n3, n4 /
np(N) nodes of only load /n1,n2,n5*n24/
L all lines /l1*l38/
L_ex(L)  existing-Lines /l1*l9/
L_can(L)   candidate-Lines /l10*l38/
Nm(np) nodes that have microgrid / n20*n24 /
nodePars node parameters /Psmax1,Qsmax1,NIC1,Psmax2,Qsmax2,NIC2,
Vmax,Vmin,Psmax,Qsmax,PsmaxR,QsmaxR,NRC,state /
LinePars Line parameters /fbus,tbus,R,X,PFmax,QFmax,Rr,Xr,PFRmax,QFRmax,LRC,R1,X1,
PFI1max,QFI1max,LIC1,R2,X2,PFI2max,QFI2max,LIC2,R3,X3,PFI3max,QFI3max,LIC3,InitState,Length/
DG_par /Pmin,Pmax,IC,MC,SUC,SDC,RU,RD,m /
ES_par /PR,IC,RC,RDC,Einit,Emin,Emax /
PV_par /cap,IC/ ;
alias (t,j);
Alias(tt,t);

Parameters
number_of_nodes(y) /
y1       11
y2       12
y3       13
y4       15
y5       18
y6       20
y7       22
y8       23
y9       23
y10      23 /
numn(n) /
n1  1
n2  2
n3  3
n4  4
n5  5
n6  6
n7  7
n8  8
n9  9
n10 10
n11 11
n12 12
n13 13
n14 14
n15 15
n16 16
n17 17
n18 18
n19 19
n20 20
n21 21
n22 22
n23 23
n24 24
/
YY(y) numerical value of each year
        / y1      1
          y2      2
          y3      3
          y4      4
          y5      5
          y6      6
          y7      7
          y8      8
          y9      9
          y10     10 /
node_data(n,nodePars)
Line_data(L,LinePars)
LoadP_data(n,y)
LoadQ_data(n,y)
A(L,N)
DG_param(i,DG_par)
ES_param(e,ES_par)
PV_param(pv,PV_par)
Ps(y,d,t) output power of pv cell equal to efficiency times global irradiance at time t
LMP(y,d,t) price of energy at each hour t
Pd(y,d,t) Load demand at time t

$ call gdxxrw node_data.xlsx par node_data rng=sheet1!A1:N25 rdim=1 cdim=1
$GDXIN node_data.gdx
$load node_data
$GDXIN
$ call gdxxrw Line_data.xlsx par Line_data rng=sheet1!A1:AC39 rdim=1 cdim=1
$GDXIN Line_data.gdx
$load Line_data
$GDXIN
$ call gdxxrw LoadP.xlsx par LoadP_data rng=sheet1!A1:K25 rdim=1 cdim=1
$GDXIN LoadP.gdx
$load LoadP_data
$GDXIN
$ call gdxxrw LoadQ.xlsx par LoadQ_data rng=sheet1!A1:K25 rdim=1 cdim=1
$GDXIN LoadQ.gdx
$load LoadQ_data
$GDXIN
$ call gdxxrw Inc_matrix.xlsx par A rng=sheet1!A1:Y39 rdim=1 cdim=1
$GDXIN Inc_matrix.gdx
$load A
$GDXIN

$ call gdxxrw LMPs.xlsx par LMP rng=sheet3!A1:AA41 rdim=2 cdim=1
$ gdxin LMPs.gdx
$ load LMP
$ gdxin
$ call gdxxrw Solar.xlsx par Ps rng=Psolar!A1:AA41 rdim=2 cdim=1
$ gdxin Solar.gdx
$ load Ps
$ gdxin
$ call gdxxrw Demand.xlsx par Pd rng=Sheet1!A1:AA41 rdim=2 cdim=1
$ gdxin Demand.gdx
$ load Pd
$ gdxin
$ call gdxxrw DG_param.xlsx par DG_param rng=sheet1!A1:J7 rdim=1 cdim=1
$ gdxin DG_param.gdx
$ load DG_param
$ gdxin
$ call gdxxrw ES_param.xlsx par ES_param rng=sheet1!A1:H7 rdim=1 cdim=1
$ gdxin ES_param.gdx
$ load ES_param
$ gdxin
$ call gdxxrw PV_param.xlsx par PV_param rng=sheet1!A1:C7 rdim=1 cdim=1
$ gdxin PV_param.gdx
$ load PV_param
$ gdxin

Scalar
LCC Load Curtailment Cost /10000/
Ce Charge effiency / 0.95 /
DCe Discharge effiency / 0.95 /
M big number /1000/
Vref KV /20/
Lf Loss factor /0.35/
LossCost_C dolar per Mwh /60/
LOLEcost value of loss load dolar per Mwh /30000/
SubOC Subsationg operation cost dolar per MVA /2000/
roi rate of interest /.1/
LoadFactor /0.8/  ;

Variables
DNObjFunc
InvCost(y)
OprCost(y)
LossCost(y)
PF(L,y)
QF(L,y)
Pprod(n,y)
Qprod(n,y)
rp(n,y)
rq(n,y)
V(n,y)
MGinvcost(Nm,y) Annual cost of investments yearly
MGIncome(Nm,y) from selling power to Network yearly
MGOpCost(Nm,y) Operation cost yearly
MGNPV(Nm) Net Present Value of total costs
 ;
positive variables
Pprod(n,y)
Qprod(n,y)
rp(n,y)
rq(n,y)
V(n,y)
Revenue(y)
uu(l,y)
vv(l,y)
PFlin(l,y)
SU(Nm,i,y,d,t) Startup cost of unit i at time t
SD(Nm,i,y,d,t) Shutdown cost of unit i at time t
PMG(Nm,i,y,d,t) Real Power scheduled for unit i at time t
LC(Nm,y,d,t) Load Curtailment
En(Nm,e,y,d,t) Energy stored in ESS system e at time t
C(Nm,e,y,d,t) Charge power of ESS system e at time t
Dc(Nm,e,y,d,t) discharge power of ESS system e at time t
Ppur(Nm,y,d,t) P purchased from main Network at the price of LMP
Psale(Nm,y,d,t) P sold to Network at price of LMP
 ;

BINARY Variables
XI1(L_can,y)
XI2(L_can,y)
XI3(L_can,y)
XR(L_ex,y)
Z1(S_can,y)
Z2(S_can,y)
ZR(np,y)
b(l,y)
II(Nm,i,y,d,t) Commitment state of unit i at time t
SII(Nm,pv,y,d,t) commit of pv
SMG(Nm,i,y) planning state of each unit
SI(Nm,pv,y) planning state of each pv cell
ES(Nm,e,y) planning state of ESS
U(Nm,e,y,d,t) 0 or 1 for charge or discharge of ESS
;

equations
DNTotCost Total Cost
eq_InvCost
eq_OprCost
eq_LossCost
eq_lin1
eq_lin2
eq_lin3
eq_lin4

KCL

KVL1
KVL21
KVL22
KVL31
KVL32
KVL41
KVL42

eq_powerflow1
eq_powerflow2
eq_Qflow

eq_sub_ins
eq_OnlySub
eq_exSub1
eq_exSub2

eq_Line_rep
eq_only_ins
eq_Line_ins1
eq_Line_ins2
eq_Line_ins3
eq_ZOnce1
eq_ZOnce2
eq_OnlySubIns
eq_ZR
eq_ZROnce

eq_referenceV

Radiality1
Radiality2

obj_func   objective function
prim_Aninvcost        annual cost
prim_state_commit
gen_min
gen_max
*SU_eq
*SD_eq
rampup_limit
rampdn_limit
eq_pvcommit
ESS_balance
ESS_init
ESS_min
ESS_max
discharge_limit
charge_limit
discharge_rate_min
discharge_rate_max
charge_rate_min
charge_rate_max
power_balance
eq_Psale
eq_Ppur
Reserve
Income_eq
OpCost_eq

;


DNTotCost.. DNObjFunc =e= sum(y, ( InvCost(y) + OprCost(y) + LossCost(y))/((1+roi)**YY(y)))  ;

eq_InvCost(y)..  InvCost(y) =e= sum(S_can, (Z1(S_can,y)-Z1(S_can,y-1))*node_data(S_can,'NIC1') + (Z2(S_can,y)-Z2(S_can,y-1))*node_data(S_can,'NIC2'))
+sum(np,(ZR(np,y)-ZR(np,y-1))*node_data(np,'NRC'))+ sum(L_ex,(XR(L_ex,y)-XR(L_ex,y-1))*Line_data(L_ex,"LRC"))
+ sum(L_can,(XI1(L_can,y)-XI1(L_can,y-1))*Line_data(L_can,"LIC1"))
+ sum(L_can,(XI2(L_can,y)-XI2(L_can,y-1))*Line_data(L_can,"LIC2"))
+ sum(L_can,(XI3(L_can,y)-XI3(L_can,y-1))*Line_data(L_can,"LIC3")) ;

eq_OprCost(y).. OprCost(y) =e= sum(S_can,Z1(S_can,y)*node_data(S_can,'Psmax1')*SubOC + Z2(S_can,y)*node_data(S_can,'PSmax2')*SubOC)
+ sum(S_ex,node_data(S_ex,'Psmax')*SubOC) +  sum(n, rp(n,y)*LOLEcost) ;

eq_LossCost(y).. LossCost(y) =e= sum(L_ex, Lf*PFlin(L_ex,y)*Line_data(L_ex,'R')*LossCost_C)
+ sum(L_can, Lf*PFlin(L_can,y)*Line_data(L_can,'R')*LossCost_C) ;
*eq_Salvage.. SalvageCost =e= sum(S_ex,XR(S_ex,'y7')) + sum(S_can,XI(S_can,'y7'))
*eq_Revenue(y).. Revenue(y) =e= price_en(y)*sum(n,LoadP_data(n,y) - rp(n,y)) ;

eq_lin1(l,y).. PFlin(l,y) =e= uu(l,y) + vv(l,y);
eq_lin2(l,y).. PF(l,y) =e= uu(l,y) - vv(l,y) ;
eq_lin3(l,y).. uu(l,y) =l= M*b(l,y) ;
eq_lin4(l,y).. vv(l,y) =l= M*(1-b(l,y)) ;

KCL(n,y).. Pprod(n,y) + rp(n,y) +  sum(L$(Line_data(L,'tbus')=numn(n)),PF(L,y)) =e=  sum(L$(Line_data(L,'fbus')=numn(n)),PF(L,y)) + LoadP_data(n,y) ;

*KVL:   constraints are designed in a way to be relaxed when line is not in use.
KVL1(L_ex,y).. sum(n,A(L_ex,n)*V(n,y)) =e= -(PF(L_ex,y)*Line_data(L_ex,"R") + QF(L_ex,y)*Line_data(L_ex,"X"))/Vref ;
KVL21(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R1") + QF(L_can,y)*Line_data(L_can,"X1"))/Vref =g= (XI1(L_can,y) - 1)*M  ;
KVL22(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R1") + QF(L_can,y)*Line_data(L_can,"X1"))/Vref =l= (1 - XI1(L_can,y))*M  ;
KVL31(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R2") + QF(L_can,y)*Line_data(L_can,"X2"))/Vref =g= (XI2(L_can,y) - 1)*M  ;
KVL32(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R2") + QF(L_can,y)*Line_data(L_can,"X2"))/Vref =l= (1 - XI2(L_can,y))*M  ;
KVL41(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R3") + QF(L_can,y)*Line_data(L_can,"X3"))/Vref =g= (XI3(L_can,y) - 1)*M  ;
KVL42(L_can,y).. sum(n,A(L_can,n)*V(n,y)) + (PF(L_can,y)*Line_data(L_can,"R3") + QF(L_can,y)*Line_data(L_can,"X3"))/Vref =l= (1 - XI3(L_can,y))*M  ;

*power flow limitations
eq_powerflow1(L_ex,y).. PFlin(L_ex,y) =l= Line_data(L_ex,'PFmax') *(1 - XR(L_ex,y)) + Line_data(L_ex,'PFRmax')*XR(L_ex,y)  ;
eq_powerflow2(L_can,y).. PFlin(L_can,y) =l= Line_data(L_can,'PFI1max')*XI1(L_can,y) + Line_data(L_can,'PFI2max')*XI2(L_can,y)
+ Line_data(L_can,'PFI3max')*XI3(L_can,y) ;
eq_Qflow(L,y).. QF(L,y) =e= 0.75*PFlin(L,y) ;

*to model Substation placing. Ppmax is zero for non-substations.

eq_sub_ins(S_can,y).. Pprod(S_can,y) =l=  Z1(S_can,y)*node_data(S_can,'Psmax1') + Z2(S_can,y)*node_data(S_can,'Psmax2') ;
eq_OnlySub(S_can,y).. Z1(S_can,y) + Z2(S_can,y) =l= 1 ;
eq_exSub1(np,y).. Pprod(np,y) =l= (1 - ZR(np,y))*node_data(np,'Psmax') + ZR(np,y)*node_data(np,'PsmaxR')  ;
eq_exSub2(n,y).. Qprod(n,y) =e= 0.75*Pprod(n,y)  ;

*to ensure only one replacement or installation
eq_Line_rep(L_ex,y).. XR(L_ex,y) =g= XR(L_ex,y-1) ;
eq_Line_ins1(L_can,y).. XI1(L_can,y) =g= XI1(L_can,y-1) ;
eq_Line_ins2(L_can,y).. XI2(L_can,y) =g= XI2(L_can,y-1) ;
eq_Line_ins3(L_can,y).. XI3(L_can,y) =g= XI3(L_can,y-1) ;
eq_only_ins(L_can,y).. XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y) =l= 1 ;

eq_ZOnce1(S_can,y).. Z1(S_can,y) =g= Z1(S_can,y-1) ;
eq_ZOnce2(S_can,y).. Z2(S_can,y) =g= Z2(S_can,y-1) ;
eq_OnlySubIns(S_can,y).. Z1(S_can,y) + Z2(S_can,y) =l= 1 ;

eq_ZR(np,y).. ZR(np,y) =l= node_data(np,'state') ;
eq_ZROnce(np,y).. ZR(np,y) =g= ZR(np,y-1) ;


*Node Voltage limitations
eq_referenceV(S,y).. V(S,y) =e= Vref ;
*eq_Vmax(np,y).. V(np,y) =l= Node_data(np,"Vmax") ;
*eq_Vmin(np,y).. V(np,y) =g= Node_data(np,"Vmin") ;


Radiality1(y).. card(L_ex) + sum(L_can, XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y)) =e= number_of_nodes(y) -card(S_ex) ;
Radiality2(n,y).. sum(L_can$(Line_data(L_can,"tbus")=ord(n)), XI1(L_can,y) + XI2(L_can,y) + XI3(L_can,y))
+ sum(L_ex$(Line_data(L_ex,"tbus")=ord(n)), Line_data(L_ex,'InitState')) =l= 1 ;
******************************************************************************************
obj_func(Nm).. NPV(Nm) =e= sum(y,(MGInvcost(y) + MGOpCost(y) - MGIncome(y))/((1+roi)**YY(y))) ;
prim_Aninvcost(nm,y).. Aninvcost(nm,y) =e= sum(i, DG_param(i,'IC')*SMG(nm,i)) + sum(pv, PV_param(pv,'IC')*SI(nm,pv)) + sum(e,ES(nm,e)*ES_param(e,'IC')) ;

prim_state_commit(nm,i,y,d,t).. II(nm,i,y,d,t) =l= SMG(nm,i) ;

gen_min(nm,i,y,d,t).. PMG(nm,i,y,d,t) =g= DG_param(i,'Pmin')*II(nm,i,y,d,t) ;
gen_max(nm,i,y,d,t).. PMG(nm,i,y,d,t) =l= DG_param(i,'Pmax')*II(nm,i,y,d,t) ;

SU_eq(nm,i,y,d,t).. SU(nm,i,y,d,t) =g= DG_param(i,'SUC')*(II(nm,i,y,d,t) - II(nm,i,y,d,t-1)) ;
SD_eq(nm,i,y,d,t).. SD(nm,i,y,d,t) =g= DG_param(i,'SDC')*(II(nm,i,y,d,t-1) - II(nm,i,y,d,t)) ;

eq_pvcommit(Nm,pv,y,d,t).. SII(Nm,pv,y,d,t) =l= SI(Nm,pv) ;
*rampup_limit(Nm,i,y,d,t).. P(i,y,d,t) - P(i,y,d,t-1) =l= DG_param(i,'RU') ;
*rampdn_limit(Nm,i,y,d,t).. P(i,y,d,t-1) - P(i,y,d,t) =l= DG_param(i,'RD') ;

ESS_init(Nm,e,y,d,t)$(ord(t) = 1).. En(Nm,e,y,d,t) =e= ES(Nm,e)*ES_param(e,'Einit') ;
ESS_balance(Nm,e,y,d,t).. En(Nm,e,y,d,t+1) =e= En(Nm,e,y,d,t) + Ce*C(Nm,e,y,d,t) - Dc(Nm,e,y,d,t)/DCe ;
ESS_min(Nm,e,y,d,t).. En(Nm,e,y,d,t) =g= ES_param(e,'Emin')*ES(Nm,e) ;
ESS_max(Nm,e,y,d,t).. En(Nm,e,y,d,t) =l= ES_param(e,'Emax')*ES(Nm,e) ;

discharge_limit(Nm,e,y,d,t).. Dc(Nm,e,y,d,t) =l= U(Nm,e,y,d,t)*ES_param(e,'PR') ;
charge_limit(Nm,e,y,d,t).. C(Nm,e,y,d,t) =l= (1 - U(Nm,e,y,d,t))*ES_param(e,'PR') ;

discharge_rate_min(Nm,e,y,d,t).. Dc(Nm,e,y,d,t) - Dc(Nm,e,y,d,t-1) =l= ES_param(e,'RDC') ;
discharge_rate_max(Nm,e,y,d,t).. Dc(Nm,e,y,d,t) - Dc(Nm,e,y,d,t-1) =g= -ES_param(e,'RDC') ;
charge_rate_min(Nm,e,y,d,t).. C(Nm,e,y,d,t) - C(Nm,e,y,d,t-1) =l= ES_param(e,'RC') ;
charge_rate_max(Nm,e,y,d,t).. C(Nm,e,y,d,t) - C(Nm,e,y,d,t-1) =g= -ES_param(e,'RC') ;

eq_Psale(Nm,y,d,t).. Psale(Nm,y,d,t) =l= 20;
eq_Ppur(Nm,y,d,t).. Ppur(Nm,y,d,t) =l= 20;
power_balance(Nm,y,d,t).. sum(i,PMG(Nm,i,y,d,t)) + sum(e,Dc(Nm,e,y,d,t) - C(Nm,e,y,d,t)) +
 sum(pv, Cap(pv)*Ps(y,d,t)*SII(Nm,pv,y,d,t))+ Ppur(Nm,y,d,t) =e= Pd(Nm,y,d,t) + Psale(Nm,y,d,t);

Income_eq(Nm,y).. Income(Nm,y) =e= 91.25*sum((d,t),LMP(y,d,t)*Psale(Nm,y,d,t)) ;
OpCost_eq(Nm,y).. OpCost(Nm,y) =e= 91.25*(sum((d,t), sum(i,DG_param(i,'MC')*II(Nm,i,y,d,t)
+ LMP(y,d,t)*Ppur(Nm,y,d,t) + SU(Nm,i,y,d,t) + SD(Nm,i,y,d,t) + DG_param(i,'m')*PMG(Nm,i,y,d,t)))) ;



option resLim = 10000000 ;
option optcr = 0.01 ;