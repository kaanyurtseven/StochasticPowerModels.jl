%% MATPOWER Case Format : Version 2
function mpc = case5_acdc()
mpc.version = '2';


%% stochastic data
%column_names%  dst         pa      pb
mpc.sdata = [
                'Normal'    0.0     0.0;	
                'Beta'		0.4134	1.7147;	
];

%% RES data
%column_names%  RES_bus dst_id 
mpc.RES = [ 2 2;
				3 2;
				4 2;
				5 2;
];


%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;
%% bus data
%    bus_i    type    Pd    Qd    Gs    Bs    area    Vm    Va    baseKV    zone    Vmax    Vmin
mpc.bus = [
	1	3	0.0	0.0	0.0	0.0	1	1.00	0.0	345.0	1	1.1	0.9				
	2	2	20.0	10.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	3	1	45.0	15.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	4	1	40.0	5.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
	5	1	60.0	10.0	0.0	0.0	1	1.0	0.0	345.0	1	1.1	0.9				
];

%% generator data
%    bus    Pg    Qg    Qmax    Qmin    Vg    mBase    status    Pmax    Pmin    Pc1    Pc2    Qc1min    Qc1max    Qc2min    Qc2max    ramp_agc    ramp_10    ramp_30    ramp_q    apf
mpc.gen = [
	1	0.0	0.0	500.0	-500.0	1.06	100.0	1	250.0	00.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
	2	40.0	0.0	300.0	-300.0	1.0	100.0	1	300.0	00.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0				
];

%% branch data
%    f_bus    t_bus    r    x    b    rateA    rateB    rateC    ratio    angle    status    angmin    angmax
mpc.branch = [
	1	2	0.02	0.06	0.06	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	1	3	0.08	0.24	0.05	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	3	0.06	0.18	0.04	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	4	0.06	0.18	0.04	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
	2	5	0.04	0.12	0.03	100.0	100.0	100.0	0	0	1	-59.99999999999999	59.99999999999999								
];

%%-----  OPF Data  -----%%
%% cost data
%    1    startup    shutdown    n    x1    y1    ...    xn    yn
%    2    startup    shutdown    n    c(n-1)    ...    c0
mpc.gencost = [
	2	0.0	0.0	2	100.0	0.0
	2	0.0	0.0	2	300.0	0.0
	%2	 0.0	 0.0	 3	   0.000000	  34.804734	   0.000000; % COW
	%2	 0.0	 0.0	 3	   0.000000	  24.844643	   0.000000; % COW
];

%column_names% μ dst_id λvmax σ λvmin 
mpc.bus_data = {
	0.0	0	1.645	0.0	1.645
	20.0	0	1.645	2.0	1.645
	45.0	1	1.645	4.5	1.645
	40.0	1	1.645	4.0	1.645
	60.0	0	1.645	6.0	1.645
};

%column_names% λqmax λpmax λpmin λqmin 
mpc.gen_data = {
	1.645	1.645	1.645	1.645
	1.645	1.645	1.645	1.645
};

%column_names% λcmax c_rating_a 
mpc.branch_data = {
	%1.645	100.0
	%1.645	100.0
	%1.645	100.0
	%1.645	100.0
	%1.645	100.0
	1.645	100.0
	1.645	100.0
	1.645	100.0
	1.645	100.0
	1.645	100.0
};
