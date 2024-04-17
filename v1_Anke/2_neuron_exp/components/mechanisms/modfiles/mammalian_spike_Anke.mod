TITLE HH style channels for spiking retinal ganglion cells
:
: Modified from Coggan et al (2010)


NEURON {
	SUFFIX mammalian_spike_Anke
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnatbar, gnapbar, gkbar 
	RANGE m_inf, h_inf, p_inf, n_inf, c_inf
	RANGE tau_m, tau_h, tau_p, tau_n, tau_c
    RANGE ik, inap, inat
}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	: These conductances are from the axon - they are overwritten in the HOC file for other regions
	gnatbar	= 1.5	(mho/cm2)
	gkbar	= 1.6 (mho/cm2)
	gnapbar = 0.002
	g_pas   = 0.04
	: Equilibrium potentials
	ena 	= 55	(mV)
	ek  	= -77	(mV)
	e_pas   = -70 (mV)
}

STATE {
	m h p n
}

ASSIGNED {
	v (mV)
	celsius (degC)
	inat	(mA/cm2)
	ik	(mA/cm2)
	inap  (mA/cm2)
	m_inf h_inf p_inf n_inf c_inf
	tau_m tau_h tau_p tau_n tau_c
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	inat = gnatbar * m*m*m*h * (v - ena)
	inap = gnapbar * p*p*p*(v - ena)
    idrk = gkbar * n*n*n*(0.9 + 0.1*c) * (v - ek)
}

DERIVATIVE states {	
	trates(v)
	m' = (m_inf - m)/tau_m
	h' = (h_inf - h)/tau_h
	n' = (n_inf - n)/tau_n
	p' = (p_inf - p)/tau_p
	c' = (c_inf - c)/tau_c
}


PROCEDURE trates(vm) { 
	LOCAL a,b
	
:NAT m
	a = (1.76 * (v+21.4)) / (1 - (exp(-(v+21.4)/10.3)))
	b = 0.13 * (-(v+18.7))/(1 - (exp((v + 18.7)/9.16)))
	tau_m = 1 / (a + b)
	m_inf = a * tau_m

:NAT h
	a = 0.062 * (-(v +114))/(1 - (exp((v+114)/11)))
	b = 1.7 / ( 1 + exp((v + 31.8)/13.4))
	tau_h = 1 / (a + b)
	h_inf = a * tau_h

:NAP p
	a = (0.01 * (v+27)) / (1 - (exp(-(v+27)/10.2)))
	b = 0.00025 * (-(v+34))/(1 - (exp((v + 34)/10)))
	tau_p
	p_inf

:K n (non-inactivating, delayed rectifier)
	a = 0.2120*exp(0.04)
	b = 0.1974*exp(0)
	tau_n = 1 / (a + b)
	n_inf = a * tau_n

:K c
	a = 0.00713*exp(-0.1942)
	b = 0.0935*exp(0.0058)
	tau_c = 1/ (a + b)
	c_inf = a *tau_c

}
