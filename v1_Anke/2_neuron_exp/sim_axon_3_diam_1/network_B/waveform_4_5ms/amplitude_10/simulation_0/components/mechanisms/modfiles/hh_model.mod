TITLE HH model test as general axon model


NEURON {
	SUFFIX hh_model
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gnabar, gkbar
	RANGE m_inf, h_inf, n_inf
	RANGE tau_m, tau_h, tau_n
    RANGE idrk, ina
}


UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	: These conductances are from the axon
	gnabar	= 0.040	(mho/cm2)
	gkbar	= 0.035 (mho/cm2)
	g_pas   = 0.003
	: Equilibrium potentials
	ena 	= 55	(mV)
	ek  	= -77	(mV)
	e_pas   = -65 (mV)
}

STATE {
	m h n c 
}

INITIAL {
: The initial values were determined at a resting value of -65.02 mV (unchanged compared to 5 channel model)
    m = 0.0353
    h = 0.9054
    n = 0.0677
}

ASSIGNED {
	v (mV)
	celsius (degC)
	ina	(mA/cm2)
	ik	(mA/cm2)
	m_inf h_inf n_inf
	tau_m tau_h tau_n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m*m*m*h * (v - ena)
    ik = gkbar * n*n*n*n * (v - ek)
}

DERIVATIVE states {	
	trates(v)
	m' = (m_inf - m)/tau_m
	h' = (h_inf - h)/tau_h
	n' = (n_inf - n)/tau_n
}


PROCEDURE trates(vm) { 
	LOCAL a,b
	
:NA m
	a = (0.182* (v+35)) / (1 - (exp(-1*(v+35)/9)))
	b = -0.124 * (v+35) / (1 - (exp((v+35)/9)))
	tau_m = 1 / (a + b)
	m_inf = a * tau_m

:NA h
	a = 0.25 * (exp((-1*(v+90))/12))
	b = 0.25 * (exp(((v+62))/6))/(exp((v+90)/12))
	tau_h = 1 / (a + b)
	h_inf = a * tau_h

:K n (non-inactivating, delayed rectifier)
	a = (0.02* (v-25)) / (1 - (exp(-1*(v-25)/9)))
	b = -0.002 * (v-25) / (1 - (exp((v-25)/9)))
	tau_n = 1 / (a + b)
	n_inf = a * tau_n
}
