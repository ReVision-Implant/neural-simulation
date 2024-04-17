/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mech_api.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__mammalian_spike_Anke
#define _nrn_initial _nrn_initial__mammalian_spike_Anke
#define nrn_cur _nrn_cur__mammalian_spike_Anke
#define _nrn_current _nrn_current__mammalian_spike_Anke
#define nrn_jacob _nrn_jacob__mammalian_spike_Anke
#define nrn_state _nrn_state__mammalian_spike_Anke
#define _net_receive _net_receive__mammalian_spike_Anke 
#define states states__mammalian_spike_Anke 
#define trates trates__mammalian_spike_Anke 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define gnatbar _p[0]
#define gnatbar_columnindex 0
#define gkbar _p[1]
#define gkbar_columnindex 1
#define gnapbar _p[2]
#define gnapbar_columnindex 2
#define ik _p[3]
#define ik_columnindex 3
#define inat _p[4]
#define inat_columnindex 4
#define inap _p[5]
#define inap_columnindex 5
#define m_inf _p[6]
#define m_inf_columnindex 6
#define h_inf _p[7]
#define h_inf_columnindex 7
#define p_inf _p[8]
#define p_inf_columnindex 8
#define n_inf _p[9]
#define n_inf_columnindex 9
#define c_inf _p[10]
#define c_inf_columnindex 10
#define tau_m _p[11]
#define tau_m_columnindex 11
#define tau_h _p[12]
#define tau_h_columnindex 12
#define tau_p _p[13]
#define tau_p_columnindex 13
#define tau_n _p[14]
#define tau_n_columnindex 14
#define tau_c _p[15]
#define tau_c_columnindex 15
#define m _p[16]
#define m_columnindex 16
#define h _p[17]
#define h_columnindex 17
#define p _p[18]
#define p_columnindex 18
#define n _p[19]
#define n_columnindex 19
#define c _p[20]
#define c_columnindex 20
#define ena _p[21]
#define ena_columnindex 21
#define ek _p[22]
#define ek_columnindex 22
#define Dm _p[23]
#define Dm_columnindex 23
#define Dh _p[24]
#define Dh_columnindex 24
#define Dp _p[25]
#define Dp_columnindex 25
#define Dn _p[26]
#define Dn_columnindex 26
#define Dc _p[27]
#define Dc_columnindex 27
#define ina _p[28]
#define ina_columnindex 28
#define v _p[29]
#define v_columnindex 29
#define _g _p[30]
#define _g_columnindex 30
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 static Datum* _extcall_thread;
 static Prop* _extcall_prop;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_trates(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_mammalian_spike_Anke", _hoc_setdata,
 "trates_mammalian_spike_Anke", _hoc_trates,
 0, 0
};
 /* declare global and static user variables */
#define e_pas e_pas_mammalian_spike_Anke
 double e_pas = -70;
#define g_pas g_pas_mammalian_spike_Anke
 double g_pas = 0.04;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "e_pas_mammalian_spike_Anke", "mV",
 "gnatbar_mammalian_spike_Anke", "mho/cm2",
 "gkbar_mammalian_spike_Anke", "mho/cm2",
 "ik_mammalian_spike_Anke", "mA/cm2",
 "inat_mammalian_spike_Anke", "mA/cm2",
 "inap_mammalian_spike_Anke", "mA/cm2",
 0,0
};
 static double c0 = 0;
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double p0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "g_pas_mammalian_spike_Anke", &g_pas_mammalian_spike_Anke,
 "e_pas_mammalian_spike_Anke", &e_pas_mammalian_spike_Anke,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(NrnThread*, _Memb_list*, int);
static void nrn_state(NrnThread*, _Memb_list*, int);
 static void nrn_cur(NrnThread*, _Memb_list*, int);
static void  nrn_jacob(NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(NrnThread*, _Memb_list*, int);
static void _ode_matsol(NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"mammalian_spike_Anke",
 "gnatbar_mammalian_spike_Anke",
 "gkbar_mammalian_spike_Anke",
 "gnapbar_mammalian_spike_Anke",
 0,
 "ik_mammalian_spike_Anke",
 "inat_mammalian_spike_Anke",
 "inap_mammalian_spike_Anke",
 "m_inf_mammalian_spike_Anke",
 "h_inf_mammalian_spike_Anke",
 "p_inf_mammalian_spike_Anke",
 "n_inf_mammalian_spike_Anke",
 "c_inf_mammalian_spike_Anke",
 "tau_m_mammalian_spike_Anke",
 "tau_h_mammalian_spike_Anke",
 "tau_p_mammalian_spike_Anke",
 "tau_n_mammalian_spike_Anke",
 "tau_c_mammalian_spike_Anke",
 0,
 "m_mammalian_spike_Anke",
 "h_mammalian_spike_Anke",
 "p_mammalian_spike_Anke",
 "n_mammalian_spike_Anke",
 "c_mammalian_spike_Anke",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 31, _prop);
 	/*initialize range parameters*/
 	gnatbar = 1.5;
 	gkbar = 1.6;
 	gnapbar = 0.002;
 	_prop->param = _p;
 	_prop->param_size = 31;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _mammalian_spike_Anke_reg() {
	int _vectorized = 1;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 1);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 31, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 mammalian_spike_Anke mammalian_spike_Anke.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "HH style channels for spiking retinal ganglion cells";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int trates(_threadargsprotocomma_ double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[5], _dlist1[5];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {int _reset = 0; {
   trates ( _threadargscomma_ v ) ;
   Dm = ( m_inf - m ) / tau_m ;
   Dh = ( h_inf - h ) / tau_h ;
   Dn = ( n_inf - n ) / tau_n ;
   Dp = ( p_inf - p ) / tau_p ;
   Dc = ( c_inf - c ) / tau_c ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
 trates ( _threadargscomma_ v ) ;
 Dm = Dm  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_m )) ;
 Dh = Dh  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_h )) ;
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_n )) ;
 Dp = Dp  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_p )) ;
 Dc = Dc  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_c )) ;
  return 0;
}
 /*END CVODE*/
 static int states (double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) { {
   trates ( _threadargscomma_ v ) ;
    m = m + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_m)))*(- ( ( ( m_inf ) ) / tau_m ) / ( ( ( ( - 1.0 ) ) ) / tau_m ) - m) ;
    h = h + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_h)))*(- ( ( ( h_inf ) ) / tau_h ) / ( ( ( ( - 1.0 ) ) ) / tau_h ) - h) ;
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_n)))*(- ( ( ( n_inf ) ) / tau_n ) / ( ( ( ( - 1.0 ) ) ) / tau_n ) - n) ;
    p = p + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_p)))*(- ( ( ( p_inf ) ) / tau_p ) / ( ( ( ( - 1.0 ) ) ) / tau_p ) - p) ;
    c = c + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_c)))*(- ( ( ( c_inf ) ) / tau_c ) / ( ( ( ( - 1.0 ) ) ) / tau_c ) - c) ;
   }
  return 0;
}
 
static int  trates ( _threadargsprotocomma_ double _lvm ) {
   double _la , _lb ;
 _la = ( 1.76 * ( v + 21.4 ) ) / ( 1.0 - ( exp ( - ( v + 21.4 ) / 10.3 ) ) ) ;
   _lb = 0.13 * ( - ( v + 18.7 ) ) / ( 1.0 - ( exp ( ( v + 18.7 ) / 9.16 ) ) ) ;
   tau_m = 1.0 / ( _la + _lb ) ;
   m_inf = _la * tau_m ;
   _la = 0.062 * ( - ( v + 114.0 ) ) / ( 1.0 - ( exp ( ( v + 114.0 ) / 11.0 ) ) ) ;
   _lb = 1.7 / ( 1.0 + exp ( - ( v + 31.8 ) / 13.4 ) ) ;
   tau_h = 1.0 / ( _la + _lb ) ;
   h_inf = _la * tau_h ;
   _la = ( 0.01 * ( v + 27.0 ) ) / ( 1.0 - ( exp ( - ( v + 27.0 ) / 10.2 ) ) ) ;
   _lb = 0.00025 * ( - ( v + 34.0 ) ) / ( 1.0 - ( exp ( ( v + 34.0 ) / 10.0 ) ) ) ;
   tau_p = 1.0 / ( _la + _lb ) ;
   p_inf = _la * tau_p ;
   _la = 0.2120 * exp ( 0.04 * v ) ;
   _lb = 0.1974 * exp ( 0.0 * v ) ;
   tau_n = 1.0 / ( _la + _lb ) ;
   n_inf = _la * tau_n ;
   _la = 0.00713 * exp ( - 0.1942 * v ) ;
   _lb = 0.0935 * exp ( 0.0058 * v ) ;
   tau_c = 1.0 / ( _la + _lb ) ;
   c_inf = _la * tau_c ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
   double* _p; Datum* _ppvar; Datum* _thread; NrnThread* _nt;
   if (_extcall_prop) {_p = _extcall_prop->param; _ppvar = _extcall_prop->dparam;}else{ _p = (double*)0; _ppvar = (Datum*)0; }
  _thread = _extcall_thread;
  _nt = nrn_threads;
 _r = 1.;
 trates ( _p, _ppvar, _thread, _nt, *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 5;}
 
static void _ode_spec(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 (_p, _ppvar, _thread, _nt);
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 5; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt) {
  int _i; double _save;{
  c = c0;
  h = h0;
  m = m0;
  n = n0;
  p = p0;
 
}
}

static void nrn_init(NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel(_p, _ppvar, _thread, _nt);
  }
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   ik = gkbar * n * n * n * ( 0.9 + 0.1 * c ) * ( v - ek ) ;
   inat = gnatbar * m * m * m * h * ( v - ena ) ;
   inap = gnapbar * p * p * p * ( v - ena ) ;
   ina = inat + inap ;
   }
 _current += ina;
 _current += ik;

} return _current;
}

static void nrn_cur(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}
 
}

static void nrn_jacob(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}
 
}

static void nrn_state(NrnThread* _nt, _Memb_list* _ml, int _type) {
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
 {   states(_p, _ppvar, _thread, _nt);
  }  }}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = m_columnindex;  _dlist1[0] = Dm_columnindex;
 _slist1[1] = h_columnindex;  _dlist1[1] = Dh_columnindex;
 _slist1[2] = n_columnindex;  _dlist1[2] = Dn_columnindex;
 _slist1[3] = p_columnindex;  _dlist1[3] = Dp_columnindex;
 _slist1[4] = c_columnindex;  _dlist1[4] = Dc_columnindex;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "mammalian_spike_Anke.mod";
static const char* nmodl_file_text = 
  "TITLE HH style channels for spiking retinal ganglion cells\n"
  ":\n"
  ": Modified from Coggan et al (2010)\n"
  "\n"
  "\n"
  "NEURON {\n"
  "	SUFFIX mammalian_spike_Anke\n"
  "	USEION na READ ena WRITE ina\n"
  "	USEION k READ ek WRITE ik\n"
  "	RANGE gnatbar, gnapbar, gkbar \n"
  "	RANGE m_inf, h_inf, p_inf, n_inf, c_inf\n"
  "	RANGE tau_m, tau_h, tau_p, tau_n, tau_c\n"
  "    RANGE ik, inap, inat\n"
  "}\n"
  "\n"
  "\n"
  "UNITS {\n"
  "	(molar) = (1/liter)\n"
  "	(mM) = (millimolar)\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "}\n"
  "\n"
  "PARAMETER {\n"
  "	: These conductances are from the axon - they are overwritten in the HOC file for other regions\n"
  "	gnatbar	= 1.5	(mho/cm2)\n"
  "	gkbar	= 1.6 (mho/cm2)\n"
  "	gnapbar = 0.002\n"
  "	g_pas   = 0.04\n"
  "	: Equilibrium potentials\n"
  "	ena 	= 55	(mV)\n"
  "	ek  	= -77	(mV)\n"
  "	e_pas   = -70 (mV)\n"
  "}\n"
  "\n"
  "STATE {\n"
  "	m h p n c\n"
  "}\n"
  "\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC)\n"
  "	ik	(mA/cm2)\n"
  "	ina (mA/cm2)\n"
  "		inat	(mA/cm2)\n"
  "		inap  (mA/cm2)\n"
  "	m_inf h_inf p_inf n_inf c_inf\n"
  "	tau_m tau_h tau_p tau_n tau_c\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "    ik = gkbar * n*n*n*(0.9 + 0.1*c) * (v - ek)\n"
  "		inat = gnatbar * m*m*m*h * (v - ena)\n"
  "		inap = gnapbar * p*p*p*(v - ena)\n"
  "		ina = inat + inap\n"
  "}\n"
  "\n"
  "DERIVATIVE states {	\n"
  "	trates(v)\n"
  "	m' = (m_inf - m)/tau_m\n"
  "	h' = (h_inf - h)/tau_h\n"
  "	n' = (n_inf - n)/tau_n\n"
  "	p' = (p_inf - p)/tau_p\n"
  "	c' = (c_inf - c)/tau_c\n"
  "}\n"
  "\n"
  "\n"
  "PROCEDURE trates(vm) { \n"
  "	LOCAL a,b\n"
  "	\n"
  ":NAT m\n"
  "	a = (1.76 * (v+21.4)) / (1 - (exp(-(v+21.4)/10.3)))\n"
  "	b = 0.13 * (-(v+18.7))/(1 - (exp((v + 18.7)/9.16)))\n"
  "	tau_m = 1 / (a + b)\n"
  "	m_inf = a * tau_m\n"
  "\n"
  ":NAT h\n"
  "	a = 0.062 * (-(v +114))/(1 - (exp((v+114)/11)))\n"
  "	b = 1.7 / ( 1 + exp(-(v + 31.8)/13.4))\n"
  "	tau_h = 1 / (a + b)\n"
  "	h_inf = a * tau_h\n"
  "\n"
  ":NAP p\n"
  "	a = (0.01 * (v+27)) / (1 - (exp(-(v+27)/10.2)))\n"
  "	b = 0.00025 * (-(v+34))/(1 - (exp((v + 34)/10)))\n"
  "	tau_p = 1/ (a+b)\n"
  "	p_inf = a*tau_p\n"
  "\n"
  ":K n (non-inactivating, delayed rectifier)\n"
  "	a = 0.2120*exp(0.04*v)\n"
  "	b = 0.1974*exp(0*v)\n"
  "	tau_n = 1 / (a + b)\n"
  "	n_inf = a * tau_n\n"
  "\n"
  ":K c\n"
  "	a = 0.00713*exp(-0.1942*v)\n"
  "	b = 0.0935*exp(0.0058*v)\n"
  "	tau_c = 1/ (a + b)\n"
  "	c_inf = a *tau_c\n"
  "\n"
  "}\n"
  ;
#endif
