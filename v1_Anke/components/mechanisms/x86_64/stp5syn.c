/* Created by Language version: 7.7.0 */
/* VECTORIZED */
#define NRN_VECTORIZED 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
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
 
#define nrn_init _nrn_init__stp5syn
#define _nrn_initial _nrn_initial__stp5syn
#define nrn_cur _nrn_cur__stp5syn
#define _nrn_current _nrn_current__stp5syn
#define nrn_jacob _nrn_jacob__stp5syn
#define nrn_state _nrn_state__stp5syn
#define _net_receive _net_receive__stp5syn 
#define state state__stp5syn 
 
#define _threadargscomma_ _p, _ppvar, _thread, _nt,
#define _threadargsprotocomma_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt,
#define _threadargs_ _p, _ppvar, _thread, _nt
#define _threadargsproto_ double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 /* Thread safe. No static _p or _ppvar. */
 
#define t _nt->_t
#define dt _nt->_dt
#define e _p[0]
#define tau_1 _p[1]
#define tau_r0 _p[2]
#define tau_FDR _p[3]
#define a_FDR _p[4]
#define tau_D _p[5]
#define a_D _p[6]
#define tau_i _p[7]
#define a_i _p[8]
#define tau_f _p[9]
#define a_f _p[10]
#define pbtilde _p[11]
#define i _p[12]
#define Pmax _p[13]
#define n _p[14]
#define p _p[15]
#define tau_r _p[16]
#define D _p[17]
#define pb _p[18]
#define g _p[19]
#define Dn _p[20]
#define Dp _p[21]
#define Dtau_r _p[22]
#define DD _p[23]
#define Dpb _p[24]
#define Dg _p[25]
#define v _p[26]
#define _g _p[27]
#define _tsav _p[28]
#define _nd_area  *_ppvar[0]._pval
 
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
 /* declaration of user functions */
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

 extern Prop* nrn_point_prop_;
 static int _pointtype;
 static void* _hoc_create_pnt(_ho) Object* _ho; { void* create_point_process();
 return create_point_process(_pointtype, _ho);
}
 static void _hoc_destroy_pnt();
 static double _hoc_loc_pnt(_vptr) void* _vptr; {double loc_point_process();
 return loc_point_process(_pointtype, _vptr);
}
 static double _hoc_has_loc(_vptr) void* _vptr; {double has_loc_point();
 return has_loc_point(_vptr);
}
 static double _hoc_get_loc_pnt(_vptr)void* _vptr; {
 double get_loc_point_process(); return (get_loc_point_process(_vptr));
}
 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _extcall_prop = _prop;
 }
 static void _hoc_setdata(void* _vptr) { Prop* _prop;
 _prop = ((Point_process*)_vptr)->_prop;
   _setdata(_prop);
 }
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 0,0
};
 static Member_func _member_func[] = {
 "loc", _hoc_loc_pnt,
 "has_loc", _hoc_has_loc,
 "get_loc", _hoc_get_loc_pnt,
 0, 0
};
 /* declare global and static user variables */
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "e", "mV",
 "tau_1", "ms",
 "tau_r0", "ms",
 "tau_FDR", "ms",
 "tau_D", "ms",
 "tau_i", "ms",
 "tau_f", "ms",
 "i", "nA",
 0,0
};
 static double D0 = 0;
 static double delta_t = 0.01;
 static double g0 = 0;
 static double n0 = 0;
 static double pb0 = 0;
 static double p0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 static void _hoc_destroy_pnt(_vptr) void* _vptr; {
   destroy_point_process(_vptr);
}
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[2]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"stp5syn",
 "e",
 "tau_1",
 "tau_r0",
 "tau_FDR",
 "a_FDR",
 "tau_D",
 "a_D",
 "tau_i",
 "a_i",
 "tau_f",
 "a_f",
 "pbtilde",
 0,
 "i",
 "Pmax",
 0,
 "n",
 "p",
 "tau_r",
 "D",
 "pb",
 "g",
 0,
 0};
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
  if (nrn_point_prop_) {
	_prop->_alloc_seq = nrn_point_prop_->_alloc_seq;
	_p = nrn_point_prop_->param;
	_ppvar = nrn_point_prop_->dparam;
 }else{
 	_p = nrn_prop_data_alloc(_mechtype, 29, _prop);
 	/*initialize range parameters*/
 	e = 0;
 	tau_1 = 10;
 	tau_r0 = 1000;
 	tau_FDR = 1000;
 	a_FDR = 0.5;
 	tau_D = 100;
 	a_D = 0.5;
 	tau_i = 100;
 	a_i = 0.5;
 	tau_f = 10;
 	a_f = 0.5;
 	pbtilde = 0.5;
  }
 	_prop->param = _p;
 	_prop->param_size = 29;
  if (!nrn_point_prop_) {
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
  }
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _net_receive(Point_process*, double*, double);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _stp5syn_reg() {
	int _vectorized = 1;
  _initlists();
 	_pointtype = point_register_mech(_mechanism,
	 nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init,
	 hoc_nrnpointerindex, 1,
	 _hoc_create_pnt, _hoc_destroy_pnt, _member_func);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 29, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "area");
  hoc_register_dparam_semantics(_mechtype, 1, "pntproc");
  hoc_register_dparam_semantics(_mechtype, 2, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 pnt_receive[_mechtype] = _net_receive;
 pnt_receive_size[_mechtype] = 1;
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 stp5syn /users/students/r0582625/Documents/neural-simulation/v1/components/mechanisms/modfiles/stp5syn.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[6], _dlist1[6];
 static int state(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {int _reset = 0; {
   Dn = ( 1.0 - n ) / tau_r ;
   Dp = ( pb - p ) / tau_f ;
   Dtau_r = ( tau_r0 - tau_r ) / tau_FDR ;
   DD = ( 1.0 - D ) / tau_D ;
   Dpb = ( pbtilde - pb ) / tau_i ;
   Dg = - g / tau_1 ;
   }
 return _reset;
}
 static int _ode_matsol1 (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
 Dn = Dn  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_r )) ;
 Dp = Dp  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_f )) ;
 Dtau_r = Dtau_r  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_FDR )) ;
 DD = DD  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_D )) ;
 Dpb = Dpb  / (1. - dt*( ( ( ( - 1.0 ) ) ) / tau_i )) ;
 Dg = Dg  / (1. - dt*( ( - 1.0 ) / tau_1 )) ;
  return 0;
}
 /*END CVODE*/
 static int state (double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) { {
    n = n + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_r)))*(- ( ( ( 1.0 ) ) / tau_r ) / ( ( ( ( - 1.0 ) ) ) / tau_r ) - n) ;
    p = p + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_f)))*(- ( ( ( pb ) ) / tau_f ) / ( ( ( ( - 1.0 ) ) ) / tau_f ) - p) ;
    tau_r = tau_r + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_FDR)))*(- ( ( ( tau_r0 ) ) / tau_FDR ) / ( ( ( ( - 1.0 ) ) ) / tau_FDR ) - tau_r) ;
    D = D + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_D)))*(- ( ( ( 1.0 ) ) / tau_D ) / ( ( ( ( - 1.0 ) ) ) / tau_D ) - D) ;
    pb = pb + (1. - exp(dt*(( ( ( - 1.0 ) ) ) / tau_i)))*(- ( ( ( pbtilde ) ) / tau_i ) / ( ( ( ( - 1.0 ) ) ) / tau_i ) - pb) ;
    g = g + (1. - exp(dt*(( - 1.0 ) / tau_1)))*(- ( 0.0 ) / ( ( - 1.0 ) / tau_1 ) - g) ;
   }
  return 0;
}
 
static void _net_receive (_pnt, _args, _lflag) Point_process* _pnt; double* _args; double _lflag; 
{  double* _p; Datum* _ppvar; Datum* _thread; _NrnThread* _nt;
   _thread = (Datum*)0; _nt = (_NrnThread*)_pnt->_vnt;   _p = _pnt->_prop->param; _ppvar = _pnt->_prop->dparam;
  if (_tsav > t){ extern char* hoc_object_name(); hoc_execerror(hoc_object_name(_pnt->ob), ":Event arrived out of order. Must call ParallelContext.set_maxstep AFTER assigning minimum NetCon.delay");}
 _tsav = t; {
     if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = g;
    double __primary = (g + _args[0] * Pmax) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( - 1.0 ) / tau_1 ) ) )*( - ( 0.0 ) / ( ( - 1.0 ) / tau_1 ) - __primary );
    g += __primary;
  } else {
 g = g + _args[0] * Pmax ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = n;
    double __primary = (n - n * p) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / tau_r ) ) )*( - ( ( ( 1.0 ) ) / tau_r ) / ( ( ( ( - 1.0 ) ) ) / tau_r ) - __primary );
    n += __primary;
  } else {
 n = n - n * p ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = p;
    double __primary = (p + a_f * ( 1.0 - p )) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / tau_f ) ) )*( - ( ( ( pb ) ) / tau_f ) / ( ( ( ( - 1.0 ) ) ) / tau_f ) - __primary );
    p += __primary;
  } else {
 p = p + a_f * ( 1.0 - p ) ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = tau_r;
    double __primary = (tau_r - a_FDR * tau_r) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / tau_FDR ) ) )*( - ( ( ( tau_r0 ) ) / tau_FDR ) / ( ( ( ( - 1.0 ) ) ) / tau_FDR ) - __primary );
    tau_r += __primary;
  } else {
 tau_r = tau_r - a_FDR * tau_r ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = D;
    double __primary = (D - a_D * p * n * D) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / tau_D ) ) )*( - ( ( ( 1.0 ) ) / tau_D ) / ( ( ( ( - 1.0 ) ) ) / tau_D ) - __primary );
    D += __primary;
  } else {
 D = D - a_D * p * n * D ;
     }
   if (nrn_netrec_state_adjust && !cvode_active_){
    /* discon state adjustment for cnexp case (rate uses no local variable) */
    double __state = pb;
    double __primary = (pb - a_i * pb) - __state;
     __primary += ( 1. - exp( 0.5*dt*( ( ( ( - 1.0 ) ) ) / tau_i ) ) )*( - ( ( ( pbtilde ) ) / tau_i ) / ( ( ( ( - 1.0 ) ) ) / tau_i ) - __primary );
    pb += __primary;
  } else {
 pb = pb - a_i * pb ;
     }
 } }
 
static int _ode_count(int _type){ return 6;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
     _ode_spec1 (_p, _ppvar, _thread, _nt);
 }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
	double* _p; Datum* _ppvar;
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 6; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 (_p, _ppvar, _thread, _nt);
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   double* _p; Datum* _ppvar; Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
 _ode_matsol_instance1(_threadargs_);
 }}

static void initmodel(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt) {
  int _i; double _save;{
  D = D0;
  g = g0;
  n = n0;
  pb = pb0;
  p = p0;
  tau_r = tau_r0;
 {
   n = 1.0 ;
   p = 1.0 ;
   tau_r = tau_r0 ;
   D = 1.0 ;
   pb = pbtilde ;
   g = 0.0 ;
   }
 
}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
double* _p; Datum* _ppvar; Datum* _thread;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
_thread = _ml->_thread;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _tsav = -1e20;
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
 initmodel(_p, _ppvar, _thread, _nt);
}
}

static double _nrn_current(double* _p, Datum* _ppvar, Datum* _thread, _NrnThread* _nt, double _v){double _current=0.;v=_v;{ {
   i = g * ( v - e ) ;
   Pmax = n * p * D ;
   }
 _current += i;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type) {
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
 _g = _nrn_current(_p, _ppvar, _thread, _nt, _v + .001);
 	{ _rhs = _nrn_current(_p, _ppvar, _thread, _nt, _v);
 	}
 _g = (_g - _rhs)/.001;
 _g *=  1.e2/(_nd_area);
 _rhs *= 1.e2/(_nd_area);
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

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type) {
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

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type) {
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
 {   state(_p, _ppvar, _thread, _nt);
  }}}

}

static void terminal(){}

static void _initlists(){
 double _x; double* _p = &_x;
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
 _slist1[1] = &(p) - _p;  _dlist1[1] = &(Dp) - _p;
 _slist1[2] = &(tau_r) - _p;  _dlist1[2] = &(Dtau_r) - _p;
 _slist1[3] = &(D) - _p;  _dlist1[3] = &(DD) - _p;
 _slist1[4] = &(pb) - _p;  _dlist1[4] = &(Dpb) - _p;
 _slist1[5] = &(g) - _p;  _dlist1[5] = &(Dg) - _p;
_first = 0;
}

#if defined(__cplusplus)
} /* extern "C" */
#endif

#if NMODL_TEXT
static const char* nmodl_filename = "/users/students/r0582625/Documents/neural-simulation/v1/components/mechanisms/modfiles/stp5syn.mod";
static const char* nmodl_file_text = 
  "COMMENT\n"
  "A model (# 5 using terminology of Jung Hoon Lee) of a short-term synaptic plasticity.\n"
  "The magnitude of the peak conductance is found by solving Eqs (3,4,5,8,10) in Hennig, 2013. \"Theoretical models of synaptic short term plasticity. Frontiers in computational neuroscience, 7, p.45.\"\n"
  "This file could be used to simulate simpler models by setting to zero the parameters in the unused equations(e.g, a_D,a_i,a_f).\n"
  "\n"
  "State variables:\n"
  "\n"
  "	n = fraction of available vesicles \n"
  "	D = fraction of non-desensitized receptors\n"
  "	tau_r = time constant for vesicle replenishment\n"
  "	pb = baseline vesicle release probability\n"
  "	p = vesicle release probability\n"
  "\n"
  "After each spike the magnitude of the peak conductance changes by the factor w*Pmax, where w is the static synaptic weight and Pmax is the activity-dependent factor that could be interpreted as a probability of a transmitter release by the presynaptic terminal.\n"
  "\n"
  "The post_synaptic dynamics of individual synaptic events is modeled by a single exponential synapse with the time constant tau_1.\n"
  "\n"
  "Implemented by Sergey L. Gratiy,\n"
  "\n"
  "ENDCOMMENT\n"
  "\n"
  ": Declaring parameters as RANGE allows them to change as a function of a position alogn the cell\n"
  "NEURON {\n"
  "	POINT_PROCESS stp5syn\n"
  "	RANGE e, i, tau_1, tau_r0, a_FDR, tau_FDR, a_D, tau_D, a_i, tau_i, a_f, tau_f, pbtilde, Pmax\n"
  "	NONSPECIFIC_CURRENT i\n"
  "}\n"
  "\n"
  "UNITS {\n"
  "	(nA) = (nanoamp)\n"
  "	(mV) = (millivolt)\n"
  "	(umho) = (micromho)\n"
  "}\n"
  "\n"
  ": Declaration of the default values of parameters\n"
  "PARAMETER {\n"
  "	: e = synaptic reversal potential\n"
  "	e = 0 (mV)\n"
  "\n"
  "	: tau_1 = baseline level of release probability\n"
  "	tau_1 = 10 (ms)\n"
  "\n"
  "	: tau_r0 = baseline level of tau_r0\n"
  "	tau_r0 = 1000 (ms)\n"
  "\n"
  "	: tau_FDR = time constant for tau_r relaxation\n"
  "	tau_FDR = 1000 (ms)\n"
  "\n"
  "	: a_FDR = amount of tau_r reduction after each spike\n"
  "	a_FDR = 0.5\n"
  "\n"
  "	: tau_D = relaxation time of D \n"
  "	tau_D = 100 (ms) \n"
  "\n"
  "	: a_D = amount of desentization\n"
  "	a_D = 0.5 \n"
  "\n"
  "	: tau_i = relaxation time for p0 \n"
  "	tau_i = 100 (ms) \n"
  "\n"
  "	: a_i = amount of decrease of baseline probability after each spike \n"
  "	a_i = 0.5 \n"
  "\n"
  "	: tau_f = facilitation time constant (relaxation time constant for p) \n"
  "	tau_f = 10 (ms) \n"
  "\n"
  "	: a_f = amount of facilitation  (increase of p after each spike)\n"
  "	a_f = 0.5 \n"
  "\n"
  "	: pbtilde = baseline level of p0\n"
  "	pbtilde = 0.5 \n"
  "\n"
  "}\n"
  ": Declaration of dependent and external variables that collectively are called ASSIGNED\n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	i (nA)\n"
  "	Pmax\n"
  "}\n"
  "\n"
  ": Declaration of the state variables\n"
  "STATE {\n"
  "	n\n"
  "	p\n"
  "	tau_r\n"
  "	D\n"
  "	pb\n"
  "	g\n"
  "}\n"
  "\n"
  ": Initial conditions for the state variables\n"
  "INITIAL {\n"
  "	n=1\n"
  "	p=1\n"
  "	tau_r=tau_r0\n"
  "	D=1\n"
  "	pb=pbtilde\n"
  "	g=0\n"
  "\n"
  "}\n"
  "\n"
  ": Integration method + assignment statements\n"
  "BREAKPOINT {\n"
  "	SOLVE state METHOD cnexp\n"
  "	i = g*(v - e)\n"
  "	Pmax = n*p*D\n"
  "}\n"
  "\n"
  ": Definition of teh dynamics between the presynaptic activations\n"
  "DERIVATIVE state {\n"
  "	n' = (1-n)/tau_r\n"
  "	p' = (pb-p)/tau_f\n"
  "	tau_r' = (tau_r0-tau_r)/tau_FDR\n"
  "	D' = (1-D)/tau_D\n"
  "	pb' = (pbtilde-pb)/tau_i\n"
  "	g' = -g/tau_1\n"
  "\n"
  "}\n"
  "\n"
  ": This block defines what happens to the state variables at the moment presynaptic activation\n"
  "NET_RECEIVE(weight (umho)) {\n"
  "	g = g + weight*Pmax\n"
  "	n = n - n*p\n"
  "	p = p + a_f*(1-p)\n"
  "	tau_r = tau_r - a_FDR*tau_r\n"
  "	D = D - a_D*p*n*D\n"
  "	pb = pb - a_i*pb\n"
  "\n"
  "}\n"
  "\n"
  "\n"
  ;
#endif
