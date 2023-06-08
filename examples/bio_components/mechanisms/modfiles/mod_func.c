#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _Ca_HVA_reg();
extern void _Ca_LVA_reg();
extern void _CaDynamics_reg();
extern void _exp1isyn_reg();
extern void _exp1syn_reg();
extern void _Ih_reg();
extern void _Im_reg();
extern void _Im_v2_reg();
extern void _K_P_reg();
extern void _K_T_reg();
extern void _Kd_reg();
extern void _Kv2like_reg();
extern void _Kv3_1_reg();
extern void _Nap_reg();
extern void _NaTa_reg();
extern void _NaTs_reg();
extern void _NaV_reg();
extern void _SK_reg();
extern void _stp1syn_reg();
extern void _stp2syn_reg();
extern void _stp3syn_reg();
extern void _stp4syn_reg();
extern void _stp5isyn_reg();
extern void _stp5syn_reg();
extern void _vecevent_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," Ca_HVA.mod");
fprintf(stderr," Ca_LVA.mod");
fprintf(stderr," CaDynamics.mod");
fprintf(stderr," exp1isyn.mod");
fprintf(stderr," exp1syn.mod");
fprintf(stderr," Ih.mod");
fprintf(stderr," Im.mod");
fprintf(stderr," Im_v2.mod");
fprintf(stderr," K_P.mod");
fprintf(stderr," K_T.mod");
fprintf(stderr," Kd.mod");
fprintf(stderr," Kv2like.mod");
fprintf(stderr," Kv3_1.mod");
fprintf(stderr," Nap.mod");
fprintf(stderr," NaTa.mod");
fprintf(stderr," NaTs.mod");
fprintf(stderr," NaV.mod");
fprintf(stderr," SK.mod");
fprintf(stderr," stp1syn.mod");
fprintf(stderr," stp2syn.mod");
fprintf(stderr," stp3syn.mod");
fprintf(stderr," stp4syn.mod");
fprintf(stderr," stp5isyn.mod");
fprintf(stderr," stp5syn.mod");
fprintf(stderr," vecevent.mod");
fprintf(stderr, "\n");
    }
_Ca_HVA_reg();
_Ca_LVA_reg();
_CaDynamics_reg();
_exp1isyn_reg();
_exp1syn_reg();
_Ih_reg();
_Im_reg();
_Im_v2_reg();
_K_P_reg();
_K_T_reg();
_Kd_reg();
_Kv2like_reg();
_Kv3_1_reg();
_Nap_reg();
_NaTa_reg();
_NaTs_reg();
_NaV_reg();
_SK_reg();
_stp1syn_reg();
_stp2syn_reg();
_stp3syn_reg();
_stp4syn_reg();
_stp5isyn_reg();
_stp5syn_reg();
_vecevent_reg();
}
