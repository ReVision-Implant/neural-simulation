#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _mammalian_spike_Anke_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"mammalian_spike_Anke.mod\"");
    fprintf(stderr, "\n");
  }
  _mammalian_spike_Anke_reg();
}
