/*
 * ON Demand User-Defined Functions to check *
 * on the availability of Reconstruction Gradient and Gradients *
 * for a given Solver and Solver settings: *
 * *
 * Availability of Gradients & Reconstruction Gradients depends on: *
 * 1) the selected Solver (density based or pressure based) *
 * 2) the selected Model *
 * 3) the order of discretizations *
 * 4) whether the temporary solver memory is being retained (to keep *
 * temporary memory go to solve - set -expert and type YES *
 * for "Keep temporary solver memory from being freed?") *
 * *
 * *
 * How to use showgrad: *
 * *
 * - Read in your case & data file. *
 * - Compile showgrad.c UDF. *
 * - Load library libudf. *
 * - Attach the showgrad UDF in the Execute on Demand dialog box. *
 * - Run your solution. *
 * - Click the Execute button in the Execute on Demand dialog box. *
 * *
 * A list of available Grads and Recon Grads will be displayed in the *
 * console. *
 * *
 * 2004 Laith Zori *
 */ 
#include "udf.h"
DEFINE_ON_DEMAND(showgrad)
{
    Domain *domain;
    Thread *t;
    domain = Get_Domain(1);
    if (!Data_Valid_P())
        return;
    Message0(" >>> entering show-grad: \n ");
    thread_loop_c(t, domain)
    {
        Material *m = THREAD_MATERIAL(t);
        int nspe = MIXTURE_NSPECIES(m);
        int nspm = nspe - 1;
        Message0("::::\n ");
        Message0(":::: Reconstruction Gradients :::: \n ");
        Message0("::::\n ");
        if (NNULLP(THREAD_STORAGE(t, SV_P_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of P is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_U_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of U is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_V_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of V is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_W_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of W is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_T_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of T is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_H_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of H is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_K_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of K is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_D_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of D is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_O_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of O is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_NUT_RG)))
        {
            Message0("....show-grad:Reconstruction Gradient of NUT is available \n ");
        }
        if (nspe)
        {
            int ns = 0;
            spe_loop(ns, nspm) if (NNULLP(THREAD_STORAGE(t, SV_Y_I(ns) + SV_Y_0_RG - SV_Y_0)))
            {
                Message0("....show-grad:Reconstruction Gradient of Species
                         "available \n ",
                         ns);
            }
        }
        /********************************************************************/
        /********************************************************************/
        /********************************************************************/
        /********************************************************************/
        Message0("::::\n ");
        Message0(":::: Gradients :::: \n ");
        Message0("::::\n ");
        if (NNULLP(THREAD_STORAGE(t, SV_P_G)))
        {
            Message0("....show-grad:Gradient of P is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_U_G)))
        {
            Message0("....show-grad:Gradient of U is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_V_G)))
        {
            Message0("....show-grad:Gradient of V is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_W_G)))
        {
            Message0("....show-grad:Gradient of W is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_T_G)))
        {
            Message0("....show-grad:Gradient of T is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_H_G)))
        {
            Message0("....show-grad:Gradient of H is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_K_G)))
        {
            Message0("....show-grad:Gradient of K is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_D_G)))
        {
            Message0("....show-grad:Gradient of D is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_O_G)))
        {
            Message0("....show-grad:Gradient of O is available \n ");
        }
        if (NNULLP(THREAD_STORAGE(t, SV_NUT_G)))
        {
            Message0("....show-grad:Gradient of NUT is available \n ");
        }
        if (nspe)
        {
            int ns = 0;
spe_loop (ns,nspm)
if (NNULLP(THREAD_STORAGE(t, SV_Y_I(ns)+SV_Y_0_G-SV_Y_0)))
{
Message0("....show-grad:Gradient of Species
}
        }
    }
}