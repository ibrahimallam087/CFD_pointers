/* === Header Files === */
/* These headers provide necessary macros and functions for Fluent UDFs */
#include "udf.h"          // Standard Fluent UDF functions
#include "sg.h"           // Species transport and mass transfer functions
#include "sg_mphase.h"    // Multiphase flow functions (for VOF model)
#include "flow.h"         // Flow-related data (velocity, pressure, etc.)
#include "mem.h"          // Memory allocation functions

/* ==================================================================== */
/* UDF: Compute interfacial area density (interaction between phases)   */
/* This function runs once per iteration and computes the volume fraction */
/* gradient and temperature gradient, then stores their product in UDM[0] */
/* ==================================================================== */
DEFINE_ADJUST(area_density, domain)
{
    /* === Variable Declarations === */
    Thread *t;         // Represents a thread (fluid zone)
    Thread **pt;       // Stores phase-specific threads
    cell_t c;          // Represents a cell in the domain
    Domain *pDomain = DOMAIN_SUB_DOMAIN(domain, P_PHASE);  // Get the primary phase domain
    real voidx, voidy, voidz = 0;  // Stores gradient components

    /* === Compute Volume Fraction Gradient === */
    Alloc_Storage_Vars(pDomain, SV_VOF_RG, SV_VOF_G, SV_NULL); // Allocate memory for gradients
    Scalar_Reconstruction(pDomain, SV_VOF, -1, SV_VOF_RG, NULL); // Compute reconstruction gradient
    Scalar_Derivatives(pDomain, SV_VOF, -1, SV_VOF_G, SV_VOF_RG, Vof_Deriv_Accumulate); // Compute final gradient

    /* === Compute Temperature Gradient === */
    Alloc_Storage_Vars(domain, SV_T_RG, SV_T_G, SV_NULL);
    T_derivatives(domain);  // Compute temperature gradient
    Free_Storage_Vars(domain, SV_T_RG, SV_NULL);  // Free reconstruction gradients

    /* === Loop Over Cells to Compute Interface Area Density === */
    mp_thread_loop_c(t, domain, pt)
    {
        if (FLUID_THREAD_P(t))  // Ensure thread is a fluid region
        {
            Thread *tp = pt[P_PHASE];  // Get the primary phase thread

            begin_c_loop(c, t)  // Loop over all cells
            {
                /* Store the product of VOF gradient and temperature gradient in UDM[0] */
#if RP_3D  
                C_UDMI(c, t, 0) = (C_VOF_G(c, tp)[0] * C_T_G(c, t)[0] +
                                  C_VOF_G(c, tp)[1] * C_T_G(c, t)[1] +
                                  C_VOF_G(c, tp)[2] * C_T_G(c, t)[2]);
#endif  

#if RP_2D  
                C_UDMI(c, t, 0) = (C_VOF_G(c, tp)[0] * C_T_G(c, t)[0] +
                                  C_VOF_G(c, tp)[1] * C_T_G(c, t)[1]);
#endif  
            }
            end_c_loop(c, t)
        }
    }

    /* === Free Memory (UDM still retains values) === */
    Free_Storage_Vars(pDomain, SV_VOF_RG, SV_VOF_G, SV_NULL);
    Free_Storage_Vars(domain, SV_T_G, SV_NULL);
}

/* ==================================================================== */
/* DEFINE_SOURCE: Computes Source Term for Gas Phase                    */
/* This function is used in Fluent's governing equations                 */
/* ==================================================================== */
DEFINE_SOURCE(gas, cell, thread, dS, eqn)
{
    real source;
    Thread *tm = THREAD_SUPER_THREAD(thread); // Get the mixture thread
    Thread **pt = THREAD_SUB_THREADS(tm); // Get sub-threads (phases)

    /* Compute phase conductivity weighted by volume fraction */
    real Kl = C_K_L(cell, pt[1]) * C_VOF(cell, pt[1]);
    real Kg = C_K_L(cell, pt[0]) * C_VOF(cell, pt[0]);
    real L = 1e5;  // Characteristic length scale

    /* Compute source term using stored UDM[0] value */
    source = (Kl + Kg) * C_UDMI(cell, tm, 0) / L;
    C_UDMI(cell, tm, 1) = source;  // Store in UDM[1]
    C_UDMI(cell, tm, 2) = -source * L;  // Store in UDM[2]

    dS[eqn] = 0;  // No derivative term
    return source;
}

/* ==================================================================== */
/* DEFINE_SOURCE: Computes Source Term for Liquid Phase                 */
/* ==================================================================== */
DEFINE_SOURCE(liquid, cell, thread, dS, eqn)
{
    real source;
    Thread *tm = THREAD_SUPER_THREAD(thread);

    source = -C_UDMI(cell, tm, 1);  // Liquid source = - Gas source
    dS[eqn] = 0;
    
    return source;
}

/* ==================================================================== */
/* DEFINE_SOURCE: Computes Source Term for Energy Equation              */
/* ==================================================================== */
DEFINE_SOURCE(energy, cell, thread, dS, eqn)
{
    real source;
    Thread *tm = thread;

    source = C_UDMI(cell, tm, 2);  // Retrieve stored energy source term
    dS[eqn] = 0;
    
    return source;
}

/* ==================================================================== */
/* DEFINE_INIT: Initializes the Volume Fraction Field                   */
/* Used to set initial phase distribution based on a wave function      */
/* ==================================================================== */
DEFINE_INIT(my_init_function, domain)
{
    Thread *t;
    Thread **pt;
    Thread **st;
    cell_t c;
    Domain *pDomain = DOMAIN_SUB_DOMAIN(domain, P_PHASE);  // Get primary phase domain
    Domain *sDomain = DOMAIN_SUB_DOMAIN(domain, S_PHASE);  // Get secondary phase domain
    real xc[ND_ND], y, x;

    /* === Initialize Volume Fraction for Primary Phase === */
    mp_thread_loop_c(t, domain, pt)
    {
        if (FLUID_THREAD_P(t))
        {
            Thread *tp = pt[P_PHASE];

            begin_c_loop(c, t)
            {
                C_CENTROID(xc, c, t); // Get cell centroid
                x = xc[0];
                y = xc[1];

                /* Define a wave-based initial condition for phase distribution */
                if (y < 0.00292 + 0.0006 * cos(6.283 * x / 0.0778))
                    C_VOF(c, tp) = 1;  // Liquid phase
                else
                    C_VOF(c, tp) = 0;  // Gas phase
            }
            end_c_loop(c, t)
        }
    }

    /* === Initialize Volume Fraction for Secondary Phase === */
    mp_thread_loop_c(t, domain, st)
    {
        if (FLUID_THREAD_P(t))
        {
            Thread *sp = st[S_PHASE];

            begin_c_loop(c, t)
            {
                C_CENTROID(xc, c, t);
                x = xc[0];
                y = xc[1];

                /* Define complementary wave-based initial condition */
                if (y < 0.00292 + 0.0006 * cos(6.283 * x / 0.0778))
                    C_VOF(c, sp) = 0;
                else
                    C_VOF(c, sp) = 1;
            }
            end_c_loop(c, t)
        }
    }
}
