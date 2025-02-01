#include "udf.h"

/* Constants for material properties */
#define LAMBDA_L 0.665 /* Thermal conductivity of liquid (W/m.K) */
#define LAMBDA_V 0.025 /* Thermal conductivity of vapor (W/m.K) */
#define HFG 2260000    /* Latent heat of vaporization (J/kg) */
#define T_SAT 373.15   /* Saturation temperature (K) */
#define PI_ 3.14159
/* Nucleation site constants */
#define N0 1000      /* Reference nucleation density (site/m2) */
#define TCRIT 647.15 /* Critical temperature (K) */
#define T0 298.15    /* Reference temperature (K) */
#define GAMMA 0.719  /* Dimensionless constant */
#define A_CONST -0.0002
#define B_CONST 0.122

/* Empirical parameters for nucleation */
#define P_REFR 0.1 /* Reference pressure (MPa) */
#define C_M 2.0    /* Empirical coefficient for mass transfer */
#define RHO_L 997  /* Density of liquid water (kg/m^3) */

/* User-defined function for mass source (evaporation) */
DEFINE_SOURCE(vapor_src, cell, pri_th, dS, eqn)
{
    Thread *mix_th, *sec_th;
    real m_dot = 0.0;

    mix_th = THREAD_SUPER_THREAD(pri_th);
    sec_th = THREAD_SUB_THREAD(mix_th, 1);

    real alpha_v = C_VOF(cell, pri_th);
    real alpha_l = C_VOF(cell, sec_th);
    real T_cell = C_T(cell, mix_th);
    real time = CURRENT_TIME;
    
    if (time > 1.2e-5 * 2 && T_cell > T_SAT)
    {
        real P = C_P(cell, mix_th) / 1e6;
        real delta_Tsup = T_cell - T0;
        real grad_alpha_v = (N0 * 1e-6) / (C_M * RHO_L);
        
        m_dot = ((alpha_v * LAMBDA_V + alpha_l * LAMBDA_L) * grad_alpha_v) / HFG;
        dS[eqn] = 0.0;
    }
    else
    {
        m_dot = 0.0;
        dS[eqn] = 0.0;
    }
    
    return m_dot;
}


/* User-defined function for energy source */
DEFINE_SOURCE(energy_src, cell, pri_th, dS, eqn)
{
    Thread *mix_th, *sec_th;
    real energy_source = 0.0;

    mix_th = THREAD_SUPER_THREAD(pri_th);
    sec_th = THREAD_SUB_THREAD(mix_th, 1);

    real alpha_v = C_VOF(cell, pri_th);
    real alpha_l = C_VOF(cell, sec_th);
    real T_cell = C_T(cell, mix_th);
    real time = CURRENT_TIME;
    
    if (time > 1.2e-5 * 2 && T_cell > T_SAT)
    {
        real P = C_P(cell, mix_th) / 1e6;
        real delta_Tsup = T_cell - T0;
        real grad_alpha_v = (N0 * 1e-6) / (C_M * RHO_L);
        
        energy_source = ((alpha_v * LAMBDA_V + alpha_l * LAMBDA_L) * grad_alpha_v) * HFG / HFG;
        dS[eqn] = 0.0;
    }
    else
    {
        energy_source = 0.0;
        dS[eqn] = 0.0;
    }
    
    return energy_source;
}
