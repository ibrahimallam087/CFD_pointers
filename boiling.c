
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
static real total_mdot = 0;
static real cell_count =0;

// mass source term for vapor phase 
DEFINE_SOURCE(vapor_src, cell, pri_th, dS, eqn)
{
    Thread *mix_th, *sec_th;
    real m_dot = 0.0;

    /* Thread assignments */
    mix_th = THREAD_SUPER_THREAD(pri_th);
    sec_th = THREAD_SUB_THREAD(mix_th, 1);

    real alpha_v = C_VOF(cell, pri_th); /* Volume fraction of vapor */
    real alpha_l = C_VOF(cell, sec_th); /* Volume fraction of liquid */
    
    real p_op = RP_Get_Real("operating-pressure");
    /* Get temperature */
    real T_cell = C_T(cell, mix_th);
    real time = CURRENT_TIME;
    
    if (T_cell > T_SAT)
    {
       

        /* Calculate Nucleation Density */
        //real P = C_P(cell, mix_th) / 1e6;
        real P = (C_P(cell, mix_th)+p_op); /* Convert pressure to MPa */
        real delta_Tsup = T_cell - T_SAT;    /* Superheat relative to reference */
        real f_P = 26.006 - 3.678 * exp(-2 * P) - 21.907 * exp(-P / 24.065);
        real A = A_CONST * P * P + 0.0108 * P + 0.0119;
        real B = B_CONST * P + 1.988;
        real cos_theta = (1 - cos(41.37 * PI_ / 180)) * pow((TCRIT - T_cell) / (TCRIT - T0), GAMMA);
        real Nw = N0 * cos_theta * exp(f_P) * pow(delta_Tsup, delta_Tsup * A + B);
        
        /* Compute volume fraction gradient using nucleation density */
        real grad_alpha_v = (Nw)/ (C_M * RHO_L); /* 10 scaling for proper unit conversion */
        
        // printing some vars 
        // Message(" L VOF : %g ------- V VOF : %g \n",alpha_l,alpha_v);
        // Message(" Cell Temperature : %g \n",T_cell);
        // Message(" Cell Pressure : %g \n",P);
        // Message(" NW : %g \n",Nw);
        // Message(" DOT Product term : %g \n",grad_alpha_v);
        // -------------------------------------------------//

        

        
        /* Compute volumetric mass source */
        m_dot = ((alpha_v * LAMBDA_V + alpha_l * LAMBDA_L) * grad_alpha_v) / HFG;
        total_mdot+=m_dot;
        cell_count++;

        //Message("mass transferd from Liquid to vapor in the current cell %g\n\n\n",m_dot);
     
        dS[eqn] = 0; /* Explicit source */
    }
    else
    {
        m_dot = 0.0;
        dS[eqn] = 0.0;
    }
    
    return m_dot;
}

// mass source term for liquid phase 
DEFINE_SOURCE(liq_src, cell, sec_th, dS, eqn)
{
    Thread *mix_th, *pri_th;
    real m_dot = 0.0;

    /* Thread assignments */
    mix_th = THREAD_SUPER_THREAD(sec_th);
    pri_th = THREAD_SUB_THREAD(mix_th, 0);

    real alpha_v = C_VOF(cell, pri_th); /* Volume fraction of vapor */
    real alpha_l = C_VOF(cell, sec_th); /* Volume fraction of liquid */

    /* Get temperature */
    real T_cell = C_T(cell, mix_th);
    real time = CURRENT_TIME;
    real p_op = RP_Get_Real("operating-pressure");

    if (T_cell > T_SAT)
    {
        /* Calculate Nucleation Density */
        real P = (C_P(cell, mix_th)+p_op); /* Convert pressure to MPa */
        real delta_Tsup = T_cell - T_SAT;    /* Superheat relative to reference */
        real f_P = 26.006 - 3.678 * exp(-2 * P) - 21.907 * exp(-P / 24.065);
        real A = A_CONST * P * P + 0.0108 * P + 0.0119;
        real B = B_CONST * P + 1.988;
        real cos_theta = (1 - cos(41.37 * PI_ / 180)) * pow((TCRIT - T_cell) / (TCRIT - T0), GAMMA);
        real Nw = N0 * cos_theta * exp(f_P) * pow(delta_Tsup, delta_Tsup * A + B);
        
        /* Compute volume fraction gradient using nucleation density */
        real grad_alpha_v = (Nw)/ (C_M * RHO_L); /* 1 scaling for proper unit conversion */
        
        /* Compute volumetric mass source */
        m_dot = -((alpha_v * LAMBDA_V + alpha_l * LAMBDA_L) * grad_alpha_v) / HFG;
        dS[eqn] = 0; /* Explicit source */
   
    }
    else
    {
        m_dot = 0.0;
        dS[eqn] = 0.0;
    }
    
    return m_dot;
}

// Energy source term
DEFINE_SOURCE(enrg_src,cell,mix_th,ds,eqn)
{
    Thread *pri_th,*sec_th;
    real m_dot;
    pri_th = THREAD_SUB_THREAD(mix_th,0);
    sec_th = THREAD_SUB_THREAD(mix_th,1);
    real time = CURRENT_TIME;
    real T_cell = C_T(cell, mix_th);
    real alpha_v = C_VOF(cell,pri_th);
    real alpha_l = C_VOF(cell,sec_th);
    real p_op = RP_Get_Real("operating-pressure");
    if( T_cell>T_SAT){
         /* Calculate Nucleation Density */
        real P = (C_P(cell, mix_th)+p_op); /* Convert pressure to MPa */
        real delta_Tsup = T_cell - T_SAT;    /* Superheat relative to reference */
        real f_P = 26.006 - 3.678 * exp(-2 * P) - 21.907 * exp(-P / 24.065);
        real A = A_CONST * P * P + 0.0108 * P + 0.0119;
        real B = B_CONST * P + 1.988;
        real cos_theta = (1 - cos(41.37 * PI_ / 180)) * pow((TCRIT - T_cell) / (TCRIT - T0), GAMMA);
        real Nw = N0 * cos_theta * exp(f_P) * pow(delta_Tsup, delta_Tsup * A + B);
        
        /* Compute volume fraction gradient using nucleation density */
        real grad_alpha_v = (Nw)/ (C_M * RHO_L); /* 1 scaling for proper unit conversion */
        
        /* Compute volumetric mass source */
        m_dot = ((alpha_v * LAMBDA_V + alpha_l * LAMBDA_L) * grad_alpha_v) / HFG;
        ds[eqn] = 0; /* Explicit source */
       


    }
    else{
        ds[eqn] = 0;
        m_dot = 0;
    }
    return HFG*m_dot;


}
DEFINE_EXECUTE_AT_END(report_func){
    if(cell_count>0)
    {
        real avg_mdot =  total_mdot/cell_count;
        Message(" Average mass transferd to vapor is %d\n",avg_mdot);
        Message(" cell count %d\n",cell_count);
        cell_count = 0;
        total_mdot = 0;
        
    }
    else 
    {
        Message("Thers is no evaporation in this iteration \n");

    }
}
