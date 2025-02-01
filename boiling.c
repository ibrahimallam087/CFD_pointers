    #include "udf.h"  
    #include "sg.h"  
    #include "sg_mphase.h"  
    #include "flow.h"  
    #include "mem.h"  
      
      
      
    /**************************************************************/  
    /* UDF for specifying an interfacail area density             */  
    /**************************************************************/  
      
      
    DEFINE_ADJUST(area_density, domain)  
    {  
      Thread *t;  
      Thread **pt;  
      cell_t c;  
      Domain *pDomain = DOMAIN_SUB_DOMAIN(domain,P_PHASE);  
      real voidx, voidy, voidz=0;  
        
         
      {  
      Alloc_Storage_Vars(pDomain,SV_VOF_RG,SV_VOF_G,SV_NULL);  
      Scalar_Reconstruction(pDomain, SV_VOF,-1,SV_VOF_RG,NULL);  
      Scalar_Derivatives(pDomain,SV_VOF,-1,SV_VOF_G,SV_VOF_RG,  
                 Vof_Deriv_Accumulate);  
       }  
                   
       {  
          Alloc_Storage_Vars(domain, SV_T_RG, SV_T_G,  SV_NULL);  
          T_derivatives(domain);  
          Free_Storage_Vars(domain, SV_T_RG, SV_NULL);  
        }  
                   
          mp_thread_loop_c (t,domain,pt)  
        if (FLUID_THREAD_P(t))  
          {  
            Thread *tp = pt[P_PHASE];  
      
            begin_c_loop (c,t)  
              {  
      
    #if RP_3D  
        C_UDMI(c,t,0) = (C_VOF_G(c,tp)[0]*C_T_G(c,t)[0]+  
        C_VOF_G(c,tp)[1]*C_T_G(c,t)[1]+C_VOF_G(c,tp)[2]*C_T_G(c,t)[2]);  
    #endif  
      
    #if RP_2D  
        C_UDMI(c,t,0) = (C_VOF_G(c,tp)[0]*C_T_G(c,t)[0]+  
        C_VOF_G(c,tp)[1]*C_T_G(c,t)[1]);  
    #endif  
      
              }  
              end_c_loop (c,t)  
             }  
               
     Free_Storage_Vars(pDomain,SV_VOF_RG,SV_VOF_G,SV_NULL);  
     Free_Storage_Vars(domain, SV_T_G, SV_NULL);  
    }  
      
    DEFINE_SOURCE(gas, cell, thread, dS, eqn)  
    {  
      
      real x[ND_ND];  
      real source;  
      Thread *tm = THREAD_SUPER_THREAD(thread);  
      Thread **pt = THREAD_SUB_THREADS(tm);  
         
      real Kl = C_K_L(cell, pt[1])*C_VOF(cell, pt[1]),  
      Kg = C_K_L(cell, pt[0])*C_VOF(cell, pt[0]);  
      real L = 1e5;  
      
      source = (Kl+Kg)*C_UDMI(cell,tm,0) / L;  
        
       
        
      C_UDMI(cell, tm, 1) = source;  
        
      
      C_UDMI(cell, tm, 2) = -source*L;  
      
        
        
      dS[eqn] =0;  
        
      return source;  
    }  
      
    DEFINE_SOURCE(liquid, cell, thread, dS, eqn)  
    {  
      real x[ND_ND];  
      real source;  
      Thread *tm = THREAD_SUPER_THREAD(thread);  
      Thread **pt = THREAD_SUB_THREADS(tm);  
      
      source = -C_UDMI(cell, tm, 1);  
        
      dS[eqn] = 0;  
        
      return source;  
    }  
      
    DEFINE_SOURCE(energy, cell, thread, dS, eqn)  
    {  
      real x[ND_ND];  
      real source;  
      Thread *tm = thread;  
        
        
      
      source = C_UDMI(cell, tm, 2);  
      dS[eqn] = 0;  
        
      return source;  
    }  
      
    /***********************************************************************/  
    /* UDF for initializing flow field variables                           */  
    /***********************************************************************/  
      
      
    DEFINE_INIT(my_init_function, domain)  
    {  
      Thread *t;  
      Thread **pt;  
      Thread **st;  
      cell_t c;  
      Domain *pDomain = DOMAIN_SUB_DOMAIN(domain,P_PHASE);  
      Domain *sDomain = DOMAIN_SUB_DOMAIN(domain,S_PHASE);  
        
      real xc[ND_ND], y, x;  
      
          mp_thread_loop_c (t,domain,pt)  
        if (FLUID_THREAD_P(t))  
          {  
            Thread *tp = pt[P_PHASE];  
      
            begin_c_loop (c,t)  
              {  
         C_CENTROID(xc,c,t);  
         x=xc[0];  
         y=xc[1];  
           
          if ( y < 0.00292 + 0.0006*cos(6.283*x/0.0778) )  
          C_VOF(c,tp) = 1;  
          else  
          C_VOF(c,tp) = 0;  
      
              }  
              end_c_loop (c,t)  
             }  
               
                   mp_thread_loop_c (t,domain,st)  
        if (FLUID_THREAD_P(t))  
          {  
            Thread *sp = st[S_PHASE];  
      
            begin_c_loop (c,t)  
              {  
         C_CENTROID(xc,c,t);  
         x=xc[0];  
         y=xc[1];  
           
          if ( y < 0.00292 + 0.0006*cos(6.283*x/0.0778) )  
          C_VOF(c,sp) = 0;  
          else  
          C_VOF(c,sp) = 1;  
      
              }  
              end_c_loop (c,t)  
             }  
      
    }  

