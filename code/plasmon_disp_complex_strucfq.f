       Program StructureFactor with Q  
       implicit double precision (a-h,o-z)                                                                                                                                                                                                  
       complex(8):: drl,fl,fli1,flki,u,xch,fqw,fxccp07,const,one                                                                                                                                                                                  
        parameter( pi = 3.1415926535897932384626433832795d0)
        one = cmplx(1.d0,0.d0)                                                                                                                                                              
        sig = 0.0d0                                                                                                                                                                                                                                                                                                                                                                                                         
        dW = 0.005d0 
        dQ = 0.05d0   
        dlam = 0.01d0                                                                                                                                                                                                   
        ctil = 0.0037d0                                                                                                                                                                                                     
c       m = 1                                                                                                                                                                                                             
c       hbar = 1                                                                                                                                                                                                          
                                                                                                                                                                                                                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                                                                      
c The factors of NEO kernel                                                                                                                                                                                               
       facp = ((1.d0 + sig)/2.d0)**2                                                                                                                                                                                      
       facn = ((1.d0 - sig)/2.d0)**2                                                                                                                                                                                      
       szp = (1.d0 + sig)**(2.d0/3.d0)                                                                                                                                                                                    
       szn = (1.d0 - sig)**(2.d0/3.d0)   
       facpx = (1.d0 + sig)**(4.d0/3.d0)                           
       facnx = (1.d0 - sig)**(4.d0/3.d0)                           
       facx = (facpx + facnx)/2.d0                                 
        b = (3.d0/(4.d0*pi))**(1.d0/3.d0)                          
       rs = 5.62d0        
       dens = 3.d0/(4.d0*pi*rs**3.d0)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
       akf = ((9.d0*pi/4.d0)**(1.d0/3.d0))/rs                                                                                                                                                                             
       ef = akf**2/2.d0                                                                                                                                                                                                   
       wp = dsqrt(3.d0/rs**3) 
       constr = -(3.d0*wp**2/(4.d0*pi*akf**2))       
       const = cmplx(constr,0.d0)                                                                                                                                                                                                 
                                                                                                                                                                                                                          
c Calculate the structure factor for Q, integrating the structure factor for W and lambda                                                                                                              
c We want to reproduce Figure 2 of Langreth and Perdew, Phys. Rev. B, 15, 2884 (1977)


         Wr =0.d0
         do iQ = 1, 50             
           Q = dQ*iQ                       
c           if (iQ.eq.1) omegar = wp/ef
c           if (iQ.eq.1) omegai = 0.0d0   

          
           sumw = 0.d0 
           sumwx = 0.0d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
           do jW = 1, 120001                                                                                                                                                                                                 
             Wi = dW*(jW - 1.d0)
             omegai = Wi*ef
             
             do mlam = 1, 101
                alam = dlam*(mlam -1.d0)
             
c The variables of the complex Lindhard function                             
c          z = Q/2.d0 
          ur = Wr/(4.d0*Q)                                                       
          ui = Wi/(4.d0*Q)   
          u = cmplx(ur,ui-1.d-9)                                                                                             
          Call f(Q,u,akf,rs,fl)  
          fl = const*fl      
                                                       
c         fxc = 0.d0                                                                                                                                                                                                   
c The NEO kernel      
c Density parameters                                                                          
        d2exc = ALDAxc(Q,ur,akf,rs,sq,sig)               
        cc = (3.d0*(pi**2)*((b/rs)**3.d0))**(-2.d0/3.d0) 
        clow = (-pi/2.d0)*cc*facx/d2exc                  
         z = Q                                                                                                                                                                                           
         zp = (z**2)/((clow)*szp)                                                                                                                                                                                         
         zn = (z**2)/((clow)*szn)                                                                                                                                                                                         
         fzp = 1.d0 - dexp(-zp)                                                                                                                                                                                           
         fzn = 1.d0 - dexp(-zn)                                                                                                                                                                                           
         sq = 2.d0*Q*akf                                                                                                                                                                                                       
         vcr = alam*4.d0*pi/sq**2                                                                                                                                                                                            
          cl = -vcr    
         fxc = cl*(facp*fzp + facn*fzn)                                                                                                                                                                                                 
c         fxcalx = ALDAx(Q,u,akf,rs,sq,sig) 
c          fxcalx = cmplx(fxcalx, 0.d0)
c         fxcalc = ALDAxc(Q,u,akf,rs,sq,sig)
c          fxcalc = cmplx(fxcalc, 0.d0)
c         fxcralx = rALDAx(Q,u,akf,rs,sq,sig) 
c          fxcralx = cmplx(fxcralx, 0.d0)
c         fxcralc = rALDAxc(Q,u,akf,rs,sq,sig)
c          fxcralc = cmplx(fxcralc, 0.d0)
c           Call CP07(Q,u,akf,rs,sq,sig,fxccp07) 
c        fxc = fxcalx                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
c Construct the Hxc kernel and the structure factor                                                                                                                                                
         xchr = vcr + fxc 
         xch = cmplx(xchr,0.d0)                                                                                                                                                                                                  
         fli1  = one/(one - fl*xch)
         flki = fli1*fl
         fli = real(flki)
         sf = -(1.d0/pi)*fli/dens
         sfx = -(1.d0/pi)*fl/dens        
         sumw = sumw + dlam*dW*ef*sf  
         sumwx = sumwx + dlam*dW*ef*sfx
         
c 300     write(*, 4) Wi, fl   
c 4       format(3f15.6)
      
          enddo
         enddo
        sumwm = sumw - 1.d0
        sumwxm = sumwx - 1.d0 
        sumwc = sumwm - sumwxm
200     write(*, 3) Q, sumwm, sumwxm, sumwc       
3       format(4f30.6)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
100      enddo                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                    
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                              
        end                                                                                                                                                                                                                                                                                                                                                                                                                                        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                                                                          
         Subroutine f(z,uu,akf,rs,fl) 
         implicit double precision (a-h,o-z)                                                                                                                                                         
         complex(8)::fl,fx,uu,fx1,fx2,fx3,fx4,zu1,zu2,one,half,eight                                                                                                                                                                 
          complex(8):: zz                                                                                                                                                                                      
         parameter( pi = 3.1415926535897932384626433832795d0)                                                                                                                                                   
                                                                                                                                                                                                               
c The complex frequency-dependent Lindhard function from J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954)                                                                                                                                
                                                                                             
         one = cmplx(1.d0, 0.d0)                                                                                                                                                                                                           
         half = cmplx(1.d0/2.d0, 0.d0)                                                                                                                                                                           
         eight = cmplx(8.d0, 0.d0)                                                                                                                                                                               
         zz = cmplx(z, 0.d0)    
 
                                                                                                                                                                                     
         zu1 = zz - uu                                                                                                                                                                                           
         zu2 = zz + uu                                                                                                                                                                                           
                                                                                                                                                                                                                 
         fx1 = ((one/(eight*zz))*(one -zu1**2))                                                                                                                                                                  
         fx2 = cdlog((zu1 + one)/(zu1 - one))                                                                                                                                                    
         fx3 = ((one/(eight*zz))*(one - zu2**2))                                                                                                                                                          
         fx4 = cdlog((zu2 + one)/(zu2 - one))                                                                                                                                                               
         fx = half + fx1*fx2 + fx3*fx4                                                                                                                                                                                                                                                                               
         fl = fx                                                                                                                                                                                              
         
         return                                                                                                                                                                                                                                                                                                                                                                                                                          
         end    
               
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         Double Precision Function ALDAx(Q,u,akf,rs,sq,sig)        
         implicit double precision (a-h,o-z)                                          
         parameter( pi = 3.1415926535897932384626433832795d0)                         
                                                                                      
         aac = 1.d0/4.d0                                                              
         akc = (4.d0*pi)/(akf**2)                                                                                                                                                                                                        
         alxc = akc*aac                                                                                                                                
         ALDAx = -alxc                                               
         return                                                                       
         end                                                                          
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         Double Precision Function ALDAxc(Q,u,akf,rs,sq,sig)      
         implicit double precision (a-h,o-z)                                        
         parameter( pi = 3.1415926535897932384626433832795d0)                       
                                                                                    
         aac = 1.d0/4.d0                                                             
         akc = akf**2/(4.d0*pi)
c        ecu = gu/(1.d0 + b1u*dsqrt(1.d0/rs) + b2u*(1.d0/rs))                   
c        ecp = gp/(1.d0 + b1p*dsqrt(1.d0/rs) + b2p*(1.d0/rs)                   
c        enfsig1 = (1.d0 + sig)**(4.d0/3.d0)                        
c        enfsig2 =(1.d0 - sig)**(4.d0/3.d0)-2.d0                    
c        fsig = (2.d0**(4.d0/3.d0))-2.d0                            
c The uniform electron gas correlation energy according to Perdew and Zunger, Phys. Rev. B, 23, 5076 (1981)               
         gu =-0.1423d0                                                                                                    
         b1u = 1.0529d0                                                                                                   
         b2u = 0.3334d0                                                                                                   
         cu = 0.0020d0                                                                                                    
         du = -0.0116d0                                                                                                   
         gp = -0.0843d0                                                                                                   
         b1p = 1.3981d0                                                                                                   
         b2p = 0.2611d0                                                                                                   
         cp =  0.0007d0                                                                                                   
         dp = -0.0048d0   
         b = (3.d0/(4.d0*pi))**(1.d0/3.d0)  
         sqrs = dsqrt(rs)                                            
         AA = b*b2u*(21.d0*b1u +16.d0*b2u*sqrs)                      
         BB = 5.d0*b1u+7.d0*(b1u**2)*sqrs+8.d0*b2u*sqrs              
         decdrsn = gu*b**2*(AA + BB*b/rs)                            
         decdrsd = 36.d0*(b*b1u+(b*b2u+b/rs)*sqrs)**3*(b/rs)**3                                                                                                                                                                           
         decdrs = decdrsn/decdrsd
         alrs = akc*decdrs
         aa = aac - alrs                                                       
         ALDAxc = -(1.d0/akc)*aa
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
         Double Precision Function rALDAx(Q,u,akf,rs,sq,sig)                                                                                                                                                                                                                                                                                                        
         implicit double precision (a-h,o-z)                              
         parameter( pi = 3.1415926535897932384626433832795d0)             

         aa = 1.d0/4.d0
         akc = akf/dsqrt(aa)  
         fac1 = 4.d0*pi/(akc**2)
         fac2 = 4.d0*pi/(sq**2)
         if (akc.gt.sq) alxc = fac1
         if (sq.gt.akc) alxc = fac2                                                                                                                             
         rALDAx = -alxc   
         return
         end   
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         Double Precision Function rALDAxc(Q,u,akf,rs,sq,sig)         
         implicit double precision (a-h,o-z)                                           
         parameter( pi = 3.1415926535897932384626433832795d0)                          
                                                                                       
         aa0 = 1.d0/4.d0      
         akc = akf**2/(4.d0*pi)                                                                                                                                                                                         
c        ecu = gu/(1.d0 + b1u*dsqrt(1.d0/rs) + b2u*(1.d0/rs))          
c        ecp = gp/(1.d0 + b1p*dsqrt(1.d0/rs) + b2p*(1.d0/rs)           
c        enfsig1 = (1.d0 + sig)**(4.d0/3.d0)                           
c        enfsig2 =(1.d0 - sig)**(4.d0/3.d0)-2.d0                       
c        fsig = (2.d0**(4.d0/3.d0))-2.d0      
c The uniform electron gas correlation energy according to Perdew and Zunger, Phys. Rev. B, 23, 5076 (1981)                      
         gu =-0.1423d0                                                                                                           
         b1u = 1.0529d0                                                                                                          
         b2u = 0.3334d0                                                                                                          
         cu = 0.0020d0                                                                                                           
         du = -0.0116d0                                                                                                          
         gp = -0.0843d0                                                                                                          
         b1p = 1.3981d0                                                                                                          
         b2p = 0.2611d0                                                                                                                                                                                          
         cp =  0.0007d0                                                                                                 
         dp = -0.0048d0                                                                                                 
         b = (3.d0/(4.d0*pi))**(1.d0/3.d0)                                                                              
         sqrs = dsqrt(rs)                                                                                         
         AA = b*b2u*(21.d0*b1u +16.d0*b2u*sqrs)                     
         BB = 5.d0*b1u+7.d0*(b1u**2)*sqrs+8.d0*b2u*sqrs             
         decdrsn = gu*b**2*(AA + BB*b/rs)                           
         decdrsd = 36.d0*(b*b1u+(b*b2u+b/rs)*sqrs)**3*(b/rs)**3                                                    
         decdrs = decdrsn/decdrsd 
         aa = aa0 - akc*decdrs    
         ac = akf/dsqrt(aa)               
         fac1 = 4.d0*pi/(ac**2)           
         fac2 = 4.d0*pi/(sq**2)                                                                                                                                                                 
         if (akc.gt.sq) alxc = fac1                                                    
         if (sq.gt.akc) alxc = fac2                                                    
         rALDAxc = - alxc                                              
         return                                                                        
         end           
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
        Subroutine CP07(Q,u,akf,rs,sq,sig,fxccp07)   
        implicit double precision (a-h,o-z)                                                                                     
        complex(8):: u, omc, aknn, aknd, akn, fxccp07                    
         parameter( pi = 3.1415926535897932384626433832795d0)  
        
        cfac = 4.d0*pi/(akf**2)
        dpf = 3.d0/(4.d0*pi)
        dp = (1.d0/rs**3)
c The complex frequency in CP07        
         wp = dsqrt(3.d0/rs**3)
         del = 0.000001d0
         u = wp*u
         omc = u + del
c The compressibility sum-rule         
         alph = -0.025518491d0
         b = -0.691590707d0 
         cf = 4.d0*pi  
         fxcl = cf*alph*(dpf*dp)**b          
c The high-frequency limit
         gam1 = -0.114548231d0
         gam2 = -0.614523371d0
         fxch = gam1*(dpf*dp)**gam2         
c cnf and an
         cnf = fxch/fxcl   
         an = 6.d0*dsqrt(cnf)
         xx = dsqrt(rs)         
c bn according to the parametrization of Eq. (7) of Phys. Rev. B 57, 14569 (1998)
         en =1.d0 + 2.15d0*xx + 0.435d0*xx**3.d0 
         den = 3.d0 + 1.57d0*xx + 0.409d0*xx**3.d0
         bn = en/den     
         cc = -pi/2.d0*akf
c The uniform electron gas correlation energy according to Perdew and Zunger, Phys. Rev. B, 23, 5076 (1981)
         gu =-0.1423d0
         b1u = 1.0529d0
         b2u = 0.3334d0
         cu = 0.0020d0
         du = -0.0116d0
         gp = -0.0843d0
         b1p = 1.3981d0
         b2p = 0.2611d0
         cp =  0.0007d0
         dp = -0.0048d0
         ecu = gu/(1.d0 + b1u*dsqrt(rs) + b2u*rs)
         ecp = gp/(1.d0 + b1p*dsqrt(rs) + b2p*rs)    
         enfsig1 = (1.d0 + sig)**(4.d0/3.d0)
         enfsig2 =(1.d0 - sig)**(4.d0/3.d0)-2.d0 
         fsig = (2.d0**(4.d0/3.d0))-2.d0
c        ec = ecu + fsig*(ecp - ecu)
c        rsec = rs*ec 
c The rs-dependent cn 
         rsuff = ((1.d0/2.d0)*b1u*(rs**(-1.d0/2.d0))+b2u)
         rspff = ((1.d0/2.d0)*b1p*(rs**(-1.d0/2.d0))+b2p)
         rsfu = gu/(1.d0 + b1u*dsqrt(rs) + b2u*rs) 
         gurs = gu*rs
         rssu = ((1.d0 + 0.5d0*b1u*rs**(-1.d0/2.d0) + b2u*rs)**(-2.d0))
         rsecu = rsfu - gurs*rssu*rsuff
         rsfp = gp/(1.d0 + b1p*dsqrt(rs) + b2p*rs)
         gprs = gp*rs
         rssp = ((1.d0 + 0.5d0*b1p*rs**(-1.d0/2.d0) + b2p*rs)**(-2.d0))
         rsecp = rsfp - gprs*rssp*rspff         
         cn = rsecu + fsig*(rsecp - rsecu)        
c The frequency-dependent kn wavevector 
         aknn = (1.d0 + an*(-om/wp) + cn*(-(om**2/wp**2)))
         aknd = (1.d0 - (omc**2/wp**2))        
         akn = -(fxcl/cf*bn)*(aknn/aknd)
c The CP07 kernel                                                                 
         vc = 4.d0*pi/sq**2            
         cl = vc*bn    
         zp = akn*(sq**2) 
         fxccp07 = vc*(dexp(-zp)- 1.d0)-cfac*cn/(1.d0 + 1.d0/(sq**2))                 
         return
         end                                                                
