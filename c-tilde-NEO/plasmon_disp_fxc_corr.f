       Program Real_Plasmon                                                                                                                                                                                                    
       implicit double precision (a-h,o-z)                                                                                                                                                                                
        double precision Lindhard                                                                                                                                                                                         
        parameter( pi = 3.1415926535897932384626433832795d0)                                                                                                                                                              
        sig = 0.0d0                                                                                                                                                                                                       
        dQ = 0.05d0                                                                                                                                                                                                       
        dW = 0.005d0                                                                                                                                                                                                       
        ctil = 0.264d0 
        cf = 0.24d0
        bb = 1.1d0                                                                                                                                                                                                    
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
                                                                                                                                                          
c Density parameters                                                                                                                                                                                                      
       rs = 5.62d0                                                                                                                                                                                                        
       akf = ((9.d0*pi/4.d0)**(1.d0/3.d0))/rs                                                                                                                                                                             
       ef = akf**2/2.d0                                                                                                                                                                                                   
       wp = dsqrt(3.d0/rs**3) 

c Low-density c parameter from the second-derivative of ALDAxc
c This "c" garantees the satisfaction of the compressibility sum-rule for xc together
c Derivation can be found in my notes       
       d2exc = ALDAxc(Q,u,akf,rs,sq,sig)                    
       cc = (3.d0*(pi**2)*((b/rs)**3.d0))**(-2.d0/3.d0)     
       clow = (-pi/2.d0)*cc*facx/d2exc
       c1 = cf                         
                                                                                                                                                                                                                                                                                                                                                                                                                                   
c Calculate the plasmon dispersion for a pair of Q and W, (Q=q/kf and W=w/ef with ef = (hbar**2*kf**2)/2)                                                                                                                 
                                                                                                                                                                                                                          
       do iQ = 1, 25                                                                                                                                                                                                     
           Q = dQ*(iQ - 1.d0)                                                                                                                                                                                             
          omega = wp/ef                                                                                                                                                                                                   
          afprev = 1.d0                                                                                                                                                                                                   
          if (iQ.eq.1) W = omega                                                                                                                                                                                          
          if (iQ.eq.1) go to 200                                                                                                                                                                                          
                                                                                                                                                                                                                          
          do iW = 1, 701                                                                                                                                                                                                  
           W = dW*(iW - 1.d0)    
             
c The variables of the Lindhard function   
c The Lindhard function comes from Eqs 7-8 in J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954)                                 
          z = Q/2.d0                                                       
          u = W/(2.d0*Q)                                                   
          zu = dabs(z - u)                                                 
          if (zu.gt.1.d0) drl = Lindhard(z,u,akf,rs)                       
          if (zu.le.1.d0) go to 100                                        
                                                                                                                                                                                                            
c The NEO kernel   
         sq = Q*akf                                                                                                                                                                                                       
         zp = (z**2)/((c1)*szp)                                                                                                                                                                                         
         zn = (z**2)/((c1)*szn)                                                                                                                                                                                         
         fzp = 1.d0 - dexp(-zp)                                                                                                                                                                                           
         fzn = 1.d0 - dexp(-zn)                                                                                                                                                                                                                                                                                                                                                                                               
         vc = 4.d0*pi/sq**2                                                                                                                                                                                               
         cl = -vc    
         fxc = cl*(facp*fzp + facn*fzn)                                                                                                                                                                                                 
c         fxcalx = ALDAx(Q,u,akf,rs,sq,sig) 
c           fxcalc = ALDAxc(Q,u,akf,rs,sq,sig)
c          fxcralx = rALDAx(Q,u,akf,rs,sq,sig) 
c          fxcralc = rALDAxc(Q,u,akf,rs,sq,sig)
c          fxc = fxcralc                                                                                                                                                                    
c If RPA, then fxc = 0.d0                                                                                                                                                                                                 
         fxc = 0.d0                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
c Construct the Hxc kernel and find the zeroes of the dielectric function                                                                                                                                                 
         xch = vc + fxc                                                                                                                                                                                                   
         fqw = xch*drl - 1.d0                                                                                                                                                                                             
         afqwm = dabs(fqw)                                                                                                                                                                                                
         if (afqwm.lt.afprev) omega = W                                                                                                                                                                                   
         if (afqwm.lt.afprev) afprev = afqwm                                                                                                                                                                              
                                                                                                                                                                                                                          
100       enddo                                                                                                                                                                                                           
                                                                                                                                                                                                                          
c         write(*, 2) "Q,", "omega"                                                                                                                                                                                          
200       write(*, 2) Q, omega,clow                                                                                                                                                                                            
2        format(3f15.6)                                                                                                                                                                                                   
                                                                                                                                                                                                                          
        enddo                                                                                                                                                                                                             
                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                              
        end                                                                                                                                                                                                               
                                                                                                                                                                                                                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                                                                          
         Double Precision Function Lindhard(z,u,akf,rs)                                                                                                                                                                   
         implicit double precision (a-h,o-z)                                                                                                                                                                              
         parameter( pi = 3.1415926535897932384626433832795d0)                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                     
c re-formulated in real frequencies according to J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954)                                                                                                                            
                                                                                                                                                                                                                          
         wp = dsqrt(3.d0/(rs**3))                                                                                                                                                                                         
         wp2 = wp**2                                                                                                                                                                                                      
         const = -(3.d0*wp2/(4.d0*pi*akf**2))                                                                                                                                                                             
                                                                                                                                                                                                                          
         zu1 = z - u                                                                                                                                                                                                      
         zu2 = z + u                                                                                                                                                                                                      
                                                                                                                                                                                                                          
         fx1 = ((1.d0/(8.d0*z))*(1-zu1**2))                                                                                                                                                                               
         fx2 = dlog(dabs((zu1 + 1.d0)/(zu1 - 1.d0)))                                                                                                                                                                      
         fx3 = ((1.d0/(8.d0*z))*(1.d0 - zu2**2))                                                                                                                                                                          
         fx4 = dlog(dabs((zu2 + 1.d0)/(zu2 - 1.d0)))                                                                                                                                                                      
         fx = 1.d0/2.d0 + fx1*fx2 + fx3*fx4                                                                                                                                                                               
         Lindhard = const*fx                                                                                                                                                                                              
                                                                                                                                                                                                                          
         return                                                                                                                                                                                                           
         end                
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
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
         if (akc.gt.sq) alx = fac1
         if (sq.gt.akc) alx = fac2                                                                                                                             
         rALDAx = -alx   
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
         if (ac.gt.sq) alxc = fac1                                                    
         if (sq.gt.ac) alxc = fac2                                                    
         rALDAxc = - alxc                                              
         return                                                                        
         end         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                   
































                                     


