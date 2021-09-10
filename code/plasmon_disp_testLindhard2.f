       Program Plasmon                                                                                                                                                                                                    
       implicit double precision (a-h,o-z)  
       complex(8)::drlc,DLinc,uu,fx,fx1,fx2,fx3,fx4,zu1,zu2,one                                                                                                                                                                                                                                                                                                                                                                     
        parameter( pi = 3.1415926535897932384626433832795d0)                                                                                                                                                                                                                                                                                                                                                                    
        dQ = 0.5d0                                                                                                                                                                                                       
        dW = 0.5d0                                                                                                                                                                                                                                                                                                                                                                                                        
c       m = 1                                                                                                                                                                                                             
c       hbar = 1                                                                                                                                                                                                          
                                                                                                                                                                                                                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                                                                                                                                                                                                                          
c Calculate drl (Q=q/kf and W=w/ef with ef = (hbar**2*kf**2)/2)                                                                                                                 
                                                                                                                                                                                                                          
       do iQ = 2,5                                                                                                                                                                                                     
           Q = dQ*(iQ - 1.d0)                                                                                                                                                                                             
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
          do iW = 1,11                                                                                                                                                                                                 
           W = dW*(iW - 1.d0)    
             
c The variables of the Lindhard function                                  
          z = Q/2.d0                                                       
          u = W/(2.d0*Q)
          uu = cmplx(u, 4.d0)                                                                                                   
          drl1 = DLindhard1(z,u)                       
          drl2 = DLindhard2(z,u)  
          Call  DLindhardc(z,uu, DLinc) 
          drlc = DLinc                                    
          qoverkf = 4*z 
c why qoverkf = 4*z ?
          Wr = W
c          Wi = 0.1d09
          Wi = 4.0d0
200       write(*, 2) qoverkf,Wr,Wi,drl1,drl2,drlc  
2        format(7f10.3)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
100       enddo                                                                                                                                                                                                           
                                                                                                                                                                                                                          
c         write(*, 2) "Q, f"                                                                                                                                                                                          
                                                                                                                                                                                         
                                                                                                                                                                                                  
                                                                                                                                                                                                                          
        enddo                                                                                                                                                                                                             
                                                                                                                                                                                                                          
        stop                                                                                                                                                                                                              
        end                                                                                                                                                                                                               
                                                                                                                                                                                                                          
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                                                                                                          
         Double Precision Function DLindhard1(z,u)                                                                                                                                                                   
         implicit double precision (a-h,o-z)                                                                                                                                                                              
         parameter( pi = 3.1415926535897932384626433832795d0)                                                                                                                                                             
                                                                                                                                                                                                                          
c The imaginary frequency-dependent Lindhard function from PRB 61, 13433 (2000)                                                                                                                                           
c re-formulated in real frequencies according to J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954)                                                                                                                            
                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
          zu1r = z -u                                                                                                                                                                                                      
          zu2r = z + u
          zu1rp = zu1r+1.d0
          zu1rm = zu1r-1.d0
          zu2rp = zu2r+1.d0
          zu2rm = zu2r-1.d0 
cc I didnot understand why using if statement here.                                                                                                                                                                                                     
          if (dabs(zu1rp).lt.1.d-9.or.dabs(zu1rm).lt.1.d-9) fx2r = 0.d0
       if (dabs(zu1rp).gt.1.d-9.and.dabs(zu1rm).gt.1.d-9)
     1        fx2r=dlog(dabs((zu1r+1.d0)/(zu1r-1.d0)))   
          if (dabs(zu2rp).lt.1.d-9.or.dabs(zu2rm).lt.1.d-9) fx4r = 0.d0                 
       if (dabs(zu2rp).gt.1.d-9.and.dabs(zu2rm).gt.1.d-9)
     1        fx4r=dlog(dabs((zu2r+1.d0)/(zu2r-1.d0)))                                                                                                                                                                                                             
          fx1r = ((1.d0/(8.d0*z))*(1.d0 - zu1r**2))                                                                                                                                                                                                                                                                                                                                                            
          fx3r = ((1.d0/(8.d0*z))*(1.d0 - zu2r**2))                                                                                                                                                                                                                                                                                                                                              
          fxr = 1.d0/2.d0 + fx1r*fx2r + fx3r*fx4r                                                                                                                                                                               
          DLindhard1 = fxr   
                                                                                                                                                                                                                                                                                                                                                                                                                  
         return                                                                                                                                                                                                           
         end         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
         Double Precision Function DLindhard2(z,u)                                       
         implicit double precision (a-h,o-z)                                                   
         parameter( pi = 3.1415926535897932384626433832795d0)                                  
                                                                                               
c The imaginary frequency-dependent Lindhard function from PRB 61, 13433 (2000)                
c re-formulated in real frequencies according to J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954) 
                                                                                               
                                                                                                                                                                                                         
          zu1r = z - u                                                                          
          zu2r = z + u    
          pic = pi/(8.d0*z)   
          
          if (dabs(zu1r).gt.1.d0) f2 = 0.d0
          if (zu2r.lt.1.d0) f2 = (pi/2.d0)*u
          if (zu2r.gt.1.d0.and.dabs(zu1r).lt.1.d0) f2=pic*(1.d0-zu1r**2)                                                                                                                                                                                       
          DLindhard2 = f2                                                                 
                                                                                                   
          return                                                                                             
          end              
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        Subroutine DLindhardc(z,uu, DLinc)  
         implicit double precision (a-h,o-z)                                                                                  
         complex(8)::DLinc,fx,uu,fx1,fx2,fx3,fx4,zu1,zu2,one,half,eight 
          complex(8):: zz                                                       
c         parameter( pi = 3.1415926535897932384626433832795d0)                                                                 
                                                                                                                              
c The complex frequency-dependent Lindhard function from PRB 61, 13433 (2000)                                                 
c re-formulated in real frequencies according to J. Lindhard Dan. Mat. Fys. Medd. 28, 8 (1954)                                
                                                                                                                              
         one = cmplx(1.d0, 0.d0)
         half = cmplx(1.d0/2.d0, 0.d0)  
         eight = cmplx(8.d0, 0.d0) 
         zz = cmplx(z, 0.d0)                                                                                                                                                                                                                                                                                                                                                                                                        
         zu1 = zz - uu                                                                                                          
         zu2 = zz + uu                                                                                                          
                                                                                                                              
         fx1 = ((one/(eight*zz))*(one -zu1**2))                                                                    
         fx2 = cdlog(zu1 + one) - cdlog(zu1 - one)                                                                                 
         fx3 = ((one/(eight*zz))*(one - zu2**2))                                                                               
         fx4 = cdlog(zu2 + one) - cdlog(zu2 - one)                                                                                 
         fx = half + fx1*fx2 + fx3*fx4                                                                                               
         DLinc = fx                                                                                                        
                                                                                                                              
         return                                                                                                               
         end                                                                                                                  
                                                                                                                              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc                                                                






                                                                                                                 
                               















  
  
  
  
  
  
  
  
  
  

                                                                               









                                                                                                                                                                                                                                                      








                                       









                                                                                                





         
         
         
         
         
         
         
         
                                                                                         



 
  




                                                                        





































































                                     


