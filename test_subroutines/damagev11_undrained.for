C----+-----------------------------------------------------------------+
c     This is an version of VUMAT or rocks with frictional hyteresis
c     not accounted for
C----+-----------------------------------------------------------------+
C	  DO NOT FORGET TO CHANGE THE PATH FOR THE LINE BELOW
C	  IT IS IN THE VFRIC ROUTINE AT THE VERY END OF THE FILE
C            open(UNIT=121,
C     +			file='/home/bhat/research/Damage/run1/rtips.rpt')
C----+-----------------------------------------------------------------+

      subroutine vumat(
     + nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     + steptime, totaltime, dt, cmname, coordmp, charlenght,
     + props, density, straininc, relspininc,
     + tempold, stretchold, defgradold, fieldold,
     + stressold, stateold, enerinternold, enerinelasold,
     +    tempnew, stretchnew, defgradnew, fieldnew,
     + stressnew, statenew, enerinternnew, enerinelasnew )

      include 'vaba_param.inc'

      character*80 cmname

C----+-----------------------------------------------------------------+
c fixed arrays
C----+-----------------------------------------------------------------+

      dimension props(nprops), density(nblock), coordmp(nblock,*),
     +  charlength(nblock), straininc(nblock,ndir+nshr),
     +  relspininc(nblock,nshr), tempold(nblock),
     +  stretchold(nblock,ndir+nshr),
     +  defgradold(nblock,ndir+nshr+nshr),
     +  fieldold(nblock,nfieldv), stressold(nblock,ndir+nshr),
     +  stateold(nblock,nstatev), enerinternold(nblock),
     +  enerinelasold(nblock), tempnew(nblock),
     +  stretchnew(nblock,ndir+nshr),
     +  defgradnew(nblock,ndir+nshr+nshr),
     +  fieldnew(nblock,nfieldv),
     +  stressnew(nblock,ndir+nshr), statenew(nblock,nstatev),
     +  enerinternnew(nblock), enerinelasnew(nblock), epst(ndir+nshr)

C----+-----------------------------------------------------------------+
c local variables
C----+-----------------------------------------------------------------+

      EXTERNAL derivs,rkqs
      PARAMETER (nmax=50,KMAXX=200)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      dimension xp(KMAXX),yp(NMAX,KMAXX)
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/stressu/stress(6),se,sm
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
      common/timeinc/ddt,xitime
      common/defg/yold(nmax),stnew(6),stold(6),dfnew(9),dfold(9)
 	    common/reg/regime
c     CHANGES FOR UNDRAINED DEFORMATION            
      common/pressure/p_f,p_int,biot
      real avg_stress, avg_strain
      dimension avg_stress(1,4), avg_strain(1,4) 

C----+-----------------------------------------------------------------+
      den    = props(1)
      xnu    = props(2)
      es     = props(3)
      xkic   = props(4)
      d0     = props(5)
      xf     = props(6)
      beta   = props(7)
      gammat = props(8)
      xfric  = props(9)
      vm     = props(10)
      xlodot = props(11)
      xxm    = props(12)
      flag   = props(13)

c     CHANGES FOR UNDRAINED DEFORMATION      
c     this is the biot coefficient      
      biot   = props(14)
c     this is the drained bulk modulus of the solid
      K_d    = props(15)
c     this is the bulk modulus of the fluid
      K_f    = props(16)
c     this is the initial porosity used for evaluation
c     of the undrained bulk modulus
      phi_0 = props(17)
c     this is the initial pore fluid pressure 
      p_0 = props(18)
C----+-----------------------------------------------------------------+
	    fitKICD = 5.0e9
	    psi = 0.5*atan(1/xfric)
	    alpha = cos(psi)
      ntens=ndir+nshr
      nvar=1
      eps=1.0e-0
      x2time=totaltime
      g=es/(2.0*(1.0+xnu))
      ddt=dt
      qzer=0.0
      kmax=0
      dxsav=dt/2.0
      pi=2.0*asin(1.0)
      third=1.0/3.0
      a=(3.0*d0/(4.0*pi*alpha**3*xf))**third
      sigf=xkic/sqrt(pi*a)

c     CHANGES FOR UNDRAINED DEFORMATION      
c     this is the solid bulk modulus
      K_s = es / (3.0 * (1 - 2.0 *xnu))
c     this is the undrained bulk modulus
      K_u = K_d  + biot * biot * K_s * K_f 
     1 / (phi_0 * K_s + (biot - phi_0) * K_f)  


C----+-----------------------------------------------------------------+

      if (totaltime.eq.0.0) then


       do i=1,nblock
c         CHANGES FOR UNDRAINED DEFORMATION          

           stressnew(i,1)=stressold(i,1)+
     1         2.0*g*(straininc(i,1)+(xnu/(1.0-2.0*xnu))*
     1         (straininc(i,1)+straininc(i,2)+straininc(i,3)))
     1         - biot * p_0
          stressnew(i,2)=stressold(i,2)+
     1        2.0*g*(straininc(i,2)+(xnu/(1.0-2.0*xnu))*
     1        (straininc(i,1)+straininc(i,2)+straininc(i,3)))
     1        - biot * p_0
          stressnew(i,3)=stressold(i,3)+
     1        2.0*g*(straininc(i,3)+(xnu/(1.0-2.0*xnu))*
     1        (straininc(i,1)+straininc(i,2)+straininc(i,3)))
     1        - biot * p_0
          stressnew(i,4)=stressold(i,4)+2.0*g*straininc(i,4)
          stressnew(i,5)=stressold(i,5)+2.0*g*straininc(i,5)
          stressnew(i,6)=stressold(i,6)+2.0*g*straininc(i,6)
       enddo

       goto 200
      endif

C----+-----------------------------------------------------------------+

      do 100 ki=1,nblock

		  do k1=1,ntens
			 strain(k1)=stateold(ki,k1)
			 xstrain(k1)=straininc(ki,k1)
			 stnew(k1)=stretchnew(ki,k1)
			 stold(k1)=stretchold(ki,k1)
			 dfnew(k1)=defgradnew(ki,k1)
			 dfold(k1)=defgradold(ki,k1)
		  enddo

      yold(1) = stateold(ki,ntens+1)
      xldot   = stateold(ki,ntens+2)
      xkold   = stateold(ki,ntens+3)
      xkdotold   = stateold(ki,ntens+4)
      xKICDold   = stateold(ki,ntens+5)
c     CHANGES FOR UNDRAINED DEFORMATION      
      p_f = stateold(ki,ntens+10)        

      if (xKICDold.lt.xkic) then
        xKICDold = xkic
      endif

      xitime=x2time
      xetime=xitime+dt

      emax=0.0

      do k1=1,ntens
        emax=max(emax,abs(straininc(ki,k1)/dt))
      enddo

      if (emax.gt.0.0) then
        dinit=(sigf/emax)/es/10.0
          else
        dinit=dt/10.0
      endif

      if (dinit.gt.dt/10.0) dinit=dt/10.0

      em = 0.0
      ex = 0.0

      
      call odeint(yold,nvar,xitime,xetime,eps,dinit,qzer,nok,
     +     nbad,derivs,rkqs)


      do k1=1,ntens
        statenew(ki,k1)=strain(k1)+straininc(ki,k1)
      enddo

      if (yold(1).ge.0.9999) yold(1)=0.9999
      if (yold(1).le.d0) yold(1)=d0

      statenew(ki,ntens+1) = yold(1)
      statenew(ki,ntens+2) = xldot

C----+-----------------------------------------------------------------+
c       calculate updated stress
C----+-----------------------------------------------------------------+

      pi=2.0*asin(1.0)
      third=1.0/3.0
      a=(3.0*d0/(4.0*pi*alpha**3*xf))**third
      g=es/(2.0*(1.0+xnu))
      d=yold(1)
      d=d/d0
      f1=pi*alpha**(1.5)*(d**third-1.0+beta/alpha)**(0.5)
      f1=sqrt(1.0-alpha**2.0)/f1

      f2=(sqrt(1.0-alpha**2.0)/(alpha**2.0))*
     +    (1.0/(1.0/d0**(2.0/3.0)-d**(2.0/3.0)))

      f3=(2.0*alpha**0.5/pi)*(d**third-1.0)**0.5

      bb=f1+f2*f3
      aa=xfric*bb + f2*f3
      cc=aa+gammat*sqrt(alpha*d**third)
      ee=bb*cc/sqrt(cc**2-aa**2)

      do i=1,ntens
        estrain(i)=statenew(ki,i)
      enddo

      call upstress

      d=d*d0
      test1=-bb/aa*se
      test2=aa*bb/(cc**2-aa**2)*se
C----+-----------------------------------------------------------------+
	    factx = sqrt(pi*d0*(1.0-xnu)/(alpha**3.0))
      aa1 = aa*factx
      bb1 = bb*factx
      facty = 3.0*(1.0-2.0*xnu)/(2.0*(1.0+xnu))
	    convexloss=1.0e6
C----
	    d=1
C----
      f10=pi*alpha**(1.5)*(d**third-1.0+beta/alpha)**(0.5)
      f10=sqrt(1.0-alpha**2.0)/f10
C----
      f20=(sqrt(1.0-alpha**2.0)/(alpha**2.0))*
     +    (1.0/(1.0/d0**(2.0/3.0)-d**(2.0/3.0)))
C----
      f30=(2.0*alpha**0.5/pi)*(d**third-1.0)**0.5
C----
      bb0=f10+f20*f30
      aa0=xfric*bb0 + f20*f30
C----
      aa10 = aa0*factx
      bb10 = bb0*factx
C----
      convexlossD = (facty+0.5*facty*bb1**2*+0.5*aa1**2.0)
	    convexloss0 = (facty+0.5*facty*bb10**2*+0.5*aa10**2.0)
C----
	    convexloss = convexloss0/convexlossD

      if (sm.le.test1) then
         xk=0.0
      else
       	if (sm.lt.test2) then
          xk=sqrt(pi*a)*(aa*sm+bb*se)
       	else
          xk=sqrt(pi*a*(cc**2*sm**2+ee**2*se**2))
       	endif
      endif

      statenew(ki,ntens+3) = xk
      statenew(ki,ntens+4) = xkdot
      statenew(ki,ntens+5) = xKICD
      statenew(ki,ntens+6) = regime
      statenew(ki,ntens+7) = aa
      statenew(ki,ntens+8) = bb 
      statenew(ki,ntens+9) = convexloss
c     CHANGES FOR UNDRAINED DEFORMATION      
c     update the pore fluid pressure for next time step    
      p_f_update = p_0 - (K_u - K) * (estrain(1) 
     1                     + estrain(2) + estrain(3) )/ biot
c     at this point p_f_update should be equal to p_int

c     avoid the pressure to drop to negative values
c     in this case replace by the previous value     
      if (p_f_update.le.0) then                        
        statenew(ki,ntens+10) = stateold(ki,ntens+10)
      else
        statenew(ki,ntens+10) = p_f_update
      endif


      do k1=1,ntens
          stressnew(ki,k1)=stress(k1)
      enddo


 100   continue


 200   continue


C----+-----------------------------------------------------------------+
c       calculate the average stress and strain
C----+-----------------------------------------------------------------+ 

c       do ki = 1, nblock 
c        if (ki.eq.1) then
cc         Get the average stress and strain in the sample          
c          call average_stress_strain(stressnew, statenew , nblock,
c     1    ntens, nstatev, avg_stress, avg_strain )
cc         update state variables          
c          do i = 1, ntens
c          statenew(ki,ntens + 10 + i) = avg_stress(1,i)
c          statenew(ki,ntens + 14 + i) = avg_strain(1,i)
c          enddo
c        else 
cc         use the average stress and strain field
cc         calculated from the first element to update the otehr
cc         elements
c
cc         update state variables          
c          do i = 1, ntens
c          statenew(ki,ntens + 10 + i) = statenew(1,ntens + 10 + i)
c          statenew(ki,ntens + 14 + i) = statenew(ki,ntens + 10 + i)
c          enddo
c
c        endif
c      enddo

      do ki = 1, nblock 
        if (ki.eq.1) then
c         update state variables          
          do i = 1, ntens
          statenew(ki,ntens + 10 + i) = 100
          statenew(ki,ntens + 14 + i) = 200
          enddo
        else 
c         use the average stress and strain field
c         calculated from the first element to update the otehr
c         elements

c         update state variables          
          do i = 1, ntens
          statenew(ki,ntens + 10 + i) = statenew(1,ntens + 10 + i)
          statenew(ki,ntens + 14 + i) = statenew(1,ntens + 10 + i)
          enddo

        endif
      enddo


      return
      end
C----+-----------------------------------------------------------------+
      subroutine derivs(x,y,dydx)
      include 'vaba_param.inc'
      parameter (nmax=50)
      dimension y(nmax),dydx(nmax)
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/timeinc/dt,xitime
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/stressu/stress(6),se,sm
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
c     CHANGES FOR UNDRAINED DEFORMATION            
      common/pressure/p_f,p_int,biot




      pi=2.0*asin(1.0)
      third=1.0/3.0
      a=(3.0*d0/(4.0*pi*alpha**3*xf))**third
      g=es/(2.0*(1.0+xnu))

 	    d=y(1)

c     keep the damage parameters in the proper bounds
	    if (d.ge.0.9999) d=0.9999
	    if (d.le.d0) d=d0

c     normalize the damage parameter
	    d=d/d0

      f1=pi*alpha**(0.5)*(d**third-1.0+beta/alpha)**(0.5)
      f1=sqrt(1.0-alpha**2.0)/f1

      f2=(sqrt(1.0-alpha**2.0)/(alpha**2.0))*
     +    (1.0/(1.0/d0**(2.0/3.0)-d**(2.0/3.0)))

      f3=(2.0*alpha**0.5/pi)*(d**third-1.0)**0.5

      bb=f1+f2*f3
      aa=xfric*bb + f2*f3
      cc=aa+gammat*sqrt(alpha*d**third)
      ee=bb*cc/sqrt(cc**2-aa**2)

c     at each sub-time step in the RK solver, evalutate the
c     new strain. (e.g. divide the strain increment for RK2)
      do i=1,4
      	estrain(i)=strain(i)+(x-xitime)/dt*xstrain(i)
      enddo
c     CHANGES FOR UNDRAINED DEFORMATION            
c     update the pore fluid pressure
      p_int = p_0 - (K_u - K) * (estrain(1) 
     1                     + estrain(2) + estrain(3) )/ biot


      call upstress

      d=d*d0
      test1=-bb/aa*se
      test2=aa*bb/(cc**2-aa**2)*se

      if (sm.le.test1) then
         xk=0.0
      else
       if (sm.lt.test2) then
          xk=sqrt(pi*a)*(aa*sm+bb*se)
       else
          xk=sqrt(pi*a*(cc**2*sm**2+ee**2*se**2))
       endif
      endif

      xkdot = ((xk-xkold))/(xetime-xitime)
C----+-----------------------------------------------------------------+
C			CRACK GROWTH LAW
C----+-----------------------------------------------------------------+
      if (flag.eq.1) then
C----+-----------------------------------------------------------------+
C			DYNAMIC CRACK GROWTH LAW
C----+-----------------------------------------------------------------+

          if (xldot.lt.1.0e-2) then

             xKICD = xkic*(1.0+xkdot/fitKICD)

             if (xk.ge.xKICD) then
                call getrupturevelocity
             else
                xldot = 0.0
             end if

          else

             xKICD = xkic*(1.0+xkdot/fitKICD)

             if (xk.ge.xKICD) then
                call getrupturevelocity
             else
                xldot = 0.0
             end if

          end if

	  else if (flag.eq.2) then
C----+-----------------------------------------------------------------+
C			DYNAMIC CHARLES CRACK GROWTH LAW
C----+-----------------------------------------------------------------+
          if (xk.lt.0.0) call xit('derivs',2)

          if (xk.ge.xkic) then
C             Paliwal and Ramesh JMPS 2008 eqn. 31
			  xldot = 0.92*sqrt(g/den)*((xk-xkic)/(xk-0.5*xkic))
		  else
			  xldot = 0.0
		  end if


      else
C----+-----------------------------------------------------------------+
C			CHARLES CRACK GROWTH LAW
C----+-----------------------------------------------------------------+
          if (xk.lt.0.0) call xit('derivs',2)

            temp=xk/xkic
            temp1=(sqrt(g/den))**(1.0/xxm)

            if (xlodot**(1.0/xxm)*temp.gt.temp1) then
               xldot=sqrt(g/den)
            else
               xldot=xlodot*temp**xxm
            end if

      end if
C----+-----------------------------------------------------------------+


      dydx(1)=3.0*d**(2.0/3.0)*(xldot/a)*d0**third/alpha

      if (d.ge.0.9999) dydx(1)=0.0

      return
      end
C----+-----------------------------------------------------------------+
C----+	SOLVE FOR RUPTURE VELOCITY (NEWTON-RAPHSON)
C----+-----------------------------------------------------------------+
      subroutine getrupturevelocity
      include 'vaba_param.inc'
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
      common/timeinc/ddt,xitime
      common/ntensor/nndir,nnshr,ntens


C----
C----
      g=es/(2.0*(1.0+xnu))
C----
      CS = sqrt(g/den)
      CR = 0.92*CS
      a00 = 1.0/vm**5.0
	    b00 = xk/(xKIC*CR)
      c00 = (1.0-xk/xKIC)
C----
      DL = 1.0E-4
C----
      x00 = 0.5*vm
      DX = 1.0
C----
      DO 164  WHILE (ABS(DX).GT.DL)
		    FX0  = a00*x00**5+b00*x00+c00
		    DFX0 = 5.0*a00*x00**4+b00
        x10 = x00 - FX0/DFX0
        DX = x10 - x00
        x00 = x10
  164 END DO
C----
	    if (x10.gt.vm) then
	     xldot = vm
      else if (x10.lt.0.0) then
         xldot = 0.0
      else
         xldot = x10
      end if
C----
C----
      return
      end
C----+C----+-----------------------------------------------------------------+

      subroutine upstress
      include 'vaba_param.inc'
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
      common/reg/regime
c     CHANGES FOR UNDRAINED DEFORMATION            
      common/pressure/p_f,p_int,biot



      g=es/(2.0*(1.0+xnu))
      pi=2.0*asin(1.0)

      em=0.0
      do i=1,3
       em=em+estrain(i)
      enddo

      ex=0.0
      do i=1,3
       ex=ex+(estrain(i)-em/3.0)**2
      enddo

      do i=1,3
       ex=ex+2.0*estrain(i+3)**2
      enddo

      ex=sqrt(2.0/1.0*ex)


      fac = pi*d0*(1.0-xnu)/(alpha**3)
      ee = bb*cc/sqrt(cc**2-aa**2)
C----
      test1=-3.0*ex*bb*(1.0-2.0*xnu)/(2.0*aa*(1.0+xnu))
C----
      test2=3.0*(1.0-2.0*xnu)/(2.0*(1.0+xnu))+fac*cc**2/2.0
      test2=test2/(1.0+fac*ee**2/2.0)
      test2=test2*aa*bb/(cc**2-aa**2)*ex

      if (em.le.test1) then
       regime = 1.0
       call hooke
      else
       if (em.lt.test2) then
       	  regime = 2.0
          call ashbys
       else
       	  regime = 3.0
          call quadk
       endif
      endif

      return
      end
C----+-----------------------------------------------------------------+
      subroutine hooke
      include 'vaba_param.inc'
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/stressu/stress(6),se,sm
c     CHANGES FOR UNDRAINED DEFORMATION            
      common/pressure/p_f,p_int,biot

      g=es/(2.0*(1.0+xnu))

          
      sm=2.0*(1+xnu)*g/(3.0*(1.0-2.0*xnu))*em
      se=g*ex

      do i=1,3
c     CHANGES FOR UNDRAINED DEFORMATION            
        stress(i) = 2.0 * g * estrain(i) 
     1            + 3.0 * xnu * sm / (1.0 + xnu) 
     1            - biot * (p_int )
      enddo

      do i=1,3
       stress(i+3)=2.0*g*estrain(i+3)
 	    enddo

      return
      end
C----+-----------------------------------------------------------------+
      subroutine quadk
      include 'vaba_param.inc'
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/stressu/stress(6),se,sm
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
c     CHANGES FOR UNDRAINED DEFORMATION            
      common/pressure/p_f,p_int,biot


      g=es/(2.0*(1.0+xnu))
      pi=2.0*asin(1.0)
C----
      fac=3.0*(1.0-2.0*xnu)/(2.0*g*(1.0+xnu))
      fac=fac+pi*d0*(1.0-xnu)*cc**2/(2.0*alpha**3*g)
      sm=em/fac
C----
      fac1=1.0/g+pi*d0*(1.0-xnu)*ee**2/(2.0*alpha**3*g)
      se=ex/fac1
C----
      fac1 = sqrt(pi*d0*(1.0-xnu)/(alpha**3))
      cc1 = cc*fac1
      ee1 = ee*fac1
C----
	    temp0 = 1.0+ee1**2/2.0
	    temp1 = 3.0*xnu/(1.0+xnu) + ee1**2/2.0 - cc1**2/3.0

C----
      do i=1,3
c     CHANGES FOR UNDRAINED DEFORMATION            
       stress(i)=(2.0*g*estrain(i) + temp1*sm)/temp0 - biot*(p_int)
      enddo
C----
      do i=1,3
       stress(i+3)=(2.0*g*estrain(i))/temp0
 	    enddo

      return
      end
C----+-----------------------------------------------------------------+
      subroutine ashbys
      include 'vaba_param.inc'
      common/ceramic/den,xnu,es,xkic,d0,
     +        xf,alpha,beta,gammat,xfric,vm,xlodot,xxm,flag,K_d,K_u
      common/state/strain(6),xstrain(6),estrain(6),em,ex
      common/stressu/stress(6),se,sm
      common/params/aa,bb,cc,ee,xldot,xKICD,xkold,xk,xkdot
      common/params2/fitKICD
c     CHANGES FOR UNDRAINED DEFORMATION      
      common/pressure/p_f,p_int,biot


      g=es/(2.0*(1.0+xnu))
      pi=2.0*asin(1.0)

      fac3 = sqrt(pi*d0*(1.0-xnu)/(alpha**3))
	    aa1 = aa*fac3
	    bb1 = bb*fac3

      det=3.0*(1.0-2.0*xnu)/(2.0*(1.0+xnu))+aa1**2/2.0+
     1    3.0*(1.0-2.0*xnu)/(4.0*(1.0+xnu))*bb1**2

      tempa1 = (1+bb1**2.0/2.0)/det
      tempb1 = -(aa1*bb1/2.0)/det
      tempb2 = (aa1**2/2 + 3.0*(1.0-2.0*xnu)/(2.0*(1.0+xnu)))/det

      sm = g*(tempa1*em+tempb1*ex)
      se = g*(tempb1*em+tempb2*ex)

	    tempfacsig = (3.0*xnu*se/(1.0+xnu) + aa1*bb1*sm/(2.0))
	    tempfacsig = tempfacsig - aa1**2*se/3.0 + bb1**2*se/2.0
	    tempfacsig = tempfacsig*sm

	    tempfactau = aa1*bb1*se**2/3.0

	    tempfacs = (se+aa1*bb1*sm/(2.0) + bb1**2*se/2.0)


      do i=1,3
c     CHANGES FOR UNDRAINED DEFORMATION            
       stress(i) = 2.0*g*se*estrain(i) + tempfacsig - tempfactau
       stress(i) = stress(i)/tempfacs - biot*(p_int)
      enddo

      do i=1,3
       stress(i+3) = 2.0*g*se*estrain(i+3)/tempfacs
 	  enddo

      return
      end

C----+-----------------------------------------------------------------+
C			RK ODE Solver from Numerical Recipes
C----+-----------------------------------------------------------------+

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *rkqs)
      include 'vaba_param.inc'
      dimension ystart(nvar)
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=1000000,NMAX=50,KMAXX=200,TINY=1.e-15)
      dimension dydx(NMAX),xp(KMAXX),y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
 11    continue

      if (kmax.gt.0) xsav=x-2.*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
 12      continue
        if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
 13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
 14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
 15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) call
     *  xit('stepsize smaller than minimum in odeint',1)
        h=hnext
 16    continue

      call xit('too many steps in odeint',2)
      return
      END

C----+-----------------------------------------------------------------+

      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
      include 'vaba_param.inc'
      dimension dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      dimension yerr(NMAX),ytemp(NMAX)
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)

      h=htry
 1     call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.
      do 11 i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
 11    continue
      errmax=errmax/eps
      if(errmax.gt.1.)then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        xnew=x+h
        if(xnew.eq.x) call xit('stepsize underflow in rkqs',3)
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.*h
        endif
        hdid=h
        x=x+h
        do 12 i=1,n
          y(i)=ytemp(i)
 12      continue
        return
      endif
      END

C----+-----------------------------------------------------------------+

      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
      include 'vaba_param.inc'
      dimension dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      dimension ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     *ytemp(NMAX)
      PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     *B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     *B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     *B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378.,
     *C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648.,
     *DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336.,
     *DC6=C6-.25)

      do 11 i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
 11    continue
      call derivs(x+A2*h,ytemp,ak2)
      do 12 i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
 12    continue
      call derivs(x+A3*h,ytemp,ak3)
      do 13 i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
 13    continue
      call derivs(x+A4*h,ytemp,ak4)
      do 14 i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
 14    continue
      call derivs(x+A5*h,ytemp,ak5)
      do 15 i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     *B65*ak5(i))
 15    continue
      call derivs(x+A6*h,ytemp,ak6)
      do 16 i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
 16    continue
      do 17 i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+DC6*
     *ak6(i))
 17    continue
      return
      END

C----+-----------------------------------------------------------------+
      subroutine xit (errmsg,xicode)
      include 'vaba_param.inc'
      character*(*) errmsg
      parameter (nrout=6)
      write (nrout,*) ' (',errmsg,') error code=',xicode
      stop
      end
C----+-----------------------------------------------------------------+
	subroutine rotat
	include 'vaba_param.inc'
	parameter (nmax=10)
	common/defg/yold(nmax),stnew(6),stold(6),dfnew(9),dfold(9)
	common/state/strain(6),xstrain(6),estrain(6),em,ex
	dimension fnew(3,3),fold(3,3),rnew(3,3),rold(3,3),
     1  snew(3,3),sold(3,3),drot(3,3),ept(3,3),epp(3,3)
	dimension eptnew(3,3),eppnew(3,3)
	dimension imap(3,3),imap1(3,3)

	imap(1,1)=1
	imap(1,2)=4
	imap(1,3)=9
	imap(2,1)=7
	imap(2,2)=2
	imap(2,3)=5
	imap(3,1)=6
	imap(3,2)=8
	imap(3,3)=3

	imap1(1,1)=1
	imap1(1,2)=4
	imap1(1,3)=6
	imap1(2,1)=4
	imap1(2,2)=2
	imap1(2,3)=5
	imap1(3,1)=6
	imap1(3,2)=5
	imap1(3,3)=3

	do k1=1,3
	do k2=1,3
	fnew(k1,k2)=dfnew(imap(k1,k2))
	fold(k1,k2)=dfold(imap(k1,k2))
	snew(k1,k2)=stnew(imap1(k1,k2))
	sold(k1,k2)=stold(imap1(k1,k2))
	ept(k1,k2)=strain(imap1(k1,k2))
	epp(k1,k2)=yold(imap1(k1,k2))
	rold(k1,k2)=0.0
	rnew(k1,k2)=0.0
	drot(k1,k2)=0.0
	eptnew(k1,k2)=0.0
	eppnew(k1,k2)=0.0
	enddo
	enddo


	call kinv(snew,3,3)
	call kinv(sold,3,3)

	do i=1,3
	do j=1,3
	do m=1,3
	rnew(i,j)=rnew(i,j)+fnew(i,m)*snew(m,j)
	rold(i,j)=rold(i,j)+fold(i,m)*sold(m,j)
	enddo
	enddo
	enddo

	do k1=1,3
	do k2=1,3
	drot(k1,k2)=rnew(k1,k2)-rold(k1,k2)
	enddo
	enddo

    	do k1=1,3
	drot(k1,k1)=drot(k1,k1)+1.0
	enddo


	do i=1,3
	do j=1,3
	do m=1,3
	do l=1,3
	eptnew(i,j)=eptnew(i,j)+drot(m,i)*ept(m,l)*drot(l,j)
	eppnew(i,j)=eppnew(i,j)+drot(m,i)*epp(m,l)*drot(l,j)
	enddo
	enddo
	enddo
	enddo

   	do k1=1,3
	do k2=1,3
	strain(imap1(k1,k2))=eptnew(k1,k2)
	yold(imap1(k1,k2))=eppnew(k1,k2)
	enddo
	enddo


	return
	end



      subroutine kinv(a,n,ndim)
      include 'vaba_param.inc'
      dimension a(3,3),is(60),js(60)
      do 100 k=1,n
         d=0.
         do 10 i=k,n
         do 10 j=k,n
               if(abs(a(i,j)).gt.d)then
                  d=abs(a(i,j))
                  is(k)=i
                  js(k)=j
               endif
 10      continue
         do 30 j=1,n
            t=a(k,j)
            a(k,j)=a(is(k),j)
            a(is(k),j)=t
 30      continue
         do 40 i=1,n
            t=a(i,k)
            a(i,k)=a(i,js(k))
            a(i,js(k))=t
 40      continue
         a(k,k)=1./a(k,k)
         do 50 j=1,n
            if(j.ne.k)then
               a(k,j)=a(k,j)*a(k,k)
            endif
 50      continue
         do 70 i=1,n
            if(i.ne.k)then
               do 60 j=1,n
                  if(j.ne.k)then
                     a(i,j)=a(i,j)-a(i,k)*a(k,j)
                  endif
 60            continue
            endif
 70      continue
         do 80 i=1,n
            if(i.ne.k)then
               a(i,k)=-a(i,k)*a(k,k)
            endif
 80      continue
 100  continue
c
      do 130 k=n,1,-1
         do 110 j=1,n
            t=a(k,j)
            a(k,j)=a(js(k),j)
            a(js(k),j)=t
 110     continue
         do 120 i=1,n
            t=a(i,k)
            a(i,k)=a(i,is(k))
            a(i,is(k))=t
 120     continue
 130  continue
      return
      end
C----+-----------------------------------------------------------------+
C 							User subroutine VFRIC
C----+-----------------------------------------------------------------+

      subroutine vfric (
C Write only -
     *     fTangential,
C Read/Write -
     *     statev,
C Read only -
     *     kStep, kInc, nContact, nFacNod, nSlvNod, nMstNod,
     *     nFricDir, nDir, nStateVar, nProps, nTemp, nPred, numDefTfv,
     *     jSlvUid, jMstUid, jConSlvid, jConMstid, timStep, timGlb,
     *     dTimPrev, surfInt, surfSlv, surfMast, lContType,
     *     dSlipFric, fStickForce, fTangPrev, fNormal, frictionWork,
     *     shape, coordSlv, coordMst, dircosSl, dircosN, props,
     *     areaSlv, tempSlv, preDefSlv, tempMst, preDefMst )
C
      include 'vaba_param.inc'
C
      dimension props(nProps), statev(nstateVar,nSlvNod),
     1     dSlipFric(nDir,nContact),
     2     fTangential(nFricDir,nContact),
     3     fTangPrev(nDir,nContact),
     4     fStickForce(nContact), areaSlv(nSlvNod),
     5     fNormal(nContact), shape(nFacNod,nContact),
     6     coordSlv(nDir,nSlvNod), coordMst(nDir,nMstNod),
     7     dircosSl(nDir,nContact), dircosN(nDir,nContact),
     8     jSlvUid(nSlvNod), jMstUid(nMstNod),
     9     jConSlvid(nContact), jConMstid(nFacNod,nContact),
     1     tempSlv(nContact), tempMst(numDefTfv),
     2     preDefSlv(nContact, nPred),
     3     preDefMst(numDefTfv, nPred)
C
      character*8 surfInt, surfSlv, surfMast
      parameter ( zero = 0.d0, one = 1.d0 )


C----+-----------------------------------------------------------------+
C     user defined state variable, statev(nstatevar, nslvnod)
C     two state variables will be used in this analysis
C     statev(1, nslvnod) will define the slip
C     statev(2, nslvnod will define the velocity
C     Variables to be used in vfric:
C     slip is defined as slip = slip + dSlip, length ncontact
C     dSlipNslv is the same as dSlipFric, but with lenght nslvnod
C     jCon is a list of node numbers on slave surface in contact
C----+-----------------------------------------------------------------+

      dimension slip(nContact)
      dimension dSlipNslv(nSlvNod)
      dimension jCon(nContact)


      integer lend
      integer rend

      n=1
      p=1
      surfnum = 1

C     Initialize variables
      dSlipNslv(:) = zero
      slip(:) = zero
      tol = 1e-12


C----+-----------------------------------------------------------------+
C     Set values of friction coefficients and dc given by input file
C----+-----------------------------------------------------------------+

      xmus = props(1)
      xmud = props(2)
      dc = props(3)

C----+-----------------------------------------------------------------+
C     create jCon with user defined node numbers to replace jConslvid
C     (jSlvUid-1) is added to jConSlvid to change internal node numbers
C     on slave surface from 1 to nSlvNod to the user defined numbers
C     pass global nod number to jCon(ncontact)
C----+-----------------------------------------------------------------+

C    Update State variables -
C
       if (kinc .eq. 0) then
          do 13 kSlv = 1, nSlvNod
            statev(surfnum,kSlv) =
     + 			  statev(surfnum,kSlv) + dSlipNslv(kSlv)
 13       end do
       end if


       do 25  kcon = 1, nContact
           jCon(kcon) = jSlvUid(JConSlvid(kcon))
           Slip(kcon) =
     +         statev(surfnum,jConSlvid(kcon))+  dSlipFric(1,kcon)
           statev(surfnum,jConSlvid(kcon)) = Slip(kcon)
 25    end do


C----+-----------------------------------------------------------------+
CC                  Main Part of VFRIC
C----+-----------------------------------------------------------------+
CC    define the tangential force that should be applied to the nodes
CC    in contact based on slip weakenning law
CC
CC    if the nodes have not slipped yet, the tangential force is set
CC    to the peak shear force or the force required to prevent
CC    any acceleration at that node
CC
CC    if the nodes have slipped, but the slip is less than Dc, the
CC    critical slip displacement,
CC
CC
C----+-----------------------------------------------------------------+
C     define fTangential for slip from previous increment at nodes that
C     are in contact during this increment
C----+-----------------------------------------------------------------+
        lend = 460
        rend = 540

       do kcon = 1, nContact

             fn = fNormal(kcon)
             fs = fStickForce(kcon)

             if ((kcon.ge.lend).and.(kcon.le.rend)) then
             	xmus=0.3
             	xmud=0.01
             else
                xmus = props(1)
	        xmud = props(2)
	     end if

             fp = xmus*fn
             fr = xmud*fn
             xmu= xmud+(xmus-xmud)*(1.0-(Slip(kcon))/dc)

             if (ABS(Slip(kcon)).LT.tol) then
                ft = min(fp,fs)
             else
              	if (Slip(kcon) .GT. dc) then
              	   ft= min(fr,fs)
             	 else
              	   ft= min(xmu*fn,fs)
              	end if
             end if

             fTangential(1,kcon)= -ft
      end do

      return
      end


C----+-----------------------------------------------------------------+
      subroutine average_stress_strain(stressnew, statenew , nblock,
     1    ntens, nstatev, avg_stress, avg_strain )
      include 'vaba_param.inc'
c      integer nblock, ntens , nstatev
c      real, intent(in), dimension(nblock,ntens) :: stressnew
c      real, intent(in), dimension(nblock,nstatev) :: statenew
c      real, intent(out), dimension(1,ntens) :: avg_stress
c      real, intent(out), dimension(1, ntens) :: avg_strain

      integer nblock, ntens, nstatev
      real stress new, statenew, avg_stress, avg_strain
      dimension stressnew(nblock,ntens), statenew(nblock,nstatev),
     1          avg_stress(1,ntens), avg_strain(1,ntens) 


C----+-----------------------------------------------------------------+
C                  SUBROUTINE: FOR EVALUATION OF 
C                  THE AVERAGE STRESS AND STRAIN
C                  IN THE SAMPLE
C----+-----------------------------------------------------------------+
      
c     initialize the average stress and strain
      do i = 1, ntens
        avg_stress(1, i) = 0.0
        avg_strain(1,i) = 0.0
      enddo       

      do i = 1, ntens
        do ki = 1, nblock
          avg_stress(1,i) = avg_stress(1,i) + stressnew(ki,i)
          avg_strain(1,i) = avg_strain(1,i) + statenew(ki,i)
        enddo
        avg_stress(1,i) = avg_stress(1,i) / nblock
        avg_strain(1,i) = avg_strain(1,i) / nblock
      enddo

      return
      end