module surch

    implicit real*8(a-h,o-x)
    implicit character*15(y)
    implicit complex*16(z)

    contains

!===========================================================================================
!===========================================================================================

subroutine main(zp,pmass,zt,amass,emin,emax,edel,lmin,lmax,w0d,ri0d,aid,v0d,rr0d,ard,&
	     wd,rid,ad,rc0,kout,dr,rmax,yy,iampxs,thdmin,thdmax,delth,ivread,fnrv,icoul,&
           iwread,fnrw,flinv,flinw,flout,iwavef,rmax2,kpot,ex_deg,ex_ruthr,icount,&
           nri0,nai,amc12,pm,hc,alf,irdim,ldim,ithdim,i1,i2,i3,i4,i5,i6,&
	     w0max,ri0max,aimax,v0max,rr0max,armax,error_min)
    
      implicit double precision (a-h,o-z)
      allocatable error_min(:,:,:,:,:,:)
      complex*16 q(800),y(0:800),s(181),fnc(181)
      dimension v(800),w(800),vc(800)
      !dimension axsec(181),rxsec(181),ruthr(181),thetad(181)
      dimension rv(800),vext(800),rw(800),wext(800)
      dimension ex_deg(400),ex_ruthr(400)
      character*15 flinv,flinw,flout,y_outfile
      integer iemax

      allocatable axsec(:),rxsec(:),ruthr(:),thetad(:)

      n=int((thdmax-thdmin)/delth)

      allocate(axsec(n),rxsec(n),ruthr(n),thetad(n))


      q=(0.d0,0.d0)
      y=(0.d0,0.d0)
      s=(0.d0,0.d0)
      fnc=(0.d0,0.d0)
      v=0.d0
      w=0.d0
      vc=0.d0
      thetad=0.d0
      axsec=0.d0
      rxsec=0.d0
      ruthr=0.d0
      rv=0.d0
      vext=0.d0
      rw=0.d0
      wext=0.d0


      !==================
      v0=v0max-dble(i1-1)*v0d
      rr0=rr0max-dble(i2-1)*rr0d
      ar=armax-dble(i3-1)*ard

      
      !------------------- 12/1 add -----------------------------

      w0=w0max-dble(i4-1)*w0d
      ri0=ri0max-dble(i5-1)*ri0d
      ai=aimax-dble(i6-1)*aid
    
      !   =======================================================================
      !  >>>> icoul = 1    coulomb potential is calcuated in this program
      !  >>>>       = 0    coulomb potential is included in the external input
      !   =======================================================================
      
      irmax=int(rmax/dr+1.001)
      nthmax=int((thdmax-thdmin)/delth)
      !nthmax=int((thdmax-thdmin)/delth+1.001)

      !print*,'nthmax=',nthmax
      
      !   ========================================================================
      !  >>>> irmax is radial counter
      !  >>>> nthmax is angular counter
      !   ========================================================================
      
      if(iwavef.eq.1) irmax2=int(rmax2/dr+1.001)
      
      open(4,file=flout)
      
      
      ! ---- dimension check --------------------
      
      if(irmax.gt.irdim.or.lmax.gt.ldim.or.nthmax.gt.ithdim) stop
      
      ! -----------------------------------------
      !
      ! ---- external input of real potential ----
      !
      if(ivread.eq.1) then
            open(2,file=flinv)
      !
      !          read(2,*) ivmax,drv
      !  -----newly defined 98. 1.10---------------
            read(2,899) drv,ivmax
            read(2,898) idummya
        899 format(30x,f10.0,9x,i5) 
        898 format(i1)
      !      write(6,899) drv,ivmax
      !      write(6,898) idummya
            rvmax=drv*(ivmax-1)
      !  ---newly defined 98. 9.30----------
            do ir=1,ivmax
      !         read(2,*) vext(ir)
                rv(ir)=drv*(ir-1)
      !          write(6,*) rv(ir),vext(ir)
            end do  
            read(2,868) (vext(ir),ir=1,ivmax)
        868 format(1p6d12.5)
      !  -----------------------------------
      !  ------------------------------------------
      !
      !     --- interporation ---
            do ir=1,irmax+3
                r0=dr*ir
                if(r0.le.rvmax) then
                    v(ir)=fnrv*suphod(r0,rv,vext,ivmax,irdim)
      !             v(ir)=vext(ir+1)
                else
                    v(ir)=0.0
                end if
            end do
      !
      !     --- warning message -------
      !     ---------------------
      !
      end if
      ! ------------------------------------------
      !
      !
      ! ---- external input of imaginary potential ----
      !
      if(iwread.eq.1) then
            open(3,file=flinw)
      !
      !      read(3,*) iwmax,drw
      !  ------newly defined 98. 1.10------------
            read(3,899) drw,iwmax
            read(3,898) idummyb
            rwmax=drw*(iwmax-1)
            do ir=1,iwmax
                read(3,*) wext(ir)
                rw(ir)=drw*(ir-1)
            end do
      !  ----------------------------------------
      !
      !     --- interporation ---
            do ir=1,irmax+3
                r0=dr*ir
                if(r0.le.rwmax) then
                    w(ir)=fnrw*suphod(r0,rw,wext,iwmax,irdim)
      !             w(ir)=wext(ir+1)
                else
                    w(ir)=0.0
                end if
            end do
      
      end if
      ! ------------------------------------------
      !
      ! ---- in case of identical particles ------
      isame=0
      eps=1.0d-10
      zsa=dabs(zp-zt)
      asa=dabs(pmass-amass)
      if(zsa.lt.eps.and.asa.lt.eps) isame=1
            ldel=1
      if(isame.eq.1) ldel=2
      ! ------------------------------------------
      !
            rr =rr0
            rrc=1.3d0*16.d0**(1.d0/3.d0)
            ri =ri0
            rd =rid
      !
      !
            if(ivread.ne.1) then
                do i=1,irmax+3
                    r=dr*i
                    v(i)=ws(r,rr,ar,v0)
                    !v(i)=0.d0
      !              write(6,*) r,v(i),rr,ar
                end do
            end if
      !
            if(iwread.ne.1) then
            do i=1,irmax+3
                r=dr*i
                w(i)=ws1(r,ri,ai,w0)!+wsd(r,rd,ad,wd)
            end do
            end if
      
      ! ---- nuclear and coulomb potentials ---------------
      !
            call cpot(zp,zt,rrc,irmax,hc,alf,dr,vc)
      !
      ! ---------------------------------------------------
      !
      !
      !
      ! ---  do loop of energy -------------------
      !
            iemax=int((emax-emin)/edel+1.0001)
      !
      !
            do ie=1,iemax
                elab=edel*(ie-1)+emin
      
                ecm=elab*amass/(amass+pmass)
                fmu=amass*pmass/(amass+pmass)
      !
                fk=sqrt(2.0*fmu*pm*ecm)/hc
                rho1=fk*dr*irmax
                rho2=fk*dr*(irmax+1.0)
                drho=fk*dr
                eta=zp*zt/alf*sqrt(fmu*pm/2.0/ecm)
      !
      !     write(6,80) ecm,fk,eta
      !  80 format(1h ,2x,'ecm=',f10.5,2x,'fk=',f10.5,2x,'eta=',f10.5)
      !
      ! ---------------------------------------------------
      !
      !
      !-----Loop for the partial waves-----
            do l=lmin,lmax,ldel
      
      !   -------------------------------------------------------
            call qfact(icoul,l,ecm,irmax,dr,v,w,vc,fmu,pm,hc,q)
      !   -------------------------------------------------------
      
            do i=1,irmax
                r=dr*i
                vcfugal=l*(l+1)*hc**2/(2.d0*fmu*pm*r**2)
                ur=v(i)+vc(i)
                if(kpot.eq.1) write(8,*)r,ur+vcfugal
                if(kpot.eq.1) write(9,*)r,w(i)
            end do
      !
            if(l.eq.1) then
                qy=2.0d0*yy/dr**2
            else
                qy=0.0
            end if

            !print*,"check"
      !
      !
      !===Fox Goodwin method===
      !   -------------------------------------------------------
            call scat(yy,qy,dr,irmax,q,y)
      !   ===========================================
      !   >>> scat is differential calculation
      !   ===========================================
      !   -------------------------------------------------------
      !
      !
      !===Cal. of the Coul. w.f.===
      !   -------------------------------------------------------
            call coulfg(l,eta,rho1,rho2,drho,f1,g1,f2,g2,fp1,gp1)
      !   -------------------------------------------------------
      !
      !
      !===Cal. of the S-matrix s===
      !   ---------------------------------------------------
            call smat(l,irmax,y,dr,ecm,fk,f1,g1,fp1,gp1,s,ps,iwavef,kout)
      !   ---------------------------------------------------
      !
      !
      !===Out put of wave function===
            if(iwavef.eq.1) then
      !   -------------------------------------------------------------
            call wavef(l,irmax,irmax2,y,dr,fk,eta,f1,g1,fp1,gp1,s,ps)
      !   -------------------------------------------------------------
            end if
      !
      !
      !
            end do
      !-----End loop of the partial wave L-----
      !
            !print*,'check1'
      !
      !===Diff. cross section===
            if(iampxs.eq.1) then
      ! ---------------------------------------------------------
            do nth=1,nthmax
                thetad(nth)=delth*(nth-1)+thdmin !calculation degree
            end do

            !print*,"check2"
      
      !
            call ampxsc(isame,lmin,lmax,fk,eta,s,nthmax,thetad,ecm,fnc,axsec,rxsec,ruthr)
      
            end if
      ! ---------------------------------------------------------
            call auto_surch(thetad,ruthr,ex_deg,ex_ruthr,icount,n,error)

            y_outfile='plot_rdepth.dat'
            ND=44
            !call out_data(y_outfile,thetad,ruthr,n,ND)
      
            !if(error_min > error) then
      
                  !iiimin=iii
                  !jjjmin=jjj
                  !kkkmin=kkk
                  !error_min=error
          
            !end if
      
            error_min(i1,i2,i3,i4,i5,i6)=error
      
            !print*,error
      
      !----------------------------------------------------------------
      
      !
      !===Reaction cross section 
      !-----------------------------------------
            !call absx(lmin,lmax,ldel,fk,s,ecm,absout)
      !-----------------------------------------
      !
            end do
      
            !end do !11/21 add
      
      !------------- 12/1 add -----------------------------------------
            !end do
      
            !print*,iii,w0,error!_min(iii,1,1) !12/1 add
            !100.d0/nw0*
      
            !print*,iii,thetad(120),ruthr(120)

        !end do
    
end subroutine

!===========================================================================================

subroutine auto_surch(calc_deg,calc_ruthr,ex_deg,ex_ruthr,icount,n,error)
    
    implicit real*8(a-h,o-x,z)
    implicit character*15(y)

    dimension ex_deg(400),ex_ruthr(400),con_ruthr(400)
    dimension calc_deg(181),calc_ruthr(181)

    do i=1,icount-1
      !print*,"ex_data=",ex_deg(i),ex_ruthr(i)
    end do

    do i=1,165
      !print*,calc_deg(i),calc_ruthr(i)
    end do

    !print*,"icount=",icount

    do i=1,icount-1

      do j=1,n

          if(calc_deg(j) > ex_deg(i)) then
              if(calc_deg(j-1) <= ex_deg(i)) then
                  deldeg = ex_deg(i) - calc_deg(j-1)
                  druthr = (calc_ruthr(j) - calc_ruthr(j-1)) / (calc_deg(j) - calc_deg(j-1))
                  con_ruthr(i) = calc_ruthr(j-1) + deldeg*druthr
                  !print*,"b=",con_ruthr(i)
              end if
              !exit
          end if

      end do

  end do

  call sum_solidangle(ex_ruthr,ex_deg,icount,ex_sum)
  call sum_solidangle(con_ruthr,ex_deg,icount,calc_sum)

  !print*,"calc=",calc_sum,"ex=",ex_sum


    error=0.d0
    do i=1,icount-1

      f=(ex_ruthr(i)-con_ruthr(i))**2.d0 /sqrt(ex_ruthr(i))/dble(icount-1)
      error = error + f

    end do

    !error = abs(ex_sum - calc_sum)

    !print *, energy
    !print *, error
    !print*, icount-1

end subroutine

!===========================================================================================

subroutine out_data(y_outfile,x,f,n,ND)

    implicit real*8(a-h,o-x)
    implicit character*15(y)
    implicit complex*16(z)

    allocatable x(:),f(:)

    y_file=y_outfile

    open(ND,file=y_file,position='append')

    do i=1,n
        write(ND,*) x(i),f(i)
    end do

    write(ND,'(/)')
    write(ND,'(/)')

    close(ND)

end subroutine

!===========================================================================================

subroutine sum_solidangle(base,high,icount,sum)

    implicit real*8(a-h,o-x)
    implicit character*15(y)
    implicit complex*16(z)

    dimension base(400),high(400),deg(400)

    n = icount-1
    pi = 4.d0*datan(1.d0)

    !allocate(ex_deg(icount),ex_ruthr(icount))

    sum = 0.d0
    do i=1,n-1
      deg(i) = high(i) / 180.d0*pi
      sum = sum + (base(i)+base(i+1))*(high(i+1)-high(i))*0.5*sin(deg(i))
      !print*,high(i),sin(deg(i))
    end do

end subroutine

!===========================================================================================
!===========================================================================================

end module surch

program auto
    !
    !     *****************
    !     *  fox-goodwin  *
    !     *****************
    ! this progaram can calculate external input real(nucler+coul,nucler),
    ! imaginary potential scattering.(pot is fixed in energy)
    ! input data is scatein.data 
    !
    !! the format of the input potential is changed. the format is adjusted to 
    !! the format of the him code.    
    !! 
    !! This code was revised for the studies of the fusion reaction. 
    !!
          !$ use omp_lib

    use surch
    
          implicit double precision (a-h,o-z)
          allocatable error_min(:,:,:,:,:,:),pcount(:)
          dimension ex_deg(400),ex_ruthr(400) 
          character*15 flinv,flinw,flout
          character*15 y_file_ex
          integer absout,nthread

          !allocatable axsec(:),rxsec(:),ruthr(:),thetad(:)
          call cpu_time(time1)

    
    !
    ! ---- physical constants ---
    !     pm=931.49432
    !     hc=197.3216
    !     alf=137.03604
    ! ---------------------------
    !
    !   === same as himut ===
    !      pm=931.200154 d0
    !      hc=197.3216 d0
    !      alf=137.0388 d0
    ! ---------------------------
    !
    !   === same as Three body time dependent===
    !       pm=931.494320d0
    !       hc=197.3270530d0
    !       alf=137.03598950d0
    ! ---------------------------
    !   === same as Variational method ===
          amc12=931.4941D0
          pm=1.00727647D0*amc12
          hc=197.3271D0
          alf=137.035982D0
    ! ---------------------------
    !
          irdim=800
          ldim=120
          ithdim=181
    !
    !
    ! ----- read from 1 -----
          open(1,file='scatein2.dat')
          absout=7
          !open(7,file='absxsec')
          open(8,file='realpot')
          open(9,file='imagpot')
    !
          read(1,11) zp,pmass,zt,amass !
          read(1,12) emin,emax,edel,lmin,lmax !
          read(1,13) v0,rr0,ar !
          read(1,13) w0,ri0,ai !
          read(1,13) wd,rid,ad !
          read(1, *) rc0 !
          read(1, *) kout !
          read(1,15) dr,rmax !
          read(1,14) yy !
          read(1,18) iampxs,thdmin,thdmax,delth !
          read(1,16) ivread,fnrv,icoul !
          read(1,19) iwread,fnrw !
          read(1,17) flinv,flinw,flout !
          read(1,19) iwavef,rmax2
          read(1,2019) kpot
    !
       11 format(5x,4f8.3)
       12 format(3f14.0,2i5)
       13 format(5x,3f8.3)
       14 format(5x,d15.5)
       15 format(5x,2f8.3)
       16 format(i5,f10.3,i5)
       17 format(3a15)
       18 format(i5,3f8.3)
       19 format(i5,f10.3)
    2019  format(i5)
    
        print*,"=============================================="
        y_file_ex='44Ti_ex8-52.txt'
        call read_exfile(y_file_ex,ex_deg,ex_ruthr,icount)
        print*,"read file is ",y_file_ex
        print*,"icount=",icount
        print*,""

    !====================================================================
    ! >>> real pot parameter
    !====================================================================

        v0max=v0
        v0min=180.d0
        v0d=1.d0
        nv0=int((v0max-v0min)/v0d+0.1d0)

        rr0max=rr0
        rr0min=4.0d0
        rr0d=0.1d0
        nrr0=int((rr0max-rr0min)/rr0d+0.1d0)

        armax=ar
        armin=1.0d0
        ard=0.1d0
        nar=int((armax-armin)/ard+0.1d0)
    
    !====================================================================
    ! >>> imag pot parameter
    !====================================================================
    
        w0max=w0
        w0min=0.d0
        w0d=1.d0
        nw0=int((w0max-w0min)/w0d+0.1d0)
    
        ri0max=ri0
        ri0min=1.0d0
        ri0d=0.4d0
        nri0=int((ri0max-ri0min)/ri0d+0.1d0)

        aimax=ai
        aimin=1.0d0
        aid=0.4d0
        nai=int((aimax-aimin)/aid+0.1d0)

      !=========== number of parameters ===============
        npara=6
      !================================================


      !=========== thread number ======================
        nthread=12
      !================================================

        print*,"=============================================="
        print*,"max roop "
        print*,"nv0=",nv0,"nrr0=",nrr0,"nar=",nar
        print*,"nw0=",nw0,"nri0=",nri0,"nai=",nai
        print*,"=============================================="
    
    !----------------------------------------------------------
        allocate(error_min(nv0,nrr0,nar,nw0,nri0,nai),pcount(npara))
    
    !------ initial value errormin ----------------------------
        error_min=100.d0

        !$ write(6,*) "ok"
    
    !==================================================================================
    !==================================================================================
    ! >>>>> openMP parameter
    !==================================================================================
        !$omp parallel &
        !$omp private(i1,i2,i3,i4,i5,i6) &
        !$omp num_threads(nthread)
        !$omp do schedule(dynamic,1)


     !=================================================================================
     ! >>> i1 is real pot depth roop
     ! >>> i2 is real pot radius roop
     ! >>> i3 is real pot diff roop
     ! >>> i4 is imag pot depth roop
     ! >>> i5 is imag pot radius roop
     ! >>> i6 is imag pot diff roop 
     !=================================================================================

        do i1=1,nv0


        do i2=1,nrr0

      
        do i3=1,nar
      
        
        do i4=1,nw0


        do i5=1,nri0


        do i6=1,nai


        call main(zp,pmass,zt,amass,emin,emax,edel,lmin,lmax,w0d,ri0d,aid,v0d,rr0d,ard,&
        wd,rid,ad,rc0,kout,dr,rmax,yy,iampxs,thdmin,thdmax,delth,ivread,fnrv,icoul,&
        iwread,fnrw,flinv,flinw,flout,iwavef,rmax2,kpot,ex_deg,ex_ruthr,icount,&
        nri0,nai,amc12,pm,hc,alf,irdim,ldim,ithdim,i1,i2,i3,i4,i5,i6,&
        w0max,ri0max,aimax,v0max,rr0max,armax,error_min)

        end do


        end do


        end do


        end do


        end do

        write(6,*) 100.d0/nv0*i1


        end do
    
        !$omp end do
        !$omp end parallel
    !----------------------------------------------------------------
    
        open(23,file="44Ti_parameter",status="replace")

        pcount(:)=minloc(error_min)

        v00=v0max-dble(pcount(1)-1)*v0d
        rr00=rr0max-dble(pcount(2)-1)*rr0d
        ar0=armax-dble(pcount(3)-1)*ard
        w00=w0max-dble(pcount(4)-1)*w0d
        ri00=ri0max-dble(pcount(5)-1)*ri0d
        ai0=aimax-dble(pcount(6)-1)*aid

        !print*,'count=',pcount(:)
    
        !print*,minloc(error_min)
    !----------------------------------------------------------------

        print*,"=============================================="
        print*,'v=',v00,'rr=',rr00,'ar=',ar0 
        print*,'w=',w00,'ri=',ri00,'ai=',ai0!12/1 add
        print*,"=============================================="
        print*,'error=',minval(error_min) !11/22 add
    
        write(23,*)'v=',v00,'rr=',rr00,'ar=',ar0
        write(23,*)'w=',w00,'ri=',ri00,'ai=',ai0
        write(23,*)'error=',minval(error_min)

        !y_outfile='test.dat'
        !y_plotfile='testplot.txt'
        !y_gifname='testgif.gif'
        !ND=42
        !call create_giffile(y_outfile,y_plotfile,y_gifname,ND,nw0)
    
        close(23)
    !
        close(1)
        close(2)
        close(3)
        close(4)
        close(7)
        close(8)
        close(9)
        !stop

        call cpu_time(time2)

        print*,time2-time1

	!$ print*,(time2-time1)/dble(12)




    end program
    !process nosource

    !          ***** ampxsc *****
    !
          subroutine ampxsc(isame,lmin,lmax,fk,eta,s,nthmax,thetad,ecm,fnc,axsec,rxsec,ruthr)
    !
          implicit real*8(a-h,o-z)
    !
          complex*16 s(181),fnc(181)
          complex*16 facki,factf,ai
          complex*16 exsig(120),fcl(362)
          dimension sigma(120),thetad(181)
          dimension axsec(181),rxsec(181),ruthr(181)
          data memo1/0/
    !
          pai=3.141592653579893d0
          ldim=120
          mdim=1
          nthdim=181
          nthdm2=nthdim*2
          mmax=0
          lmax1=lmax+1
    
          call sigmal(eta,lmax1,ldim, sigma,exsig)
    
          sigm0=sigma(1)
    
          call fcoul(eta,fk,sigm0,isame,thetad,nthmax,nthdim,nthdm2,fcl)
    
          ai=(0.0d0,1.0d0)
          facki=ai/(2.0d0*fk)
          do nth=1,nthmax
          fnc(nth)=0.0
          end do
    !
    ! ---- scattering amplitudes ----
    !
          ldel=1
          if(isame.eq.1) ldel=2
    !
          pxsc=0.0
          do l=lmin,lmax,ldel
    !
          l1=l+1
          pxsc=pxsc+(2*l+1.0d0)*cdabs(1.0d0-s(l1))**2
          factf=(2*l+1.0d0)*exsig(l1)**2*(1.0d0-s(l1))*facki
    !
          do nth=1,nthmax
          thd=thetad(nth)
          fnc(nth)=fnc(nth)+factf*pl(l,thd)
          end do
    !
          end do
          pxsc=pxsc*pai/fk**2*10.0d0
    
    
    
          do nth=1,nthmax
    
          if(isame.ne.1) then
          fnc(nth)=fnc(nth)+fcl(nth)
          else
          nthp=nth+nthmax
          fnc(nth)=fnc(nth)*2+fcl(nth)+fcl(nthp)
          end if
    
          end do
    !
    ! ---- differential cross sections ----
    
          sumxsc=0.0
          thfac=pai/180.0d0
    
          do nth=1,nthmax
          axsec(nth)=cdabs(fnc(nth))**2*10.0d0
          theta=thetad(nth)*thfac
          sumxsc=sumxsc+axsec(nth)*dsin(theta)
    !
          if(eta.eq.0.0) then
      !410 format(1h ,3x,'th=',f7.2,3x,'axsec=',1pd13.4)
    !
          else
    !
          if(isame.ne.1) then
          rxsec(nth)=cdabs(fcl(nth))**2*10.0d0
          else
          nthp=nth+nthmax
          rxsec(nth)=cdabs(fcl(nth)+fcl(nthp))**2*10.0d0
          end if
    
          ruthr(nth)=axsec(nth)/rxsec(nth)
    !
          end if
    
          !print*,axsec(nth),rxsec(nth)
    !
          end do
    !
    !
    !  --- angle-integrated cross section ---
          delth=(thetad(2)-thetad(1))*thfac
          sumxsc=sumxsc*2.0*pai*delth
          if(memo1.eq.0) then
          memo1=1.0
    !     write(6,510) thetad(1),thetad(nthmax)
      !510 format(1h ,'--- angle-integrated cross section ---'/
      !   &          ,'    thmin=',f5.2,'    through    thmax=',f7.2/)
          end if
    !
    !     write(6,500) ecm,sumxsc
      !500 format(1h ,3x,'ecm=',f10.6,5x,'ai-xsc=',1pe13.4)
    !     write(6,520) ecm,pxsc
          write(4,520) ecm,pxsc
      520 format(1h ,3x,f10.6,5x,1pe13.4)
    !
    ! ------- all is over -----------
    !
          return
          end
    !process source
    !---------------------------------------------------------------------
    
            subroutine read_exfile(y_file_ex,ex_deg,ex_ruthr,icount)
    
              implicit real*8(a-h,o-x,z)
              implicit character*15(y)
    
              dimension ex_ruthr(400),ex_deg(400)
    
              NE=41
    
              open(NE,file=y_file_ex)
    
              !print*,'input file is ',y_file_ex
    
              !do while( .true. )
    
              icount=1
              do while( .true. )
    
                  read(NE, *, iostat=ios) energy,deg,ruthr
                  ex_ruthr(icount)=ruthr
                  ex_deg(icount)=deg
                  if(ios < 0) exit   ! can't use openMP
    
                  !print *, ex_deg(icount),ex_ruthr(icount)
    
                  icount=icount+1
    
              end do
    
              close(NE)
    
              !print*,"icount=",icount
    
    
            end subroutine
    !
    ! ---- legendre polynomial ----
    !
          function pl(l,thd)
    
          implicit real*8(a-h,o-z)
          dimension xpl(200)
    !
    
           xpl(1)=1.0
           pl=xpl(1)
          if(l.eq.0) return
    
           pai=3.141592653579893d0
           theta=thd*pai/180.0d0
           x=dcos(theta)
           xpl(2)=x
           pl=xpl(2)
          if(l.eq.1) return
    !
          do il=2,l
           rl=il
           xpl(il+1)=( (2.0d0*rl-1.0d0)*x*xpl(il)-(rl-1.0)*xpl(il-1) )/rl
          end do
          pl=xpl(l+1)
    !
          return
          end
    !*process   source
    !
    !
    ! ***** wspot *****
    !
    !
    ! --- ws ---
    !
          function ws(r,r0,a0,v0)
    !
          implicit real*8(a-h,o-z)
    !      write(6,*) r,r0,a0,v0
          rou=4.5d0
          x2=(r-r0)/(a0) !10/26 act
          !x2=(r-r0*40.d0**(1.d0/3.d0))/a0 !10/26 inact
          x3=-(r/rou)*(r/rou)
          !xx=(r-r0)/a0
          Ea=15.d0
          alpha=3.625d0-0.0105d0*Ea
          if(x2.ge.150.0) then
          ws=0.0
          else
          !ws1=-v0/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))
    !      write (6,*) r,v0
    !      ws1=0.d0
          !ws2=-(v0*alpha*dexp(x3))/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))
    !      ws2=0.d0
    !      ws=ws1+ws2
          !ws=-v0*(1.d0+alpha*dexp(x3))/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))
          !ws=-v0/((1.d0+dexp(x2))*(1.d0+dexp(x2))) !10/26 inact
          ws=-v0/((1.d0+dexp(x2))*(1.d0+dexp(x2))) !10/26 act
    !      write(6,*) r,ws
          !open(77,file='cpot')
          !write(77,*) r,ws
          endif
          return
          end
    
          function ws1(r,r0,a0,v0)
    !
          implicit real*8(a-h,o-z)
          !rou=4.5d0
    !      write(6,*) r,r0,a0,v0
          x2=(r-r0)/(a0) !10/26 act
          !x2=(r-r0*40.d0**(1.d0/3.d0))/a0 !10/26 inact
          !x3=-(r/rou)*(r/rou)
          !xx=(r-r0)/a0
          Ea=15.d0
          alpha=3.625d0-0.0105d0*Ea
          if(x2.ge.150.0) then
          ws=0.0
          else
          !ws1=-v0/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))  !10/26 inact
          ws1=-v0/((1.0d0+dexp(x2))*(1.0d0+dexp(x2))) !10/26 act
    !      write (6,*) r,v0
    !      ws1=0.d0
          !ws2=-(v0*alpha*dexp(x3))/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))
    !      ws2=0.d0
    !      ws=ws1+ws2
          !ws=-v0*(1.d0+alpha*dexp(x3))/((1.0d0+dexp(x2))*(1.0d0+dexp(x2)))
          !ws1=-v0/(1.d0+dexp(xx))
    !      write(6,*) r,ws
          !open(77,file='cpot')
          !write(77,*) r,ws
          endif
          return
          end
    !
    !
    !*process   source
    !
    ! --- wsd ---
    !
          function wsd(r,r0,a0,v0)
    !
          implicit real*8(a-h,o-z)
          x2=(r-r0)/(a0)
          !Ri=3.122d0
          !ai=0.65d0
          !w0=-25.d0
          if(x2.ge.150.0) then
          ws=0.0
          else
          x1=dexp(x2)
          wsd=-4.0d0*v0*x1/(1.0d0+x1)**2
          !Wsd=w0/((1.d0+EXP((r-Ri)/(2.d0*ai)))*(1.d0+EXP((r-Ri)/(2.d0*ai))))
    !      open(77,file='cipot')
    !      write(77,*) r,wsd
          endif
          return
          end
    !
    !
    ! ***** cpot *****
    !
          subroutine cpot(zp,zt,rr,irmax,hc,alf,dr,vc)
          implicit double precision (a-h,o-z)
          dimension vc(800)
          ic=int(rr/dr+1.001)
          do i=1,ic
          r=dr*i
          vc(i)=zp*zt*hc/(2.0*rr)/alf*(3.0-(r/rr)**2)
    !      write(6,15) r,vc(i)
    !   15 format(1h ,2x,'r=',f10.5,2x,'vc=',f10.5)
    !      write(6,*) r,vc(i)
          end do
          do i=ic+1,irmax+3
          r=dr*i
          vc(i)=zp*zt*hc/alf/r
    !      write(6,16) r,vc(i)
    !   16 format(1h ,2x,'r=',f10.5,2x,'vc=',f10.5)
    !      write(6,*) r,vc(i)
    
          end do
          return
          end
    !
    !
    ! ***** qfact *****
    !
          subroutine qfact(icoul,l,ecm,irmax,dr,v,w,vc,fmu,pm,hc,q)
          implicit double precision (a-h,o-z)
          complex*16 cai,q(800)
          dimension v(800),w(800),vc(800)
          cai=(0.0,1.0)
          do i=1,irmax+3
          r=dr*i
    !     ---------------------------------------------------------------
    !      q(i)=l*(l+1)/r**2-2.0*fmu*pm*(ecm-v(i)-cai*w(i)-vc(i))/hc**2
    !     ---- temporal -------------------------------------------------
           q(i)=l*(l+1)/r**2-2.0*fmu*pm*(ecm-v(i)-cai*w(i)-icoul*vc(i))/hc**2
    !     ---------------------------------------------------------------
    !     write(6,20) r,q(i)
    !  20 format(1h ,2x,'r=',f5.1,2x,'q=',2f10.5)
          end do
          return
          end
    !
    !
    ! ***** scat *****      Fox-Goodwin method
    !
          subroutine scat(yy,qy,dr,irmax,q,y)
          implicit double precision (a-h,o-z)
          complex*16 y(0:800),q(800)
          y(1)=yy
          y(2)=((2.0+5.0*dr**2*q(1)/6.0)*y(1)+dr**2*qy/12.0)/(1.0-dr**2*q(2)/12.0)
          do i=3,irmax+3
          y(i)=((2.0+5.0*dr**2*q(i-1)/6.0)*y(i-1)-(1.0-dr**2*q(i-2)/12.0)*y(i-2))/(1.0-dr**2*q(i)/12.0)
          r=dr*i
    !      write(6,30) r,y(i)
       !30 format(1h ,2x,'r=',f10.5,2x,'y=',2e15.6)
          end do
          return
          end
    !
    !
    ! ***** smat *****
    !
          subroutine smat(l,irmax,y,dr,ecm,fk,f1,g1,fp1,gp1,s,ps,iwavef,kout)
          implicit double precision (a-h,o-z)
          complex*16 cai,s(120),y(0:800),d
    !
          pai=3.141592653579893d0
          cai=(0.0d0,1.0d0)
          xfact=pai*10.0d0/fk**2
    !
          do i=irmax-3,irmax+3
    !      write(6,11) y(i)
       !11 format(1h ,5x,'y=',1p2e15.5)
          end do
          d=-12.0*y(irmax-3)+108.0*y(irmax-2)-540.0*y(irmax-1)&
            +12.0*y(irmax+3)-108.0*y(irmax+2)+540.0*y(irmax+1)
          d=d/(720.0*dr*y(irmax))
          s(l+1)=d*(g1-cai*f1)-fk*(gp1-cai*fp1)
          s(l+1)=s(l+1)/(d*(g1+cai*f1)-fk*(gp1+cai*fp1))
          !sd=atan(aimag(s(l+1))/real(s(l+1)))/2.d0
          ss=cdabs(s(l+1))
          pxsc=xfact*(2*l+1.0d0)*cdabs(1.0d0-s(l+1))**2
    !      ps=dreal(cdlog(s(l+1))/2.0d0/cai)
          ps=real(cdlog(s(l+1))/2.0d0/cai)
          psd=(atan(aimag(s(l+1))/real(s(l+1)))/2.d0)
          !if(ps.lt.0.0) ps=ps+pai
          !psd=ps*180.0d0/pai
    ! 
    !-------a little modified for gnu plot 98.10. 4-------
    !-------a little modified for gnu plot 04. 1.24-------
          if(iwavef.eq.1) go to 2255   
    !
    !      if(kout.eq.1)                         !10/26 inact
    !     &write(6,50) l,ecm,pxsc,s(l+1),ss,psd  !10/26 inact
    !
          if(kout.eq.2) then
              write(6,50) l,ecm,pxsc,s(l+1),1.d0-ss**2,psd
          end if
    !
          if(kout.eq.3) then
              write(6,50) l,ecm,pxsc,s(l+1),(abs(1.d0-ss))**2,psd
          end if
    !
    !   50 format(1h ,i3,1x,f11.6,1x,1pd13.4,3d13.4,0pf9.4)
       50 format(2x,i3,2x,f10.7,1pd13.4,3e13.4,0pf9.4)
    !
     2255 continue 
    !
          if(iwavef.ne.1) go to 2266 
    !
    !      if(kout.eq.1)                          !10/26 inact
    !     &write(6,250) l,ecm,pxsc,s(l+1),ss,psd  !10/26 inact
    !
          if(kout.eq.2) then
              write(6,250) l,ecm,pxsc,s(l+1),1.d0-ss**2,psd
          end if
    !
          if(kout.eq.3) then
              write(6,250) l,ecm,pxsc,s(l+1),(1.d0-ss)**2,psd
          end if
          
    !
      250 format('#',1x,i3,f10.6,1pd13.4,3d13.4,0pf9.4)
     2266 continue
    !-----------------------------------------------------
    !
    !     write(6,50) l,s(l+1),ss,psd,pxsc
    !  50 format(1h ,2x,'l=',i4,2x,'s=',1p2e15.6,2x,'ss=',1pe15.6
    !    1          ,2x,'ps=',0pf10.4,2x,'pxsc='1pd12.4)
          return
          end
    !
    !
    ! ***** wavef *****
    !
          subroutine wavef(l,irmax,irmax2,y,dr,fk,eta,f1,g1,fp1,gp1,s,ps)
          implicit double precision (a-h,o-z)
          complex*16 cai,s(100),y(0:300),cfact,cnorm,wf
    !
          cai=(0.0,1.0)
          irmax1=irmax+1
    !    -----------------------------------
    !     cfact=cdexp(-cai*ps)*cai/2.0d0
          cfact=(1.0, 0.0)
    !    -----------------------------------
    !
          cnorm=cfact*( (g1-cai*f1) - s(l+1)*(g1+cai*f1) )/y(irmax)
    !     write(6,90) cdabs(cfact)
    !  90 format(' ****  |cfact|=',1pd12.4,' ****'/)
    !
          y(0)=0.0
          do i=0,irmax
          r=dr*i
          wf=cnorm*y(i)
          wfabs=cdabs(wf)
          write(6,110) r,wfabs,wf
    !  110 format(f10.5,1p3d15.5)
      110 format(f10.5,1p3e15.5)
          end do
    !
          drh=dr*fk
          irmax1=irmax+1
    !
          do i=irmax1,irmax2
          r=dr*i
          rh1=fk*r
          rh2=rh1+drh
          call coulfg(l,eta,rh1,rh2,drh,ff1,gg1,ff2,gg2,ffp1,ggp1)
          wf=cfact*( (gg1-cai*ff1) - s(l+1)*(gg1+cai*ff1) )
          wfabs=cdabs(wf)
          write(6,110) r,wfabs,wf
          end do
    !
          return
          end
    !
    !
    !  ***** coulfg *****
    !
          subroutine coulfg(l,eta,rho1,rho2,drho,f1,g1,f2,g2,fp1,gp1)
          implicit double precision (a-h,o-z)
    !     implicit real*8(a-h,o-z)
          dimension  f(201),g(201),fp(202),gp(202)
    !
          ll=l+1 +1
          e=eta
          h=drho
          r=rho1
          rp=rho2
          te=e+e
          tf=e**2
          if(ll-50) 20,35,35
       20 elp=50.
          j=50
          go to 45
       35 elp=ll
          j=ll
       45 a=atan(e/elp)
          b=dsqrt(tf+elp**2)
          y=a*(elp-0.5d0)+e*(dlog(b)-1.0d0)-sin(a)/(12.0d0*b)&
           +sin(3.0d0*a)/(360.0d0*b**3)-sin(5.0d0*a)/(1260.0d0*b**5)&
           +sin(7.0d0*a)/(1680.0d0*b**7)-sin(9.0d0*a)/(1188.0d0*b**9)
          k=j-1
          if(j-ll)65,65,70
       65 s1=y
       70 do 100 i=1,k
          elp=elp-1.
          j=j-1
          y=y-atan(e/elp)
      100 continue
          s1=y
          del1=r-te
          rmax=dmax1(10.0d0,(tf+3.0d0+4.0d0*e)*5.0d0/12.0d0)
          del=r-rmax
          if(e-5.)280,130,130
      130 if(abs(del1)-abs(del))140,140,280
      140 del=del1
          if(del)147,145,147
      145 i=2
          go to 150
      147 i=1
      150 x=te
          t1=tf
          t2=t1**2
    !**   t3=e** .666666667
    !**   t9=e** .166666667
          t9=e**(1.0d0/6.0d0)
          t3=t9**4
          t4=t3**2
          t5=t4**2
          t6=t3*t5
          t7=t4*t6
          t8=t3*t7
          y=1.22340402d0*t9*(1.0d0+0.495957017d-1/t4-0.888888889d-2/t1+0.245519918d-2&
          /t6-0.910895806d-3/t2+0.253468412d-3/t8)
          z=-0.707881773d0/t9*(1.0d0-0.172826039d0/t3+0.317460317d-3/t1-&
          0.358121485d-2/t5+0.311782468d-3/t2-0.907396643d-3/t7)
          go to 665
      280 if(e)285,290,285
      285 if(del)310,290,290
      290 x=r
          i=2
          go to 320
      310 x=rmax
          i=1
      320 t1=tf
          t2=x+x
          t3=x-e*dlog(t2)+s1
          t4=e/t2
          ss=1.0d0
          ts=0.0
          sl=0.0
          tl=1.0d0-e/x
          sss=1.0d0
          sts=0.0
          ssl=0.0
          stl=tl
          en=0.0
          do 620 k=1,25
          t5=en+1.
          t6=t5+en
          t7=en*t5
          t8=t6*t4/t5
          t9=(t1-t7)/(t2*t5)
          t5=t8*ss-t9*ts
          ts=t8*ts+t9*ss
          ss=t5
          if(abs(ss/sss)-1.0d-10) 630,630,540
      540 t5=t8*sl-t9*tl-ss/x
          tl=t8*tl+t9*sl-ts/x
          sl=t5
          sss=sss+ss
          sts=sts+ts
          ssl=ssl+sl
          stl=stl+tl
          en=en+1.
      620 continue
      630 t8=sin(t3)
          t9=cos(t3)
          y=sss*t9-sts*t8
          z=ssl*t9-stl*t8
      665 go to (670,810),i
      670 m=1
      671 n=abs(del/h)
          if(n)675,675,680
      675 dx=del
          go to 700
      680 en=n
          dx=del/en
      700 t1=0.5d0*dx
          t2=0.25d0*t1
          t3=te
          do 805 i=1,n
          t4=dx*(t3/x-1.)*y
          x=x+t1
          t5=dx*(t3/x-1.)*(y+t1*z+t2*t4)
          x=x+t1
          t6=dx*(t3/x-1.)*(y+dx*z+t1*t5)
          y=y+dx*(z+(t4+t5+t5)/6.0d0)
          z=z+(t4+4.0d0*t5+t6)/6.0d0
      805 continue
          go to (810,828),m
      810 g(1)=y
          m=2
          del=rp-r
          w=z
          go to 671
      828 gp(1)=y
          t1=tf
          t8=sqrt(1.+t1)
          g(2)=((1./r+e)*g(1)-w)/t8
          gp(2)=((1./rp+e)*y-z)/t8
          t2=1.0d0
          t3=2.0d0
          do 910 i=3,ll
          t4=t2+t3
          t5=t2*t3
          t6=t3*sqrt(t2**2+t1)
          t7=t2*sqrt(t3**2+t1)
          g (i)=(t4*(e+t5/r )*g (i-1)-t6*g (i-2))/t7
          gp(i)=(t4*(e+t5/rp)*gp(i-1)-t6*gp(i-2))/t7
          t2=t2+1.0d0
          t3=t3+1.0d0
      910 continue
    !**   i=l+11
          i=ll+11
          n=r+r+11.0d0
          if(i-n)960,960,950
      950 n=i
      960 y=1.0d-20
          yp=y
          x=y
          xp=x
          z=0.0
          zp=z
          t2=n
     1000 t3=t2+1.0d0
          t4=t2+t3
          t5=t2*t3
          t6=t2*sqrt(t3**2+t1)
          t7=t3*sqrt(t2**2+t1)
          y =(t4*(e+t5/r )*y -t6*z )/t7
          yp=(t4*(e+t5/rp)*yp-t6*zp)/t7
          if(n-ll)1060,1060,1080
     1060 f(n)=y
          fp(n)=yp
          go to 1120
     1080 if(1.0d0-abs(y))1090,1090,1120
     1090 y=y*1.0d-20
          yp=yp*1.0d-20
          x=x*1.0d-20
          xp=xp*1.0d-20
     1120 n=n-1
          z=x
          zp=xp
          x=y
          xp=yp
          t2=t2-1.0d0
          if(n)1150,1150,1000
     1150 y=f(1)*g(2)-f(2)*g(1)
          yp=fp(1)*gp(2)-fp(2)*gp(1)
          z=1.0d0/(y*t8)
          zp=1.0d0/(yp*t8)
          do 1180 i=1,ll
          fp(i)=fp(i)*zp
     1180 f(i)=f(i)*z
          f1=f(l+1)
          f2=fp(l+1)
          g1=g(l+1)
          g2=gp(l+1)
          al1=l+1
          a=al1**2/rho1+eta
          b=al1**2+eta**2
          b=sqrt(b)
          fp1=(a*f1-b*f(l+2))/al1
          gp1=(a*g1-b*g(l+2))/al1
          return
          end
    !*process nosource
    !
    !          ***** sigmal *****
    !
          subroutine sigmal(eta,lmax1,ldim, sigma,exsig)
          implicit real*8(a-h,o-z)
          complex*16 exsig(ldim)
          dimension sigma(ldim)
    !                        ------> sigma(l+1),exp(wi*sigma(l+1))
    !     ****************************
          lmax=lmax1-1
          if(lmax1.gt.ldim) go to 9000
          if(eta.ge.10.0)  go to 10
          eta2=eta*eta
          eta2a=eta+eta
          eta6=eta2+16.0d0
          sigma0=-(eta/(12.0*eta6))*(1.0+(eta2-48.0)/(30.0*(eta6**2))&
                +((eta2-160.0)*eta2+1280.0)/(105.0*(eta6**4)))&
                -eta+(eta/2.0)*dlog(eta6)+3.5*datan(0.25*eta)&
                -(datan(eta)+datan(0.5*eta)+datan(eta/3.0))
          go to 11
       10 einv1=1.0/eta
          einv2=einv1*einv1
          einv3=einv1*einv2
          einv5=einv3*einv2
          einv7=einv5*einv2
          einv9=einv7*einv2
          sigma0=0.7853981634d0+eta*dlog(eta)-eta&
               -(0.08333333333d0*einv1+0.00277777777d0*einv3&
               +0.00079365079d0*einv5+0.00059523810d0*einv7&
               +0.00084175084d0*einv9)
       11 modtpi=sigma0/6.2831853071796d0
          sigmaz=sigma0-6.2831853071796d0*modtpi
          sigma(1)=sigma0
          exsig(1)=dcmplx(dcos(sigmaz),dsin(sigmaz))
          if(lmax.le.0) return
          do 1 ll=1,lmax
          ll1=ll+1
          sigma(ll1)=sigma(ll)+datan(eta/ll)
          modtpi=sigma(ll1)/6.2831853071796d0
          sigmaz=sigma(ll1)-6.2831853071796d0*modtpi
          exsig(ll1)=dcmplx(dcos(sigmaz),dsin(sigmaz))
        1 continue
          return
     9000 continue
          return
    !     debug subchk
          end
    !*process nosource
    !
    !              ***** fcoul *****
    
          subroutine fcoul(eta,fk,sigm0,isame,thetad,nthmax,nthdim,nthdm2,fcl)
          implicit real*8(a-h,o-z)
          complex*16 fcl(nthdm2)
          dimension thetad(nthdim)
          data pai/3.141592653589793d0/
    
          pai2=pai+pai
          dr=pai/180.0d0
    !
          do it=1,nthmax
          theta=thetad(it)*dr
          sine=dsin(theta/2.0d0)  +1.0d-15
          sine2=sine*sine
          fact=-eta/(2.0d0*fk*sine2)
          arg=-eta*dlog(sine2)+2.0d0*sigm0
          mod2pi=int(arg/pai2)
          arg=arg-mod2pi*pai2
          fcl(it)=fact*dcmplx(dcos(arg),dsin(arg))
          end do
    !
          if(isame.ne.1) go to 900
    !
          do it=1,nthmax
          itp=it+nthmax
          theta=thetad(it)*dr
          theta=pai-theta
          sine=dsin(theta/2.0d0)  +1.0d-15
          sine2=sine*sine
          fact=-eta/(2.0d0*fk*sine2)
          arg=-eta*dlog(sine2)+2.0d0*sigm0
          mod2pi=int(arg/pai2)
          arg=arg-mod2pi*pai2
          fcl(itp)=fact*dcmplx(dcos(arg),dsin(arg))
          end do
    !
      900 continue
          return
    !     debug subchk
          end
    !*process nosource
    !
    !
    !
    !
    ! ***  suphod  (real*8)
    !    this program is o.k. only when jisu=3.    ------ care ----
    !    this program is o.k. only when dx=const   ------ care ----
    !
    !
          function  suphod(xin,x,y,nrmax,nrdim)
          implicit real*8(a-h,o-z)
          dimension  x(nrdim),y(0:nrdim)
    !
          jisu=3
          izisu=jisu+1
    !
    !
          if(nrmax.ge.izisu) go to 4000
          write(6,400) nrmax,jisu
      400 format(1h0,' nrmax,jisu='2i5,10x,'stop in suphod')
          stop
     4000 continue
    !
          nnr=int((xin-x(1))/(x(2)-x(1))+1)
    !
          if(nnr.ge.1.and.nnr.le.nrmax) go to 5000
          write(6,500) nnr
      500 format(1h0,' nnr='i5,10x,'stop in hokan --- xin < x(1) ---'&
                 ,1x,' or  --- xin > x(nrmax) ---')
          stop
     5000 continue
    !
          jmin=nnr-2
          jmax=nnr+1
          if(jmin.le.1) jmin=1
          if(jmax.ge.nrmax) jmin=nrmax-3
    !
    !
          x0=x(jmin)
          x1=x(jmin+1)
          x2=x(jmin+2)
          x3=x(jmin+3)
    !
          y0=y(jmin)
          y1=y(jmin+1)
          y2=y(jmin+2)
          y3=y(jmin+3)
    !
          xin0=xin-x0
          xin1=xin-x1
          xin2=xin-x2
          xin3=xin-x3
    !
          x01=x0-x1
          x02=x0-x2
          x03=x0-x3
    !
          x10=x1-x0
          x12=x1-x2
          x13=x1-x3
    !
          x20=x2-x0
          x21=x2-x1
          x23=x2-x3
    !
          x30=x3-x0
          x31=x3-x1
          x32=x3-x2
    !
          xx0=xin1/x01*xin2/x02*xin3/x03
    !
          xx1=xin0/x10*xin2/x12*xin3/x13
    !
          xx2=xin0/x20*xin1/x21*xin3/x23
    !
          xx3=xin0/x30*xin1/x31*xin2/x32
    !
          suphod=xx0*y0+xx1*y1+xx2*y2+xx3*y3
    !
          return
          end
    !
    !
    !==============================================
          subroutine absx(lmin,lmax,ldel,fk,s,ecm,absout)
    !==============================================
    !
          implicit real*8 (a-h,o-z)
          complex*16 s(120)
          integer absout
    !
          pai=3.141592653579893d0
    !
          abxsec=0.d0
          do l=lmin,lmax,ldel
          abxsec=abxsec+dble(2*l+1)*( 1.d0-abs(s(l+1))**2 )
    !      abxsec=abxsec+( 1.d0-cdabs(s(l+1))**2 )
          end do
    !
          abxsec=abxsec*pai/fk**2 *10.d0 ! in a unit of mb
    !
          !write(absout,*) ecm,abxsec
    !
          return
          end
    
