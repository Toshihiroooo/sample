module scattering_ws_mod

    implicit real*8(a-h,o-x)
    implicit character*256(y)
    implicit complex*16(z)

    contains

!===================================================================================================================

    subroutine variational_og(theta,N,bmin,bmax,l,s,ZE,Z_V)

        IMPLICIT real*8(a-h,o-x)
        IMPLICIT character*256(y)
        IMPLICIT complex*16(z)
        complex,parameter::zi=(0.0,1.d0),zr=(1.d0,0.0)

        ALLOCATABLE a(:),b(:),Aij(:,:),ZAij(:,:),ZTij(:,:),ZVij(:,:),ZHij(:,:),ZA(:),ZVRSUM(:),ZNOR(:),zrr(:)
        ALLOCATABLE ZW(:),ZWORK(:),ZALPHA(:),ZBETA(:),RWORK(:),ZVR(:,:),ZVL(:,:),ZE(:),Z_VR(:,:),Z_V(:,:)
        allocatable rgauss(:),amp(:)

     !------- setting parameter --------------------------------------------------------------------------------
        z=exp(zi*theta)
        dj=dble(l)+s
        sl=1.d0/2.d0*(dj*(dj+1)-dble(l)*dble(l+1)-s*(s+1.d0))
        amcc=931.4941d0
        hc=197.3271d0
        A1=40.d0
        A2=1.d0
        aucc=amcc*(A1*A2/(A1+A2))

        ND=50

        !print*,theta,N,bmin,bmax,l,zi
        !print*,s

     !------ Gaussian base expansion ----------------------------
        open(ND,file="checkin")

        read(ND,*) y_file,irmax,dr

        close(ND)
        ND=ND+1

        open(ND,file=y_file)

        read(ND,'(33x,I1)')lam
        read(ND,'(25x,I2)')nn

        nmax=nn

        allocate(rgauss(nn),amp(nn))

        read(ND,*) (rgauss(i),i=1,nmax)
        read(ND,'()')
        read(ND,*) (amp(i),i=1,nmax)

        close(ND)
        ND=ND+1

     !--------------------------------------------------------------------------------------------------------
        ALLOCATE(a(N),b(N),Aij(N,N),ZAij(N,N),ZTij(N,N),ZVij(N,N),ZHij(N,N),ZA(N),ZVRSUM(N),ZNOR(N),zrr(nn))
        
     !------- potential parameter ----------------------------------------------------------------------------
        zv=-42.d0*(1.d0+0.03d0*zi)

        do i=1,nmax
            rgauss(i)=1.d0/rgauss(i)/rgauss(i)
        end do
        
     !-------- initialize an array -----------------------------------------------------------------------------
        a=0.d0
        b=0.d0
        ZAij=(0.d0,0.d0)
        ZTij=(0.d0,0.d0)
        ZVij=(0.d0,0.d0)
        ZHij=(0.d0,0.d0)
        ZVRSUM=(0.d0,0.d0)

     !-------- setting width -------------------------------------------------------------------------------------
        DO i=1,N

            b(i)=bmin*(bmax/bmin)**(DBLE(i-1)/DBLE(N-1))
            a(i)=1.d0/(2.d0*b(i)**2.d0)

        END DO

     !-------- calculation of matrix element ----------------------------------------------------------------------
        !ND=50
        !open(ND,file='matrix',status='replace')
        DO i=1,N

            DO j=1,N
            !----- norm -----------------------------------------------------------------------------------------    
                ZAij(i,j)=(4.d0*a(i)*a(j)/(a(i)+a(j))**2.d0)**(0.5d0*(dble(l)+1.5d0))*zr
               
            !----- kinetic energy -------------------------------------------------------------------------------  
                ZTij(i,j)=0.5d0*hc*hc/aucc*(DBLE(l)+1.5d0)*(4.d0*a(i)*a(j)/(a(i)+a(j)))*ZAij(i,j)/(z**2.d0)
    
            !----- potential ------------------------------------------------------------------------------------
                do k=1,nmax
                    ZVij(i,j)=ZVij(i,j)+&
                    zv*amp(k)*((4.d0*a(i)*a(j)/(a(i)+a(j)+rgauss(k)*z**2.d0)**2.d0)**(dble(l)+1.5d0))**0.5d0
                end do
                
            !----- hamiltonian ----------------------------------------------------------------------------------
                ZHij(i,j)=ZTij(i,j)+ZVij(i,j)

                !print*,i,j,ZAij(i,j)

                !write(ND,*)i,j,ZVij(i,j)
    
            END DO
    
        END DO

        do i=1,N
            !print*,i,ZVij(i,1)
        end do

        !close(ND)
        !ND=ND+1
     !-------- setting lapack --------------------------------------------------------------------------------------
        NWR=8*N
        imax=N
        LDA=N
        LDB=N
        LDVL=N
        LDVR=N
        LWORK=33*N

        ALLOCATE(ZW(N),ZWORK(LWORK),ZALPHA(N),ZBETA(N),RWORK(NWR),ZVR(N,N),ZE(N),Z_VR(N,N),Z_V(N,N),ZVL(N,N))
           
     !------- using lapack ---------------------------------------------------------------------------------------
        CALL ZGGEV("V","V",N,ZHij,LDA,ZAij,LDB,ZALPHA,ZBETA,ZVL,LDVL,ZVR,LDVR,ZWORK,LWORK,RWORK,INFO)
        PRINT *,'INFO=',INFO
        
     !------- eigenenergy ------------------------------------------------------------------------------------------
        DO i=1,N
            ZW(i)=ZALPHA(i)/ZBETA(i)
            !print*,i,ZW(i)
            !print*,i,ZVR(1,i)
        END DO

        nm=N
        call bubble(ZW,ZE,ZVR,Z_VR,nm)

        do i=1,N
            !print*,i,ZE(i)
            !print*,i,Z_VR(1,i)
        end do
        
     !------- normalization -----------------------------------------------------------------------------------------
        !ND=50
        !open(ND,file='vector',status='replace')
        DO j=1,N
            DO i=1,N
                do k=1,N
                    ZAij(i,k)=(4.d0*a(i)*a(k)/(a(i)+a(k))**2.d0)**(0.5d0*(DBLE(l)+1.5d0))*zr
                    ZVRSUM(j)=ZVRSUM(j)+(Z_VR(i,j)*Z_VR(k,j)*ZAij(i,k))
                END DO
            END DO
        
              ZA(j)=(ZVRSUM(j))**0.5d0
              !print*, j,ZA(j)
     !------- coefficient ------------------------------------------------------------------------------------------
            DO i=1,N
                Z_V(i,j)=Z_VR(i,j)/ZA(j)
                !print*, i,Z_V(i,1)

                !write(ND,*)i,j,Z_V(i,j)

            END DO
            !print*, j,Z_V(j,1)
        END DO

        !close(ND)
        !ND=ND+1

        do j=1,N
            do i=1,N
                do k=1,N
                    !ZNOR(j)=ZNOR(j)+Z_V(i,j)*Z_V(k,j)*ZAij(i,k)
                end do
            end do
        
            !print*,j,ZNOR(j)
        end do

        return

    end subroutine variational_og

!===================================================================================================================

    subroutine bubble(ZW,Z,ZV,Z_V,N)

        IMPLICIT real*8(a-h,o-x)
        IMPLICIT character*256(y)
        IMPLICIT complex*16(z)

        ALLOCATABLE Z(:),ZW(:),ZV(:,:),Z_V(:,:)
        !ALLOCATE(Z(N),Z_V(N,N))

        DO i=1,N
            Z(i)=ZW(i)
            DO j=1,N
                Z_V(i,j)=ZV(i,j)
            END DO
        END DO

        DO i=1,N
            DO j=i,N
                if(real(Z(i))>real(Z(j))) then
                zbox=Z(i)
                Z(i)=Z(j)
                Z(j)=zbox

                DO k=1,N
                    zvrbox=Z_V(k,i)
                    Z_V(k,i)=Z_V(k,j)
                    Z_V(k,j)=zvrbox
                END DO

                end if
            END DO
        END DO

        return

    end subroutine bubble

!===================================================================================================================

    subroutine fn_gamma(x,f)

        IMPLICIT real*8(a-h,o-y)

        pi=4.d0*datan(1.d0)
        nx=int(x)
    
        !　--- half integer
        if(x-dble(nx).eq.0.5d0)then
    
            m=int(x-0.5d0)
            ff=1.d0
    
            if(m.eq.0)then
                f=pi**0.5d0
            else
                do i=1,2*m-1,+2
                    ff=ff*dble(i)
                end do
    
                f=ff*(pi**0.5d0)/2.d0**dble(m)
    
            end if
    
        ! --- integer
        else if(x-dble(nx).eq.0.d0)then
            f=1.d0
    
            if(x.eq.0.d0)then
                f=1.d0
            else
                do i=1,nx-1
                    f=f*dble(i)
                end do
            end if
    
        end if

        return
    
    end subroutine fn_gamma

!====================================================================================================================
!====================================================================================================================

    end module scattering_ws_mod

!====================================================================================================================
!====================================================================================================================

program scattering

    use scattering_ws_mod

    implicit real*8(a-h,o-x)
    implicit character*256(y)
    implicit complex*16(z)

    allocatable a(:),b(:),cn(:),zbl(:),ZE(:),Z_V(:,:),zd(:)
    allocatable rgauss(:),amp(:),zrr(:),zx(:),zbj0(:),zbj1(:)

 !------- setting parameter ---------------------------
    zi=(0.d0,1.d0)
    zr=(1.d0,0.d0)

    !print*,"input angular CSM : degree"
    !read*,deg

    !print*,"input l"
    !read*,l

    !print*,"input spin"
    !read*,s

    !print*,"input angle of scattering"
    !read*,deg2

    deg=15.d0
    l=0
    s=0.5d0
    deg2=0.d0

    print*,"================================================"
    print*,"CSM angular is ",deg
    print*,"quantum number l is ",l
    print*,"qyantum number s is ",s
    print*,"scattering angular is ",deg2
    print*,"================================================"

    N=40
    !s=1.d0/2.d0
    dj=s+dble(l)
    amcc=931.4941d0
    hc=197.3271d0
    A1=40.d0
    A2=1.d0
    aucc=amcc*(A1*A2/(A1+A2))

    ND=50
    NE=30
    NF=40

 !------ Gaussian base expansion ----------------------------
    open(ND,file="checkin")

    read(ND,*) y_file,irmax,dr

    close(ND)
    ND=ND+1

    open(ND,file=y_file)

    read(ND,'(33x,I1)')lam
    read(ND,'(25x,I2)')nn

    nmax=nn

    allocate(rgauss(nn),amp(nn))

    read(ND,*) (rgauss(i),i=1,nmax)
    read(ND,'()')
    read(ND,*) (amp(i),i=1,nmax)

    close(ND)
    ND=ND+1


    allocate(a(N),b(N),cn(N),zbl(N),zd(N),zrr(nn))
    allocate(zx(nn),zbj0(nn),zbj1(nn))

 !------- setting parameter ---------------------------
    pi=4.d0*datan(1.d0)
    theta=deg*pi/180.d0
    theta2=deg2*pi/180.d0
    z=exp(theta*zi)

    bmin=0.1d0
    bmax=41.d0

 !------- setting width --------------------
    do i=1,N
        b(i)=bmin*(bmax/bmin)**(dble(i-1)/dble(N-1))
        a(i)=1.d0/(2.d0*b(i)**2.d0)
        !print*,b(i)
    end do

 !------- nomalization ----------------
    call fn_gamma(dble(l+1.5d0),fg)

    do i=1,N
        cn(i)=sqrt(2.d0*(2.d0*a(i))**(dble(l)+1.5d0)/fg)
        !print*,cn(i)
    end do

 !------ potential parameter ---------------
    zv=-42.d0*(1.d0+0.03d0*zi)

        do i=1,nmax
            rgauss(i)=1.d0/rgauss(i)/rgauss(i)
        end do

 !------ caluclation b ---------------

    call variational_og(theta,N,bmin,bmax,l,s,ZE,Z_V)

    open(ND,file='energy_sc',status='replace')


    do i=1,N
        print*,i,real(ZE(i)),aimag(ZE(i))
        write(ND,*)real(ZE(i)),aimag(ZE(i))
    end do

    close(ND)
    ND=ND+1

    sl=0.5d0*(dj*(dj+1.d0)-dble(l)*dble(l+1)-s*(s+1.d0))
    e=15.d0
    es=0.d0!10.d0
    imax=5000
    de=e/dble(imax)

    write(y_l,"(i1)")l

    if(s==0.5) then

        y_s='u'

    elseif(s==-0.5) then

        y_s='d'
        
    end if

    open(ND,file='ap_c_l='//trim(y_l)//'_s='//trim(y_s),status='replace')
    !open(NE,file='ps_c_l='//trim(y_l)//'_s='//trim(y_s),status='replace')
    open(NE,file="ld_and_trance",status="replace")
    open(NF,file='width_l='//trim(y_l)//'_s='//trim(y_s),status='replace')

  do k=1,imax

    ee=es+dble(k)*de

    dk=sqrt(2.d0*aucc/(hc**2.d0)*ee)

    !print*,dk

    zbl=(0.d0,0.d0)

    !print*,ee
    do i=1,N

        do j=1,nmax

            zbl(i)=zbl(i)+cn(i)*dk*(z**1.5d0)*(sqrt(pi)/4.d0)*((z*0.5d0*dk)**dble(l))*&
            exp((-1.d0)*z**2.d0/(4.d0*(a(i)+rgauss(j)*(z**2.d0)))*dk**2.d0)*&
            zv*amp(j)/((a(i)+rgauss(j)*(z**2.d0))**(dble(l)+1.5d0))

        end do

        !print*,i,zbl(i)
    end do


 !------ calculation d -----------------------------

    zd=(0.d0,0.d0)
    do j=1,N

        do i=1,N

            zd(j)=zd(j)+exp((-1.d0)*zi*theta/2.d0)*zbl(i)*Z_V(i,j)

        end do
        !print*,j,zd(j)

    end do

 !------ scattering amplitude ---------------------
    zfsl=(0.d0,0.d0)

    do i=1,N

        zfsl=zfsl+(-1.d0)*z*2.d0*aucc/hc/hc/dk/dk*zd(i)*zd(i)/(ee-ZE(i))

        !zfsl=zfsl+(-1.d0)*z*2.d0*aucc/hc/hc/dk/dk*zd(i)*zd(i)!/(ee-ZE(i))

        !zfsl=zfsl+(-1.d0)*z*2.d0*aucc/hc/hc/dk/dk/(ee-ZE(i))

    end do
    !print*,zfsl

 !----- borun term -------------

    do i=1,nn
        zx(i)=zi*(dk**2.d0)/(2.d0*rgauss(i))
    end do

    do i=1,nn
        zbj0(i)=sin(zx(i))/zx(i)
        zbj1(i)=(-1.d0)*(zx(i)*cos(zx(i))-sin(zx(i)))/(zx(i)*zx(i))
    end do
 !------ level density -----------

    U=8.362d0+ee
    delta=12.8d0*(A1+A2)**(-0.5d0)
    di=dble(l)+s
    aa=(A1+A2)/5.5d0
    t=(-1.d0+(1.d0+4.d0*(U+delta))**0.5d0)/(2.d0*aa)
    sigma=0.015d0*((A1+A2)**(5.d0/3.d0))*t

    !print*,aa,U,delta

    rho=1.d0/(24.d0*(2.d0**0.5d0))*(2.d0*di+1.d0)/((sigma**3.d0)*(aa**0.25d0))*&
    exp(2.d0*(aa*(U-delta))**0.5d0-di*(di+1.d0)/(2.d0*(sigma**2.d0)))/&
    (U-delta+t)**(5.d0/4.d0)

    !rho=2.d0*(aa*(U-delta))!**0.5d0!-di*(di+1.d0)/(2.d0*(sigma**2.d0))
    !print*,rho

  !-------------------------------------------------------------------
    if(l==0) then

    zfb0=(0.d0,0.d0)
    do i=1,nn
        zfb0=zfb0+zv*(dk**2.d0)*(amp(i)*(2.d0*rgauss(i))**(-1.5d0))*&
        exp(-2.d0*(dk**2.d0)/(4.d0*rgauss(i)))*zbj0(i)
    end do

        zfb0=-2.d0*aucc/(hc**2.d0*dk**2.d0)*sqrt(pi/2.d0)*zfb0
        !print*,zfb0

        zfs0=zfb0+zfsl
        !zfs0=zfsl
        !zfs0=zfb0

        write(ND,*) ee,real(zfs0),aimag(zfs0),real(zfs0)**2.d0+aimag(zfs0)**2.d0

        za0=2.d0*zi*dk*zfs0
        !print*,za0

        write(NE,*) ee,rho,(1.d0-(real(1.d0+za0)**2.d0+aimag(1.d0+za0)**2.d0))
        write(NF,*) ee,(1.d0-(real(1.d0+za0)**2.d0+aimag(1.d0+za0)**2.d0))/(pi*2.d0*rho)*1000000.d0

  !------------------------------------------------------------------
    elseif(l==1) then
    zfb1=(0.d0,0.d0)

    do i=1,nn
        zfb1=zfb1+zv*(dk**2.d0)*(amp(i)*(2.d0*rgauss(i))**(-1.5d0))*&
        exp(-2.d0*(dk**2.d0)/(4.d0)*rgauss(i))*(-zi)**dble(l)*zbj1(i)
    end do

        zfb1=-2.d0*aucc/(hc**2.d0*dk**2.d0)*sqrt(pi/2.d0)*zfb1

        zfs1=(zfb1+zfsl)*cos(theta2)

        write(ND,*) ee,real(zfs1),aimag(zfs1),real(zfs1)**2.d0+aimag(zfs1)**2.d0

        za1=(2.d0*zi*dk)/(cos(theta2))*zfs1

        write(NE,*) ee,real(1.d0+za1),aimag(1.d0+za1)

    end if


  end do

close(ND)
close(NE)
close(NF)
ND=ND+1
NE=NE+1
NF=NF+1

end program
