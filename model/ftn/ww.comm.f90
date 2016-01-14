      MODULE WW_cc

      USE CMP_COMM, ONLY:   &

         MPI_COMM_WW => COMM_local,   &

         Coupler_id,     &
         component_master_rank_local,     &
         process_rank_local,     &
         component_nprocs,     &
         ibuffer,     &

         MPI_INTEGER,MPI_STATUS_SIZE,     &

         kind_REAL,kind_alt_REAL,MPI_kind_REAL,MPI_kind_alt_REAL

      implicit none

!     integer,parameter:: kind_R=kind_REAL
!<-: if WW uses real*8; if it uses real*4 :->
      integer,parameter:: kind_R=kind_alt_REAL

      integer,parameter:: kind_SBC=kind_R
      integer MPI_kind_R,MPI_kind_SBC

!CBT
      integer,parameter:: kind_CURR=kind_R
      integer MPI_kind_CURR
      integer,parameter:: kind_wst=kind_R
      integer MPI_kind_wst
      integer,parameter:: kind_kpph=kind_R
      integer MPI_kind_kpph
!CBT.

      integer NX,NY
      real(kind=kind_R) dtc /1./       !<- Coupling period
!     real(kind=kind_R) dt                !<- WW time step
!     integer i_dtc2dtw /100/,  !<- Coupling period / WW time step
      integer n_ts /0/          !<- number of time step
      real(kind=kind_R) x0,y0,sx,sy
      real(kind=kind_R) time_s

      integer NGP,rc
      logical getPSEUDOICE,getSBC,getCUR,getKPPH
      integer itime(2)

      character*100 ch

!Controls:
!      logical dummy_C_mode /.false./
!      logical dummy_dummy_C_mode /.false./
     logical dummy_C_mode /.true./
     logical dummy_dummy_C_mode /.true./
      integer nunit_announce /6/, VerbLev /3/

      SAVE

      END MODULE WW_cc
!
!***********************************************************************
!
      MODULE SURF_cc

      USE WW_cc, ONLY: kind_SBC

      implicit none

      integer num_SBC /2/  ! # of surf. boundary cond. to be dealt with
                        ! num_SBC=2 means that only the 2 surf. wind
                        ! components are supplied

      integer,allocatable:: ice(:,:)    ! "pseudoice"
      real(kind=kind_SBC),dimension(:,:,:),target,allocatable:: SBC,SBC1
!CBT
      real(kind=kind_SBC),dimension(:,:), allocatable:: RIBN, ZLML
!CBT
      SAVE

      END MODULE SURF_cc
!
!***********************************************************************
!
!CBT
      MODULE supl_vars_def

      USE WW_cc, ONLY: kind_CURR, kind_wst, kind_KPPH

      IMPLICIT NONE

      real(kind=kind_CURR),dimension(:,:),target,allocatable:: SFCURX, &
      SFCURY, DPCURX, DPCURY 
      real(kind=kind_CURR),dimension(:,:),target,allocatable:: KPPML
      real(kind=kind_wst),dimension(:,:),allocatable::dtaux, dtauy
      real(kind=kind_wst),dimension(:,:),allocatable::alpha, gamma,  &
                                                      lamda, wbcond
      real(kind=kind_wst),dimension(:,:),allocatable::stkdfx, stkdfy
      SAVE

      END MODULE supl_vars_def

!CBT.
!
!***********************************************************************
!
      SUBROUTINE WW_CMP_START

      USE WW_cc

      implicit none

      integer WW_id /3/, WW_master_rank_local /0/
      character*20 s
!

                      !<-id of WW as a component of the coupled system
      call CMP_INIT(WW_id,1)
                          !<-"flexibility level"
      if (Coupler_id.ge.0) then
        dummy_C_mode=.false.
        VerbLev=min(VerbLev,ibuffer(4))
      end if

      dummy_dummy_C_mode=dummy_dummy_C_mode.and.dummy_C_mode

      write(s,'(2L1)') dummy_C_mode,dummy_dummy_C_mode
      call WW_ANNOUNCE('dummy_C_mode='//s(1:1)//    &
      '  dummy_dummy_C_mode='//s(2:2),1)

      call CMP_INTRO(WW_master_rank_local)

      write(s,'(i2)') VerbLev
      call WW_ANNOUNCE('back from CMP_INTRO, VerbLev='//s,2)

      if (kind_R.eq.kind_REAL) then
        MPI_kind_R=MPI_kind_REAL
      else 
        MPI_kind_R=MPI_kind_alt_REAL
      end if
      if (kind_SBC.eq.kind_REAL) then
        MPI_kind_SBC=MPI_kind_REAL
      else 
        MPI_kind_SBC=MPI_kind_alt_REAL
      end if

!CBT
      if (kind_CURR.eq.kind_REAL) then
        MPI_kind_CURR=MPI_kind_REAL
      else
        MPI_kind_CURR=MPI_kind_alt_REAL
      end if
      if (kind_wst .eq. kind_REAL) THEN
        MPI_kind_wst = MPI_kind_REAL
      else
        MPI_kind_wst = MPI_kind_alt_REAL
      endif
      if (kind_kpph .eq. kind_REAL) THEN
        MPI_kind_kpph = MPI_kind_REAL
      else
        MPI_kind_kpph = MPI_kind_alt_REAL
      endif
!CBT.

      call WW_ANNOUNCE('WW_CMP_START: returning',2)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_INIT(imod,nx_,ny_,x0_,y0_,sx_,sy_)

      USE WW_cc

      USE SURF_cc
      USE supl_vars_def

      implicit none

      integer imod
      integer nx_,ny_
      real(kind=kind_R) x0_,y0_,sx_,sy_

      logical first /.true./
      character*80 s
      SAVE
!

!     if (VerbLev.ge.3) print*,'WW_INIT entered, args: ',
!    >imod,nx_,ny_,x0_,y0_,sx_,sy_
      write (s,'(i4,2i5,1p,4e15.7)') imod,nx_,ny_,x0_,y0_,sx_,sy_
      call WW_ANNOUNCE('WW_INIT entered, args: '//s,3)

      if (.not.first) then
        call GLOB_ABORT(1,'WW_INIT called more than once: aborted',1)
      end if
      first=.false.

      NX=nx_
      NY=ny_
      x0=x0_
      y0=y0_
      sx=sx_
      sy=sy_

      NGP=NY*NX

      allocate(ice(NX,NY))
      allocate(SBC(NX,NY,num_SBC))
!CBT
      allocate(RIBN(NX,NY), ZLML(NX,NY))
!CBT
      if (dummy_dummy_C_mode)  allocate(SBC1(NX,NY,num_SBC))
!CBT
      allocate(SFCURX(NX,NY),SFCURY(NX,NY))
      allocate(DPCURX(NX,NY),DPCURY(NX,NY))
      allocate(DTAUX(NX,NY),DTAUY(NX,NY))
      allocate(ALPHA(NX,NY),GAMMA(NX,NY),LAMDA(NX,NY))
      allocate(WBCOND(NX,NY))
      allocate(KPPML(NX,NY))
      allocate(STKDFX(NX,NY),STKDFY(NX,NY))
!CBT.

      write (s,'(i4,2i5,1p,4e15.7)') imod,nx_,ny_,x0_,y0_,sx_,sy_
      call WW_ANNOUNCE('WW_INIT returning, args: '//s,2)
!     if (VerbLev.ge.2) print*,'WW_INIT returning, args: ',
!    >imod,nx_,ny_,x0_,y0_,sx_,sy_

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_RECVDTC

      USE WW_cc

      implicit none

      character*20 s
!

      if (dummy_C_mode) then
        dtc=600.
        RETURN
      end if

      if (Coupler_id.lt.0) RETURN    !   <- standalone mode

      IF (component_master_rank_local.eq.process_rank_local) THEN
        call CMP_RECV(dtc,1)
        write(s,'(1pe20.12)') dtc
        call WW_ANNOUNCE('coupling period received: '//s,1)
      END IF

      call MPI_BCAST(dtc,1,MPI_kind_REAL,    &
      component_master_rank_local,MPI_COMM_WW,rc)

      call WW_ANNOUNCE('(BP) dtc='//trim(s)//' broadcast OK',2)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_SENDGRIDS

      USE WW_cc

      implicit none

      real(kind=kind_REAL),dimension(NX,NY):: ALON,ALAT
      integer i,j,ibuf(2)
      character*20 s
!

      if (Coupler_id.lt.0) RETURN    !   <- standalone mode

      IF (component_master_rank_local.eq.process_rank_local) THEN
        ibuf(1)=NX
        ibuf(2)=NY
        call WW_ANNOUNCE('to send grid dimensions',2)
        call CMP_INTEGER_SEND(ibuf,2)
        call WW_ANNOUNCE('grid dimensions sent',1)
      END IF

      do i=1,NX
        ALON(i,:)=x0+(i-1)*sx
      end do
      do j=1,NY
        ALAT(:,j)=y0+(j-1)*sy
      end do

      call WW_ANNOUNCE('(BP) to send grid arrays',2)

      call CMP_SEND(ALON,NGP)
      call CMP_SEND(ALAT,NGP)

      call WW_ANNOUNCE('the 2 grid arrays sent',1)

      call WW_ANNOUNCE('(BP) WW_SENDGRIDS: returning',2)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_SENDSLM(ma)

      USE WW_cc

      implicit none

      integer,dimension(NY,NX):: ma  ! NB! why?...

      real(kind=kind_R) FG(NX,NY)
      integer i,j
!

      if (Coupler_id.lt.0) RETURN    !   <- standalone mode

      do i=1,NX
      do j=1,NY
        if (ma(j,i).gt.0) then  !<- "active sea points" (value 1) in WW
          FG(i,j)=1.
        else           !<- according to Arun: 0 is a land point;
          FG(i,j)=0.   ! negative values are: points covered by ice;
        end if         ! points that became dry; land points for a
      end do           ! relocatable grid
      end do

      call CMP_gnr_SEND(FG,NGP,MPI_kind_R)

      call WW_ANNOUNCE('Sea/land mask sent',1)

      return
      END
!
!***********************************************************************
!
!CBT
      SUBROUTINE WW_SENDWST(dtau_x,dtau_y,chrnk,msang,mwlen,     &
                            bcondw,stk_x,stk_y)
      USE WW_cc
      USE supl_vars_def
      IMPLICIT NONE
      INTEGER :: i, j
      REAL(KIND=kind_wst) dtau_x(ny,nx), dtau_y(ny,nx), &
          chrnk(ny,nx), msang(ny,nx), mwlen(ny,nx),     &
          bcondw(ny,nx), stk_x(ny,nx),stk_y(ny,nx)

      if (Coupler_id.lt.0) RETURN    !   <- standalone mode
      
      DO i = 1, nx
      DO j = 1, ny
         dtaux(i,j) = dtau_x(j,i)
         dtauy(i,j) = dtau_y(j,i)
         alpha(i,j) = chrnk(j,i)
         gamma(i,j) = msang(j,i)
         lamda(i,j) = mwlen(j,i)
         wbcond(i,j) = bcondw(j,i) 
         stkdfx(i,j) = stk_x(j,i)
         stkdfy(i,j) = stk_y(j,i)
      ENDDO
      ENDDO
      call CMP_gnr_SEND(dtaux,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(dtauy,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(stkdfx,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(stkdfy,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(alpha,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(gamma,NGP,MPI_kind_wst) 
      call CMP_gnr_SEND(wbcond,NGP,MPI_kind_wst)
      call CMP_gnr_SEND(lamda,NGP,MPI_kind_wst)

      call WW_ANNOUNCE('Delta(Tau) sent',1)

      RETURN
      END 
!CBT.
!
!***********************************************************************
!

      SUBROUTINE WW_ANNOUNCE(s,DbgLev)

      USE WW_cc, ONLY: nunit_announce,VerbLev,MPI_COMM_WW

      implicit none

      character*(*) s
      integer DbgLev

      integer ierr
!
      if (DbgLev.le.VerbLev) then
        if (s(1:5).eq.'(BP) ') then
          call MPI_BARRIER(MPI_COMM_WW,ierr)
        end if
        CALL CMP_ANNOUNCE(nunit_announce,'WW: '//s)
      end if

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_TSTEP_INIT(itime_)

      USE WW_cc

      USE W3TIMEMD, ONLY: DSEC21,TICK21

      implicit none

      integer itime_(2)

      integer init_time(2),ntime_dtc

      SAVE
!

      n_ts=n_ts+1
      itime=itime_
      if (n_ts.eq.1) init_time=itime 
      time_s=DSEC21(init_time,itime)
      ntime_dtc=(time_s+0.1)/dtc

      getPSEUDOICE=time_s-ntime_dtc*dtc.lt.0.1
      getSBC=getPSEUDOICE
      getCUR=getPSEUDOICE
      getKPPH=getPSEUDOICE

!zz      write(ch,'(i2," n_ts=",i6," time_s=",f10.1," ntime_dtc=",i6,   &
!zz      " itime: ",2i8," getPSEUDOICE=",L1," getSBC=",L1)')   &
!zz      process_rank_local,n_ts,time_s,ntime_dtc,itime,getPSEUDOICE,getSBC
      write(ch,'(i2," n_ts=",i6," time_s=",f10.1," ntime_dtc=",i6)')   &
      process_rank_local,n_ts,time_s,ntime_dtc
      write(ch, '(" itime: ",2i8," getPSEUDOICE=",L1," getSBC=",L1)')   &
      itime,getPSEUDOICE,getSBC

      call WW_ANNOUNCE('WW_TSTEP_INIT: '//trim(ch),3)
!     if (VerbLev.ge.3 .and. n_ts.eq.2) then
!       print*,'WW_TSTEP_INIT: '//trim(ch)
!     end if

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_CDETECT(cf)

      USE WW_cc, ONLY: Coupler_id, dummy_C_mode

      implicit none

      logical cf
!

      cf= Coupler_id.ge.0 .or. dummy_C_mode

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_RECV_SBC(RIB, ZBL)

      USE WW_cc

      USE SURF_cc

      implicit none
      REAL(kind=kind_SBC), DIMENSION(ny,nx)   :: RIB, ZBL
      REAL, PARAMETER :: RIB_L = 5.0_kind_SBC, z1_ll = 5.0_kind_SBC, &
                                               z1_ul = 50.0_kind_SBC

      integer i, j, n
!

      call WW_ANNOUNCE('WW_RECV_SBC entered '//trim(ch),3)

      if (.not. getSBC) RETURN

      if (dummy_C_mode) then
        if (.not. getPSEUDOICE) call GLOB_ABORT(1,   &
        'getPSEUDOICE not equiv. to getSBC in dummy mode, aborted',1)
        !call WW_DUMMY_SBC
!CBT
        CALL DUMMY_WND
!CBT
        RETURN
      end if

      if (Coupler_id.lt.0) RETURN     !   <- standalone mode

      if (getSBC) then
!CBT
        call CMP_gnr_RECV(RIBN(:,:),NGP,MPI_kind_SBC)
        call CMP_gnr_RECV(ZLML(:,:),NGP,MPI_kind_SBC)
!CBT
        do n=1,num_SBC
          call CMP_gnr_RECV(SBC(:,:,n),NGP,MPI_kind_SBC)
        end do
!CBT
        call MPI_BCAST(SBC,NGP*num_SBC,MPI_kind_SBC,    &
        component_master_rank_local,MPI_COMM_WW,rc)
!CBT
        call MPI_BCAST(RIBN,NGP,MPI_kind_SBC,    &
        component_master_rank_local,MPI_COMM_WW,rc)
        call MPI_BCAST(ZLML,NGP,MPI_kind_SBC,    &
        component_master_rank_local,MPI_COMM_WW,rc)
        DO j = 1, ny
        DO i = 1, nx
          IF ( ABS(RIBN(i,j)) .LT. RIB_L ) THEN
             RIB(j,i) = RIBN(i,j)
          ELSE
             RIB(j,i) = 0.0_kind_SBC
          ENDIF         
          ZBL(j,i) = MIN(z1_ul, MAX(z1_ll, ZLML(i,j)))
        ENDDO
        ENDDO 
!CBT
      end if

      return
      END
!
!***********************************************************************
!
!CBT
       SUBROUTINE WW_RECV_CUR(cdx,cdy)
       USE WW_cc
       USE supl_vars_def

       IMPLICIT NONE
       INTEGER :: i, j 
       REAL(kind=kind_CURR), DIMENSION(ny,nx) :: cdx, cdy
       REAL, PARAMETER :: sfcurmax = 5._kind_CURR                   
       
       call WW_ANNOUNCE('WW_RECV_CUR entered '//trim(ch),3) 
       
       if (.not. getCUR) RETURN 

       if (dummy_C_mode) then
        if (.not. getPSEUDOICE) call GLOB_ABORT(1,   &
        'getPSEUDOICE not equiv. to getCUR in dummy mode, aborted',1)
        call DUMMY_CUR
        RETURN
      end if

      if (Coupler_id.lt.0) RETURN     !   <- standalone mode

      if (getCUR) then
       call CMP_gnr_RECV(SFCURX,NGP,MPI_kind_CURR)
       call CMP_gnr_RECV(SFCURY,NGP,MPI_kind_CURR)
       call CMP_gnr_RECV(DPCURX,NGP,MPI_kind_CURR)
       call CMP_gnr_RECV(DPCURY,NGP,MPI_kind_CURR)
       call MPI_BCAST(SFCURX,NGP,MPI_kind_CURR,    &
       component_master_rank_local,MPI_COMM_WW,rc)
       call MPI_BCAST(SFCURY,NGP,MPI_kind_CURR,    &
       component_master_rank_local,MPI_COMM_WW,rc)
       call MPI_BCAST(DPCURX,NGP,MPI_kind_CURR,    &
       component_master_rank_local,MPI_COMM_WW,rc)
       call MPI_BCAST(DPCURY,NGP,MPI_kind_CURR,    &
       component_master_rank_local,MPI_COMM_WW,rc) 
       DO i = 1, nx
       DO j = 1, ny
         IF ( ABS(SFCURX(i,j)) .GT. sfcurmax ) THEN
            SFCURX(i,j) = 0.0_kind_CURR
         ENDIF
         IF ( ABS(SFCURY(i,j)) .GT. sfcurmax ) THEN
            SFCURY(i,j) = 0.0_kind_CURR
         ENDIF
         IF ( ABS(dpcurx(i,j)) .GT. sfcurmax ) THEN
            cdx(j,i) = 0.0_kind_CURR
         ELSE
            cdx(j,i) = dpcurx(i,j)
         ENDIF
         IF ( ABS(dpcury(i,j)) .GT. sfcurmax ) THEN
            cdy(j,i) = 0.0_kind_CURR
         ELSE
            cdy(j,i) = dpcury(i,j)
         ENDIF
       ENDDO
       ENDDO

      end if

      return
      END 
!CBT.

!
!***********************************************************************
!
!CBT

       SUBROUTINE WW_RECV_KPP(kpph)
       USE WW_cc
       USE supl_vars_def

       IMPLICIT NONE
       INTEGER :: i, j
       REAL(kind=kind_KPPH), DIMENSION(ny,nx) :: kpph
       REAL, PARAMETER :: kpphmax = 300._kind_KPPH,   &
                          kpphmin = 6._kind_KPPH

       call WW_ANNOUNCE('WW_RECV_CUR entered '//trim(ch),3)

       if (.not. getKPPH) RETURN

       if (Coupler_id.lt.0) RETURN     !   <- standalone mode

       if (getKPPH) then
       call CMP_gnr_RECV(KPPML,NGP,MPI_kind_KPPH)
       call MPI_BCAST(KPPML,NGP,MPI_kind_KPPH,    &
       component_master_rank_local,MPI_COMM_WW,rc)
       DO i = 1, nx
       DO j = 1, ny
        kpph(j,i) = max(min(kppml(i,j),kpphmax),kpphmin)
       ENDDO
       ENDDO

      end if

      return
      END

!CBT.

!
!***********************************************************************
!
      SUBROUTINE WW_DUMMY_SBC

      USE WW_cc

      USE SURF_cc

      implicit none

      real(kind=kind_R) XC,YC,XC0,YC0,SC,SDR,X,Y,UMX,RMX,RMX2,WSPD,  &
      DISTANCE,dt_dummy,t
      real(kind=kind_R),parameter:: DERA=3.1415927/180.
      integer i,j,k,k0
      real(kind=kind_R),dimension(:,:,:),pointer:: F
!

      call WW_ANNOUNCE('WW_DUMMY_SBC entered',3)

      dt_dummy=dtc

      XC0=300.
      YC0=25.
      SC=0.00005
      SDR=90.*DERA
      UMX=45.
      RMX=0.5
      RMX2=5.

      k0=0
      if (dummy_dummy_C_mode) k0=1
      F=>SBC
!
      DO k=0,k0
!
      if (k.eq.1) F=>SBC1
      t=time_s+k*dt_dummy
      XC=XC0-SC*SIN(SDR)*t
      YC=YC0-SC*COS(SDR)*t

      do j=1,NY
        do i=1,NX
          X=x0+(i-1)*sx-XC
          Y=y0+(j-1)*sy-YC
          DISTANCE=SQRT(X**2+Y**2)
          IF (DISTANCE.LT.0.01*RMX .OR. DISTANCE.GT.RMX2) THEN
            F(i,j,:)=0.
          ELSE
            IF (DISTANCE.LE.RMX) THEN
              WSPD=UMX*DISTANCE/RMX
            ELSE
              WSPD=UMX*RMX/DISTANCE
            END IF
            F(i,j,1)=-Y*WSPD/DISTANCE
            F(i,j,2)=X*WSPD/DISTANCE
          END IF
!         if (k.eq.0) then
          if (k.eq.k0) then
            IF(DISTANCE.LT.1.1*RMX2) THEN
              ice(i,j)=0
            ELSE
              ice(i,j)=1
            END IF
          end if
        end do
      end do
!
      END DO
!
      if (time_s.gt.0.1) RETURN
      if (.not.dummy_dummy_C_mode) SBC=0.
      RETURN
      END
!
!***********************************************************************
!
      SUBROUTINE WW_UPDATE_WIND(itime0,itime1,WX0,WY0,WX1,WY1)

      USE WW_cc

      USE SURF_cc

      USE W3TIMEMD, ONLY: DSEC21,TICK21

      implicit none

      integer itime0(2),itime1(2)
      real(kind=kind_SBC),dimension(NX,NY):: WX0,WY0,WX1,WY1
      real(kind=kind_R) dtime
      character*80 s
!

      if (Coupler_id.lt.0 .and. .not.dummy_C_mode) RETURN
                                !   <- standalone mode

      write (s,'(4i15)') itime0,itime1
      call WW_ANNOUNCE('WW_UPDATE_WIND: itime0,itime1 on entry: '//s,3)
!        print*,'WW_UPDATE_WIND: itime0,itime1 on entry: ',itime0,itime1
      itime0=itime
      itime1=itime0
      dtime=dtc
      call TICK21(itime1,dtime)

      WX0=SBC(:,:,1)
      WY0=SBC(:,:,2)
      if (dummy_dummy_C_mode) then
!        if (n_ts.gt.1) then  !debug
!<- commented out Dec.18 2008, should've earlier
        WX1=SBC1(:,:,1)
        WY1=SBC1(:,:,2)
!        end if   !debug
      else
        WX1=WX0
        WY1=WY0
      end if
      write (s,'(4i15)') itime0,itime1
      call WW_ANNOUNCE('WW_UPDATE_WIND: itime0,itime1 on exit: '//s,3)
!        print*,'WW_UPDATE_WIND: itime0,itime1 on exit: ',itime0,itime1

      return
      END

!
!***********************************************************************
!
!CBT
      SUBROUTINE WW_UPDATE_CUR(itime0,itime1,CX0,CY0,CX1,CY1)

      USE WW_cc
      USE supl_vars_def
      USE W3TIMEMD, ONLY: DSEC21,TICK21

      implicit none

      integer itime0(2),itime1(2)
      real(kind=kind_CURR),dimension(NX,NY):: CX0,CY0,CX1,CY1
      real(kind=kind_R) dtime
      character*80 s

       
      if (Coupler_id.lt.0 .and. .not.dummy_C_mode) RETURN
                                !   <- standalone mode

      write (s,'(4i15)') itime0,itime1
      call WW_ANNOUNCE('WW_UPDATE_CUR: itime0,itime1 on entry: '//s,3)
      itime0=itime
      itime1=itime0
      dtime=dtc
      call TICK21(itime1,dtime)

      CX0=SFCURX
      CY0=SFCURY
  
      if (dummy_dummy_C_mode) then
!        if (n_ts.gt.1) then  !debug
!<- commented out Dec.18 2008, should've earlier
        CX1=0.
        CY1=0.
!        end if   !debug
      else
        CX1=CX0
        CY1=CY0
      end if
      write (s,'(4i15)') itime0,itime1
      call WW_ANNOUNCE('WW_UPDATE_CUR: itime0,itime1 on exit: '//s,3)
!        print*,'WW_UPDATE_CUR: itime0,itime1 on exit: ',itime0,itime1


      return
      END
!CBT.
!
!***********************************************************************
!
      SUBROUTINE WW_RECV_PSEUDOICE

      USE WW_cc

      USE SURF_cc

      implicit none
!

      call WW_ANNOUNCE('WW_RECV_PSEUDOICE entered '//ch,3)
 
      if (.not. getPSEUDOICE) RETURN

      if (dummy_C_mode) then
        call WW_ANNOUNCE('WW_RECV_PSEUDOICE: dummy mode: ice(:,:) '//     &
        'assigned in WW_RECV_SBC. Returning',3)
        RETURN
      end if

      if (Coupler_id.lt.0) RETURN     !   <- standalone mode

      call CMP_gnr_RECV(ice,NGP,MPI_INTEGER)
                                       ! this array's values must be 1
                                       ! for "ice" and 0 for "open sea"
                                       ! - potential computational GPs
                                       ! that may be either actually
                                       ! computational (open sea) in WW
                                       ! or WW land or WW ice
!*** However, currently NO real WW ice is supposed to be there! ***

      call MPI_BCAST(ice,NGP,MPI_INTEGER,   &
      component_master_rank_local,MPI_COMM_WW,rc)
        
      call WW_ANNOUNCE('WW_RECV_PSEUDOICE to return',3)

      return
      END
!
!***********************************************************************
!
      SUBROUTINE WW_UPDATE_PSEUDOICE(itime2,rice)

      USE WW_cc

      USE SURF_cc

      USE W3TIMEMD, ONLY: DSEC21,TICK21

      implicit none

      integer itime2(2)
      real(kind=kind_R) rice(NX,NY)

      real(kind=kind_R) dtime
      character*30 s
!

      if (Coupler_id.lt.0 .and. .not.dummy_C_mode) RETURN
                                !   <- standalone mode

      write (s,'(2i12)') itime2
      call WW_ANNOUNCE('WW_UPDATE_PSEUDOICE: itime2 on entry: '//s,3)
!        print*,'WW_UPDATE_PSEUDOICE: itime2 on entry: ',itime2
      itime2=itime
      dtime=dtc
      call TICK21(itime2,dtime)

!        if (n_ts.gt.1) then   !debug
      rice=ice
!        end if   !debug
      write (s,'(2i12)') itime2
      call WW_ANNOUNCE('WW_UPDATE_PSEUDOICE: itime2 on exit: '//s,3)

      return
      END
!
!***********************************************************************
!
SUBROUTINE DUMMY_WND
  USE WW_cc
  USE SURF_cc
  IMPLICIT NONE
  INTEGER  :: i, j
  REAL(KIND = kind_R)   :: vmax, rmax, renv, t, dt_dummy, xcen, ycen
  REAL(KIND = kind_R)   :: dx, dy, dist, wspd
  REAL(KIND = kind_R), PARAMETER :: d2r = 3.1415927/180.
  REAL(KIND = kind_R), PARAMETER :: rearth = 6371.e3

  CALL WW_ANNOUNCE('DUMMY_WND entered',3)

  vmax = 45.0
  rmax = 37000.
  renv = 232000.

  dt_dummy = dtc
  t = time_s + dt_dummy
  xcen = x0 + (nx/2)*sx -  0.75/(6*60*60)*t
  ycen = y0 + (ny/2)*sy

  DO j = 1, ny
  DO i = 1, nx
    dx = rearth * d2r * (x0 + (i-1)*sx - xcen) * cos(ycen*d2r)
    dy = rearth * d2r * (y0 + (j-1)*sy - ycen)
    dist = SQRT(dx*dx + dy*dy)
    IF ( dist < 1  .OR. dist > renv ) THEN
       sbc1(i,j,:) = 0.0
    ELSE
       IF ( dist <= rmax ) THEN
          wspd = vmax * dist / rmax
       ELSE
          wspd = vmax * rmax / dist
       END IF
       sbc1(i,j,1) = -dy * wspd / dist
       sbc1(i,j,2) = dx * wspd / dist
    ENDIF
  END DO
  END DO
 sbc = sbc1
END SUBROUTINE DUMMY_WND
!
!***********************************************************************
!
SUBROUTINE DUMMY_CUR
  USE WW_cc
  USE supl_vars_def
  IMPLICIT NONE
  INTEGER  :: i, j
  REAL(KIND = kind_R)   :: vmax, rmax, renv, t, dt_dummy, xcen, ycen
  REAL(KIND = kind_R)   :: dx, dy, dist, wspd
  REAL(KIND = kind_R), PARAMETER :: d2r = 3.1415927/180.
  REAL(KIND = kind_R), PARAMETER :: rearth = 6371.e3

  CALL WW_ANNOUNCE('DUMMY_CUR entered',3)

  vmax = 2.5
  rmax = 37000.
  renv = 232000.

  dt_dummy = dtc
  t = time_s + dt_dummy
  xcen = x0 + (nx/2)*sx -  0.75/(6*60*60)*t
  ycen = y0 + (ny/2)*sy

  DO j = 1, ny
  DO i = 1, nx
    dx = rearth * d2r * (x0 + (i-1)*sx - xcen) * cos(ycen*d2r)
    dy = rearth * d2r * (y0 + (j-1)*sy - ycen)
    dist = SQRT(dx*dx + dy*dy)
    IF ( dist < 1  .OR. dist > renv ) THEN
       sfcurx(i,j) = 0.0
       sfcury(i,j) = 0.0
    ELSE
       IF ( dist <= rmax ) THEN
          wspd = vmax * dist / rmax
       ELSE
          wspd = vmax * rmax / dist
       END IF
       sfcurx(i,j) = -dy * wspd / dist
       sfcury(i,j) = dx * wspd / dist
    ENDIF
  END DO
  END DO
END SUBROUTINE DUMMY_CUR
!
!***********************************************************************
!

