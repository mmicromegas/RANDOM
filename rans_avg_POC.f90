! proof-of-concept for extended rans_avg

program ransX_avg

  implicit none

  INTEGER, PARAMETER :: qx = 10, qy = 10, qz = 10
  INTEGER, PARAMETER :: nnuc = 25
  
  ! the 2* is due to xdotn, the 8 is due to terms in comp equations (transport, flux, variance)
  ! the +1 is due to rxx, the +1 is due to fhhx, the +1 je ddcp,+1 je ttrms, +1 ddrms, +1 sddx, +1 rdfdddivux, +1 fddx, +1 fsddx
  INTEGER, PARAMETER :: nrans = 9+2*nnuc+9*nnuc+1+1+1+1+1+1+1+1+1
  INTEGER, PARAMETER :: ix = 16
  
  integer*4 i, j, k, n, ii, jj, kk
  integer*4 ifield, rans_nnuc
  !
  real*8 velx(qx,qy,qz), vely(qx,qy,qz),velz(qx,qy,qz)
  real*8 density(qx,qy,qz),eint(qx,qy,qz),temp(qx,qy,qz),cpdat(qx,qy,qz)
  real*8 energy(qx,qy,qz),ekin(qx,qy,qz),press(qx,qy,qz),enth(qx,qy,qz) 
  real*8 divel(qx,qy,qz)

  real*8 ei(qx),dd(qx),ux(qx),uy(qx),uz(qx),xn(qx,nnuc),tt(qx),pp(qx),hh(qx),divu(qx)

  real*8 fsterad(qy,qz)
  
  real*8 fsteradjk
  
  real*8 xnuc(qx,qy,qz,nnuc)
  
  ! averages
  
  real*8 havg(4,nrans,qx)

  real*8 time,rans_tavg,rans_tend,rans_tstart,rans_end, dt
  
  character(len=4)  :: xidchar    
 
!  character(len=24), allocatable :: ransname(:)
!  allocate(ransname(nrans))
  
!  character(len=24), allocatable :: xvarname(:)
!  allocate (xvarname(ix)) 
  
!  integer*4, allocatable :: iransname(:)
!  allocate (iransname(nnuc))
 
  
! dummy xx
  real*8 xx(nnuc)
  real*8 dx  

  real*8 eps(qx,qy,qz)
  
  real*8 enuc(qx)
  
!  real*8 eh_uxf_f(qx)
!  real*8 eh_xnf_f(qx,nnuc)  
  
  ! ransX
  
  real*8 geomxM(qx,qy,qz),geomyM(qx,qy,qz),geomzM(qx,qy,qz)
  real*8 gxm(qx),gym(qx),qzm(qx)

  character(len=24) ransname(nrans)
  character(len=24) xvarname(ix)
  
  ! Index store for composition variables
  
  integer*4 iransname(nnuc)   
  
  ! Indices for specific thermodyanmic variables 
  
  integer*4 idd,iux,iuy,iuz,ienuc,iei,ipp,ihh,itt  
  
  ! Reynolds fluctuations 
  
  real*8 ddf_r(qx,qy,qz),eif_r(qx,qy,qz),enf_r(qx,qy,qz)
  real*8 uxf_r(qx,qy,qz),uyf_r(qx,qy,qz),uzf_r(qx,qy,qz)
  real*8 xnf_r(qx,qy,qz),ppf_r(qx,qy,qz),hhf_r(qx,qy,qz)

  ! Favre corrections
  
  real*8 uxf_c(qx),eh_uyf_f(qx),eh_uzf_f(qx),eh_eif_f(qx)
  real*8 hhf_c(qx),ttf_c(qx),ddf_c(qx)  
  real*8 xnf_c(qx,nnuc), uyf_c(qx),uzf_c(qx)
  
  ! Favre fluctuations 
  
  real*8 ddf_f(qx,qy,qz),eif_f(qx,qy,qz),hhf_f(qx,qy,qz)
  real*8 uxf_f(qx,qy,qz),uyf_f(qx,qy,qz),uzf_f(qx,qy,qz)  
  real*8 xnf_f(qx,qy,qz)
  
  ! Pressure fluctuations for grad P' 
  
  real*8 pf(qx)  
  
  ! Fluxes
  real*8 fh_feix(4,qx), fxnx(4,qx,nnuc), fxnxx(4,qx,nnuc)

  ! Variances  
  
  real*8 sxn(4,qx,nnuc), fsxnx(4,qx,nnuc), rfxndot(4,qx,nnuc) 
  
  ! Pressure-variance coupling terms
  
  real*8 gradxppf_r(qx,qy,qz)
  real*8 rxnfgradpf(4,qx,nnuc) 
  
  ! Reynolds average of X''
  
  real*8 rxnf_f(4,qx,nnuc)
  
  ! Xdot terms
  
  real*8 xdot(qx,qy,qz,nnuc), xnew(qx,qy,qz,nnuc) 
  real*8 xdotn(qx,nnuc), xd(qx,qy,qz)
  real*8 fxndotx(4,qx,nnuc)
  
  ! geometry terms
  real*8 galpha(4,qx,nnuc)
  
  ! Density flux equation
  
  real*8 sddx(4,qx),rdfdddivux(4,qx),fddx(4,qx),fsddx(4,qx)
  real*8 divelx(qx,qy,qz)  
 
  ! Other usefull fields
  
  real*8 rxx(4,qx), fhhx(4,qx),ddcp(4,qx),ttrms(4,qx),ddrms(4,qx)  
  
  integer*4 imode
  
  integer*4 ifieldmark 
  
  ! imode = 0 : initialize
  ! imode = 1 : update
  
  imode = 1
     
  rans_nnuc = nnuc	 
	 
  ! ransX initialize fields for new segment of time averaging (between data dumps)
  
  if(imode.eq.0) then
     rans_tstart = time
     rans_tavg   = 0.d0
     do n=1,qx
        rxx(1,i)   = 0.d0
        rxx(2,i)   = 0.d0
        rxx(3,i)   = 0.d0		
		fhhx(1,i)  = 0.d0  
		fhhx(2,i)  = 0.d0
		fhhx(3,i)  = 0.d0	
		ddcp(1,i)  = 0.d0 
		ddcp(2,i)  = 0.d0 
		ddcp(3,i)  = 0.d0 	
        ttrms(1,i) = 0.d0 		
        ttrms(2,i) = 0.d0 
        ttrms(3,i) = 0.d0
        ddrms(1,i) = 0.d0 		
        ddrms(2,i) = 0.d0 
        ddrms(3,i) = 0.d0
		sddx(1,i)  = 0.d0
 		sddx(2,i)  = 0.d0
		sddx(3,i)  = 0.d0
		fddx(1,i)  = 0.d0
 		fddx(2,i)  = 0.d0
		fddx(3,i)  = 0.d0		
		rdfdddivux(1,i) = 0.d0
		rdfdddivux(2,i) = 0.d0
		rdfdddivux(3,i) = 0.d0
		fsddx(1,i) = 0.d0
		fsddx(2,i) = 0.d0
		fsddx(3,i) = 0.d0
        do i=1,nrans
           havg(1,n,i) = 0.d0
           havg(2,n,i) = 0.d0
           havg(3,n,i) = 0.d0		   
        enddo
     enddo
  else if(imode.eq.1) then
        do i=1,qx
		   rxx(1,i)   = rxx(2,i)
		   rxx(2,i)   = 0.d0
		   fhhx(1,i)  = fhhx(2,i)
		   fhhx(2,i)  = 0.d0
		   ddcp(1,i)  = ddcp(2,i)
		   ddcp(2,i)  = 0.d0
		   ttrms(1,i) = ttrms(2,i)
		   ttrms(2,i) = 0.d0
		   ddrms(1,i) = ddrms(2,i)
		   ddrms(2,i) = 0.d0	
		   sddx(1,i)  = sddx(2,i)
		   sddx(2,i)  = 0.d0			   
		   fddx(1,i)  = fddx(2,i)
		   fddx(2,i)  = 0.d0	   
		   rdfdddivux(1,i) = rdfdddivux(2,i)
		   rdfdddivux(2,i) = 0.d0		
		   fsddx(1,i) = fsddx(2,i)
		   fsddx(2,i) = 0.d0			   
           do n=1,rans_nnuc		
              fxnx(1,i,n)    = fxnx(2,i,n)
			  fxnx(2,i,n)    = 0.
              fxnxx(1,i,n)   = fxnxx(2,i,n)
			  fxnxx(2,i,n)   = 0.			  
			  sxn(1,i,n)     = sxn(2,i,n)
			  sxn(2,i,n)     = 0.
			  fsxnx(1,i,n)   = fsxnx(2,i,n)
			  fsxnx(2,i,n)   = 0.
			  rfxndot(1,i,n) = rfxndot(2,i,n)
			  rfxndot(2,i,n) = 0.
			  fxndotx(1,i,n) = fxndotx(2,i,n)
			  fxndotx(2,i,n) = 0.
			  rxnf_f(1,i,n)  = rxnf_f(2,i,n)
		      rxnf_f(2,i,n)  = 0.
              rxnfgradpf(1,i,n) = rxnfgradpf(2,i,n)
              rxnfgradpf(2,i,n) = 0.	  
          enddo
		enddo
  endif   
   
! populate test arrays
  
!  seed = 1234
!  call srand(seed)

  do k=1,qz
     do j=1,qy  
        do i=1,qx
!		  density(i,j,k) = i*2. 
          density(i,j,k) = RAND(0)
		  velx(i,j,k) = RAND(0) 
		  vely(i,j,k) = i*4. 
		  velz(i,j,k) = i*5. 
		  energy(i,j,k) = RAND(0)
		  temp(i,j,k) = i*7.
		  press(i,j,k) = i*8.
          cpdat(i,j,k) = RAND(0)
          divelx(i,j,k) = RAND(0)
          divel(i,j,k) = RAND(0)
          geomxM(i,j,k) = RAND(0)
          geomyM(i,j,k) = RAND(0)
          geomzM(i,j,k) = RAND(0)		  
		enddo
	 enddo
  enddo


  do n=1,rans_nnuc
	do k=1,qz
      do j=1,qy  
        do i=1,qx
           xnuc(i,j,k,n) = RAND(0)
        enddo
      enddo
    enddo
  enddo	
    
! calculate ek and ei and enthalpy (specific)
  
  do k=1,qz
     do j=1,qy
        do i=1,qx
           ekin(i,j,k) = (velx(i,j,k)*velx(i,j,k) + & 
                vely(i,j,k)*vely(i,j,k) + velz(i,j,k)*velz(i,j,k))*0.5d0
           eint(i,j,k) = energy(i,j,k) - ekin(i,j,k)
		   enth(i,j,k) = eint(i,j,k) + press(i,j,k)/density(i,j,k)
        enddo
     enddo
  enddo  

  do k=1,qz
     do j=1,qy  
       fsterad(k,j) = 1.
     enddo
  enddo	 
  
  do n=1,nnuc
   xx(n) = RAND(0)
  enddo  
  
  ! dtx = 1 second   
  do k=1,qz
    do j=1,qy  

        do i=1,qx
           dd(i)   = density(i,j,k)
           tt(i)   = temp(i,j,k)		   
        enddo	
	
        do i=1,qx
           do n=1,rans_nnuc
              xn(i,n) = xnuc(i,j,k,n)
           enddo
        enddo
		
        do i=1,qx
		   xx(:) = xn(i,:)
        !  call burnzone(tt(i),dd(i),xx,1.)
		   xdot(i,j,k,:) = xx - xnuc(i,j,k,:) 
		enddo   
		
    enddo
  enddo  
  
! calculate horizonatal averages   
  
  do k=1,qz
     do j=1,qy  
        ! mesh descriptors
        fsteradjk = fsterad(j,k)	 

        do i=1,qx
           dd(i)   = density(i,j,k)
           ux(i)   = velx(i,j,k)
           uy(i)   = vely(i,j,k)
           uz(i)   = velz(i,j,k)
           ei(i)   = eint(i,j,k)
           pp(i)   = press(i,j,k)		
           tt(i)   = temp(i,j,k)
           hh(i)   = enth(i,j,k)	
		enddo
    
        do i=1,qx
           do n=1,rans_nnuc
              xn(i,n) = xnuc(i,j,k,n)
			  xdotn(i,n) = xdot(i,j,k,n)			  
           enddo
        enddo
	
        ! update current horizontally averaged fields in havg(2,:,:)
        
        ! dd  (1)
        ifield = 1
		idd = ifield
        if(imode.eq.0) ransname(ifield) = 'dd'
        do i=1,qx 
           havg(2,ifield,i) = havg(2,ifield,i)  + dd(i)*fsteradjk ! eh_dd
        enddo  
		
		! ux (2)
        ifield = ifield + 1
		iux = ifield
        if(imode.eq.0) ransname(ifield) = 'ux'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + ux(i)*fsteradjk ! eh_ux
        enddo
		
        ! uy (3)
        ifield = ifield + 1
		iuy = ifield
        if(imode.eq.0) ransname(ifield) = 'uy'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + uy(i)*fsteradjk ! eh_uy
        enddo
		
        ! uz (4)
        ifield = ifield + 1
		iuz = ifield
        if(imode.eq.0) ransname(ifield) = 'uz'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + uz(i)*fsteradjk ! eh_uz				
        enddo

        ! ei (5)
        ifield = ifield +1
		iei = ifield
        if(imode.eq.0) ransname(ifield) = 'ei'
        do i=1,qx
           havg(2,ifield,i) =  havg(2,ifield,i) + &
                ei(i)*fsteradjk
        enddo		
	
        ! pp (6)
        ifield = ifield + 1
		ipp = ifield
        if(imode.eq.0) ransname(ifield) = 'pp'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + pp(i)*fsteradjk ! eh_pp
        enddo		
		
        ! hh (7)
        ifield = ifield + 1
		ihh = ifield
        if(imode.eq.0) ransname(ifield) = 'hh'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + hh(i)*fsteradjk ! eh_hh
        enddo		
		
        ! tt (8)
        ifield = ifield + 1
		itt = ifield
        if(imode.eq.0) ransname(ifield) = 'tt'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + tt(i)*fsteradjk ! eh_tt
        enddo		
	
        ! dddivu (9)
        ifield = ifield+1
        if(imode.eq.0) ransname(ifield) = 'dddivu'
        do i=1,qx
           havg(2,ifield,i) =  havg(2,ifield,i) + &
                dd(i)*divu(i)*fsteradjk
        enddo
	
	
       ! innermass ??
        ! ddgg ( ?? )
!        ifield = ifield + 1
!        if(imode.eq.0) ransname(ifield) = 'ddgg'
!        do i=1,qx
!            havg(2,ifield,i) = havg(2,ifield,i)  &
!                -dd(i)*bigG*innermass(i)/xznl(i)**2.d0*fsteradjk
!        enddo
!        endif	
	
        do n=1,rans_nnuc   
!           if(imode.eq.0) then ! construct varnames
              9901 format(i0.4)
              write(xidchar,9901) n
              xvarname(1) = 'x' // xidchar          
              xvarname(2) = 'x' // xidchar//'sq' ! reynolds X
              xvarname(3) = 'ddx'//xidchar       ! reynolds X density
              xvarname(4) = 'ddx'//xidchar//'sq' ! for favrian sigma
              xvarname(5) = 'x'//xidchar//'ux'
              xvarname(6) = 'ddx'//xidchar//'ux'
			  xvarname(7) = 'x' // xidchar //'dot' 			  
!           endif

           ! xn 
           ifield = ifield+1                  
! the imode has to change back to O
!           if(imode.eq.0) then
		     ransname(ifield) = xvarname(1)
			 iransname(n) = ifield
!		   endif 
           do i=1,qx                                      
              havg(2,ifield,i) =  havg(2,ifield,i) + &    
                   xn(i,n)*fsteradjk                      
           enddo	
		   
           ifield = ifield+1                  
		   ransname(ifield) = xvarname(7)		   
           do i=1,qx                                      
              havg(2,ifield,i) =  havg(2,ifield,i) + &    
                   dd(i)*xdotn(i,n)*fsteradjk                      
           enddo			   
		 		   
        enddo		   
     enddo
  enddo
    	
  ifieldmark = ifield		
		
!  write(*,*) havg(2,iransname(25),:)
		
  ! get Reynolds fluctuations 

   ! ransX

  do i=1,qx
  
     ! get Reynolds fluctuations
	 
     ddf_r(i,:,:) =  density(i,:,:) - havg(2,idd,i) 
     uxf_r(i,:,:) =  velx(i,:,:)    - havg(2,iux,i) 		  
     uyf_r(i,:,:) =  vely(i,:,:)    - havg(2,iuy,i) 
     uzf_r(i,:,:) =  velz(i,:,:)    - havg(2,iuz,i)
     ppf_r(i,:,:) =  press(i,:,:)   - havg(2,ipp,i) 	
     hhf_r(i,:,:) =  enth(i,:,:)    - havg(2,ihh,i)	 
	 
	 ! get Favre corrections	

	 uxf_c(i) = (sum(density(i,:,:)*uxf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i)
	 uyf_c(i) = (sum(density(i,:,:)*uyf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i)
	 uzf_c(i) = (sum(density(i,:,:)*uzf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i)	 
	 hhf_c(i) = (sum(density(i,:,:)*hhf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i)	 
	 ddf_c(i) = (sum(density(i,:,:)*ddf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i)	 

	 ! get Favre fluctuations

	 uxf_f(i,:,:) =  uxf_r(i,:,:) - uxf_c(i)
	 uyf_f(i,:,:) =  uyf_r(i,:,:) - uyf_c(i)
	 uzf_f(i,:,:) =  uzf_r(i,:,:) - uzf_c(i)	 
	 hhf_f(i,:,:) =  hhf_r(i,:,:) - hhf_c(i)
	 ddf_f(i,:,:) =  ddf_r(i,:,:) - ddf_c(i)	 

  enddo	
  
  ! get Grad P'
  
  dx = 1. 
  do k=1,qz
     do j=1,qy    

        do i=1,qx
	     pf(i) = ppf_r(i,j,k)
	    enddo
	
        do i=2,qx-1 
         gradxppf_r(i,j,k) = 0.5d0*(pf(i+1) - pf(i-1))/dx
        enddo  
		
	 enddo
  enddo  


  do i=1,qx
     ifield = ifieldmark  ! 56
     do n=1,rans_nnuc

	 !  if(imode.eq.0) then ! construct varnames
          9902 format(i0.4)
          write(xidchar,9902) n
          xvarname(8)  = 'fx' // xidchar //'x'          
          xvarname(9)  = 'fx' // xidchar//'xx'  !
          xvarname(10) = 'sx'//xidchar       ! 
          xvarname(11) = 'fsx'//xidchar//'x' ! 
          xvarname(12) = 'rfx'//xidchar//'dot'
          xvarname(13) = 'fx'//xidchar//'dotx'
	      xvarname(14) = 'rx' // xidchar //'f_f'  
		  xvarname(15) = 'rx' // xidchar //'fgradpf'	
		  xvarname(15) = 'galpha' // xidchar 		  
     !  endif

        ! get Reynolds composition fluctuations	 
		xnf_r(i,:,:) =  xnuc(i,:,:,n) - havg(2,iransname(n),i)
		! get Favre corrections
		xnf_c(i,n)   = (sum(density(i,:,:)*xnf_r(i,:,:)*fsterad(:,:)))/havg(2,idd,i) 
        ! get Favre composition fluctuations		
		xnf_f(i,:,:) = xnf_r(i,:,:) - xnf_c(i,n)	 
		! get nuclear composition rate of change
		xd(i,:,:) = xdot(i,:,:,n)


	    ! get composition flux
        ifield = ifield+1	
        ! fxnx (57, 81, etc.) 		
		ransname(ifield) = xvarname(8)	
		fxnx(2,i,n)  = sum(xnf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get flux of composition flux
		ifield = ifield+1		
		ransname(ifield) = xvarname(9)
		fxnxx(2,i,n) = sum(xnf_f(i,:,:)*uxf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))	

        ! get composition variance
		ifield = ifield+1		
		ransname(ifield) = xvarname(10)
        sxn(2,i,n) = sum(xnf_f(i,:,:)*xnf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get flux of variance
		ifield = ifield+1		
		ransname(ifield) = xvarname(11)
		fsxnx(2,i,n) = sum(xnf_f(i,:,:)*xnf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get nuclear-variance coupling
		ifield = ifield+1		
		ransname(ifield) = xvarname(12)
		rfxndot(2,i,n) = sum(xnf_f(i,:,:)*xd(i,:,:)*density(i,:,:)*fsterad(:,:))

        ! get nuclear flux 
		ifield = ifield+1		
		ransname(ifield) = xvarname(13)
        fxndotx(2,i,n) = sum(uxf_f(i,:,:)*xd(i,:,:)*density(i,:,:)*fsterad(:,:))

        ! get Reynolds average of X''
		ifield = ifield+1		
		ransname(ifield) = xvarname(14)		
        rxnf_f(2,i,n) = sum(xnf_f(i,:,:)*fsterad(:,:)) 	
		
        ! get pressure-force coupling
		ifield = ifield+1		
		ransname(ifield) = xvarname(15)			
        rxnfgradpf(2,i,n) = sum(xnf_f(i,:,:)*gradxppf_r(i,:,:)*fsterad(:,:)) 	

	   ! get geometry terms for the flux equation
		ifield = ifield+1		
		ransname(ifield) = xvarname(16)			
        galpha(2,i,n) = sum((density(i,:,:)*xnf_f(i,:,:)*uyf_f(i,:,:)*uyf_f(i,:,:) + &
		                     density(i,:,:)*xnf_f(i,:,:)*uzf_f(i,:,:)*uzf_f(i,:,:) + &
							 density(i,:,:)*xnf_f(i,:,:)*vely(i,:,:)*vely(i,:,:) + &
							 density(i,:,:)*xnf_f(i,:,:)*velz(i,:,:)*velz(i,:,:))*fsterad(:,:)) 	
	  
	 enddo
     ! get Reynolds stress (rxx)
	 ifield = ifield+1		
     ransname(ifield) = 'rxx'	  
     rxx(2,i) = sum(uxf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

     ! get enthalpy flux (fhhx)
	 ifield = ifield + 1
     ransname(ifield) = 'fhhx'
     fhhx(2,i) = sum(uxf_f(i,:,:)*hhf_f(i,:,:)*density(i,:,:)*fsterad(:,:)) 	 

     ! heat capacity at constant pressure 
	 ifield = ifield + 1
     ransname(ifield) = 'ddcp'
     ddcp(2,i) = sum(cpdat(i,:,:)*density(i,:,:)*fsterad(:,:))  
 
     ! get Trms
	 ifield = ifield + 1
     ransname(ifield) = 'ttrms'
     ttrms(2,i) = (sum(((temp(i,:,:)-havg(2,itt,i))**2)*fsterad(:,:)))**0.5

     ! get Drms
	 ifield = ifield + 1
     ransname(ifield) = 'ddrms'
     ddrms(2,i) = (sum(((density(i,:,:)-havg(2,idd,i))**2)*fsterad(:,:)))**0.5	 

     ! get Favrian density variance (sddx)
	 ifield = ifield + 1
     ransname(ifield) = 'sddx'
     sddx(2,i) = sum(ddf_f(i,:,:)*ddf_f(i,:,:)*density(i,:,:)*fsterad(:,:)) 	 
	
    ! get density fluctuation velocity coupling
    ! divelx has to be added and constructed from divux !!	
	 ifield = ifield + 1
     ransname(ifield) = 'rdfdddivux'	
	 rdfdddivux(2,i) = sum(uxf_f(i,:,:)*ddf_f(i,:,:)*density(i,:,:)*density(i,:,:)*divelx(i,:,:)*fsterad(:,:)) 	 

    ! get density flux 
	 ifield = ifield + 1
     ransname(ifield) = 'fddx'	
	 fddx(2,i) = sum(uxf_f(i,:,:)*ddf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

    ! get variance density flux 
	 ifield = ifield + 1
     ransname(ifield) = 'fsddx'	
	 fsddx(2,i) = sum(uxf_f(i,:,:)*uxf_f(i,:,:)*ddf_f(i,:,:)*density(i,:,:)*fsterad(:,:))	 
	 	 
	enddo    
 
  write(*,*) ifield
 
!  write(*,*) xnf_c 

!   write(*,*) havg(2,idd,:)  
 
!  write(*,*) xnf_r(:,:,1)
!  write(*,*) eh_xnf_f(:,:)
  
!  write(*,*) uxf_f(:,:,1) 
!  write(*,*) ddf_r(:,:,:)	
!  write(*,*) fh_fxnx(2,:,:)
!  write(*,*) fh_sxn(2,1,:) 
   write(*,*) ttrms(2,:)
!   write(*,*) ransname(1)
 
 
  ! update RANS running average: fh_fxxx(3,:)
  
  time = 1.
  rans_tavg = 0.
  rans_end = time
  dt = 0.1
  
  if(imode.eq.1) then 
     rans_tavg = rans_tavg + dt
     rans_tend = time
     do i=1,qx 
         rxx(3,i)   = rxx(3,i)   + (rxx(2,i)   + rxx(1,n))*0.5d0*dt
         fhhx(3,i)  = fhhx(3,i)  + (fhhx(2,i)  + fhhx(1,n))*0.5d0*dt			 
         ddcp(3,i)  = ddcp(3,i)  + (ddcp(2,i)  + ddcp(1,n))*0.5d0*dt				 
         ttrms(3,i) = ttrms(3,i) + (ttrms(2,i) + ttrms(1,n))*0.5d0*dt	
         ddrms(3,i) = ddrms(3,i) + (ddrms(2,i) + ddrms(1,n))*0.5d0*dt
         sddx(3,i)  = sddx(3,i)  + (sddx(2,i)  + sddx(1,n))*0.5d0*dt				 
         fddx(3,i)  = fddx(3,i)  + (fddx(2,i)  + fddx(1,n))*0.5d0*dt
         fsddx(3,i) = fsddx(3,i) + (fsddx(2,i) + fsddx(1,n))*0.5d0*dt			 
		do n=1,rans_nnuc
          fxnx(3,i,n) = fxnx(3,i,n) + &
             (fxnx(2,i,n) + fxnx(1,i,n))*0.5d0*dt
		  fxnxx(3,i,n) = fxnxx(3,i,n) + &
             (fxnxx(2,i,n) + fxnxx(1,i,n))*0.5d0*dt
          sxn(3,i,n) = sxn(3,i,n) + &
             (sxn(2,i,n) + sxn(1,i,n))*0.5d0*dt
          fsxnx(3,i,n) = fsxnx(3,i,n) + &
             (fsxnx(2,i,n) + fsxnx(1,i,n))*0.5d0*dt
          rfxndot(3,i,n) = rfxndot(3,i,n) + &  
             (rfxndot(2,i,n) + rfxndot(1,i,n))*0.5d0*dt
          fxndotx(3,i,n) = fxndotx(3,i,n) + &
             (fxndotx(2,i,n) + fxndotx(1,i,n))*0.5d0*dt  
          rxnf_f(3,i,n) = rxnf_f(3,i,n) + &
              (rxnf_f(2,i,n) + rxnf_f(1,i,n))*0.5d0*dt
          rxnfgradpf(3,i,n) = rxnfgradpf(3,i,n) + & 		  
              (rxnfgradpf(2,i,n) + rxnfgradpf(1,i,n))*0.5d0*dt
        enddo		
     enddo	 
  else if(imode.eq.0) then
     do i=1,qx
        rxx(4,i)   = rxx(2,i) !store first instance in this averaging interval
        fhhx(4,i)  = fhhx(2,i)
        ddcp(4,i)  = ddcp(2,i)
        ttrms(4,i) = ttrms(2,i)
        ddrms(4,i) = ddrms(2,i)
        sddx(4,i)  = sddx(2,i)		
		fddx(4,i)  = fddx(2,i)
		fsddx(4,i) = fsddx(2,i)		
        rdfdddivux(4,i) = rdfdddivux(2,i) 		
		do n=1,rans_nnuc
		  fxnx(4,i,n)  = fxnx(2,i,n)
		  fxnxx(4,i,n) = fxnxx(2,i,n)
		  sxn(4,i,n)   = sxn(2,i,n)
		  fsxnx(4,i,n) = fsxnx(2,i,n)
		  rfxndot(4,i,n) = rfxndot(2,i,n)
		  fxndotx(4,i,n) = fxndotx(2,i,n)
		  rxnf_f(4,i,n) = rxnf_f(2,i,n)
		  rxnfgradpf(4,i,n) = rxnfgradpf(2,i,n)
		enddo
     enddo
  endif 
 
   ! store 
 
!  return

end 
 
!end subroutine ransX_avg
  
  
  
  
  
  