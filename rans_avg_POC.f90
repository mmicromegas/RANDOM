! proof-of-concept for extended rans_avg

program ransX_avg

  INTEGER, PARAMETER :: qx = 10, qy = 10, qz = 10
  INTEGER, PARAMETER :: nnuc = 25
  INTEGER, PARAMETER :: nrans = 5+nnuc
  INTEGER, PARAMETER :: ix = 16  
  
  integer*4 i, j, k, n, ii, jj, kk
  integer*4 ifield, rans_nnuc
  !
  real*8 velx(qx,qy,qz), vely(qx,qy,qz),velz(qx,qy,qz)
  real*8 density(qx,qy,qz),eint(qx,qy,qz),temp(qx,qy,qz)
  real*8 energy(qx,qy,qz),ekin(qx,qy,qz), press(qx,qy,qz)

  real*8 ei(qx), dd(qx), ux(qx), uy(qx), uz(qx), xn(qx,nnuc), tt(qx), pp(qx)

  real*8 fsterad(qy,qz)
  
  real*8 fsteradjk
  
  real*8 xnuc(qx,qy,qz,nnuc)
  
  ! averages
  
  real*8 havg(4,nrans,qx)

  real*8 time, rans_tavg, rans_end,  dt
  
  character(len=4)  :: xidchar    
 
!  character(len=24), allocatable :: ransname(:)
!  allocate(ransname(nrans))
  
!  character(len=24), allocatable :: xvarname(:)
!  allocate (xvarname(ix)) 
  
!  integer*4, allocatable :: iransname(:)
!  allocate (iransname(nnuc))
 
  character(len=24) ransname(nrans)
  character(len=24) xvarname(ix)
  integer*4 iransname(nnuc)  
  
  
  ! ransX

  integer*4 idd,iux,iuy,iuz,ienuc  
  
  real*8 eps(qx,qy,qz)
  
  real*8 enuc(qx)
  
  real*8 pf(qx)
  
  ! Reynolds fluctuations 
  
  real*8 ddf_r(qx,qy,qz), eif_r(qx,qy,qz)
  real*8 uxf_r(qx,qy,qz),uyf_r(qx,qy,qz),uzf_r(qx,qy,qz)
  real*8 xnf_r(qx,qy,qz), ppf_r(qx,qy,qz)

  ! Favre corrections
  
  real*8 uxf_c(qx),eh_uyf_f(qx),eh_uzf_f(qx),eh_eif_f(qx)  
  real*8 xnf_c(qx,nnuc)
  
  ! Favre fluctuations 
  
  real*8 ddf_f(qx,qy,qz), eif_f(qx,qy,qz)
  real*8 uxf_f(qx,qy,qz),uyf_f(qx,qy,qz),uzf_f(qx,qy,qz)  
  real*8 xnf_f(qx,qy,qz)
  
  ! Fluxes
  real*8 fh_feix(4,qx), fh_fxnx(4,qx,nnuc), fh_fxnxx(4,qx,nnuc), fh_xndotx(4,qx,nnuc)

  ! Variances  
  
  real*8 fh_sxn(4,qx,nnuc), fh_fsxnx(4,qx,nnuc), fh_sxndot(4,qx,nnuc) 
  
  !
  real*8 gradxppf_r(qx,qy,qz)
  real*8 eh_xnfgradpf(4,qx,nnuc), eh_xnf_f(4,qx,nnuc)
  
  ! Xdot
  
  real*8 xdot(qx,qy,qz,nnuc), xnew(qx,qy,qz,nnuc) 
  real*8 xdotn(qx,nnuc), xd(qx,qy,qz)
  real*8 fh_fxndotx(4,qx,nnuc), fh_rxx(4,qx)

  
! dummy xx
  real*8 xx(nnuc)
  
  
  ! imode = 0 : initialize
  ! imode = 1 : update
  
  imode = 1
     
  rans_nnuc = nnuc	 
	 
  ! initialize fields for new segment of time averaging (between data dumps)
  
  if(imode.eq.0) then
     rans_tstart = time
     rans_tavg   = 0.d0
     do n=1,nrans
        do i=1,qx
           havg(1,n,i) = 0.d0
           havg(2,n,i) = 0.d0
           havg(3,n,i) = 0.d0
        enddo
     enddo
  else if(imode.eq.1) then
        do i=1,qx
		   fh_feix(1,i) = fh_feix(2,i)
           fh_feix(2,i) = 0.d0
           do n=1,rans_nnuc		
              fh_fxnx(1,i,n) = fh_fxnx(2,i,n)
			  fh_fxnx(2,i,n)= 0.
          enddo
		enddo
  endif   
   
   
! populate test arrays
  
!  seed = 1234
!  call srand(seed)

  do k=1,qz
     do j=1,qy  
        do i=1,qx
		  density(i,j,k) = RAND(0) 
		  velx(i,j,k) = RAND(0) 
		  vely(i,j,k) = RAND(0) 
		  velz(i,j,k) = RAND(0) 
		  energy(i,j,k) = RAND(0) 		  
		  eps(i,j,k) = RAND(0)
		  temp(i,j,k) = RAND(0)
		  press(i,j,k) = RAND(0)
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
    
! calculate ek and ei  
  
  do k=1,qz
     do j=1,qy
        do i=1,qx
           ekin(i,j,k) = (velx(i,j,k)*velx(i,j,k) + & 
                vely(i,j,k)*vely(i,j,k) + velz(i,j,k)*velz(i,j,k))*0.5d0
           eint(i,j,k) = energy(i,j,k) - ekin(i,j,k)
        enddo
     enddo
  enddo  

  do k=1,qz
     do j=1,qy  
       fsterad(k,j) = 1.
     enddo
  enddo	 

  do n=1,nnuc
   xx(n) = 1.
  enddo  
  
  ! dtx = 1 second   
  do k=1,qz
    do j=1,qy  
      do i=1,qx
        !  call burnzone(tt(i),dd(i),xx,1.)
        do n=1,rans_nnuc
		   xnew(i,j,k,n) = xx(n)
        enddo	
        xdot(i,j,k,n) = xnew(i,j,k,n) - xnuc(i,j,k,n)
      enddo	
    enddo
  enddo
  
  
! calculate horizonatal averages   

  dtx = 1.
  
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
           tt(i)   = temp(i,j,k)
           pp(i)   = press(i,j,k)		   
		enddo
    
        do i=1,qx
           do n=1,rans_nnuc
              xn(i,n) = xnuc(i,j,k,n)
           enddo
        enddo
		
        do i=1,qx
		   do n=1,rans_nnuc
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
           havg(2,ifield,i) = havg(2,ifield,i) + ux(i)*fsteradjk ! eh_dd
        enddo
		
        ! uy (3)
        ifield = ifield + 1
		iuy = ifield
        if(imode.eq.0) ransname(ifield) = 'uy'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + uy(i)*fsteradjk ! eh_dd
        enddo
		
        ! uz (4)
        ifield = ifield + 1
		iuz = ifield
        if(imode.eq.0) ransname(ifield) = 'uz'
        do i=1,qx
           havg(2,ifield,i) = havg(2,ifield,i) + uz(i)*fsteradjk ! eh_dd				
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
 
           ! ddxdot
           ifield = ifield+1                  
		   ransname(ifield) = xvarname(7)		   
           do i=1,qx                                      
              havg(2,ifield,i) =  havg(2,ifield,i) + &    
                   dd(i)*xdotn(i,n)*fsteradjk                      
           enddo		  
        enddo		   
     enddo
  enddo
  
  ! ransX

  do i=1,qx
     ! get Reynolds fluctuations
	 ddf_r(i,:,:) =  density(i,:,:) - havg(2,idd,i) 
     uxf_r(i,:,:) =  velx(i,:,:) - havg(2,iux,i) 		  
     uyf_r(i,:,:) =  vely(i,:,:) - havg(2,iuy,i) 
     uzf_r(i,:,:) =  velz(i,:,:) - havg(2,iuz,i)
     ppf_r(i,:,:) =  press(i,:,:) - havg(2,ipp,i) 	 
	 ! get Favre corrections	
	 uxf_c(i) = (sum(density(i,:,:)*uxf_r(i,:,:)*fsterad(:,:)))/havg(2,iux,i)
	 ! get Favre fluctuations
	 uxf_f(i,:,:) =  uxf_r(i,:,:) - uxf_c(i)
  enddo	

  ! get grad P'
  
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
     do n=1,rans_nnuc
	 !  if(imode.eq.0) then ! construct varnames
          9902 format(i0.4)
          write(xidchar,9902) n
          xvarname(8) = 'fh_x' // xidchar //'x'          
          xvarname(9) = 'fh_x' // xidchar//'xx'  !
          xvarname(10) = 'fh_fsx'//xidchar       ! 
          xvarname(11) = 'fh_fsx'//xidchar//'x' ! 
          xvarname(12) = 'fh_sxn'//xidchar//'dot'
          xvarname(13) = 'fh_fxn'//xidchar//'dotx'
	      xvarname(14) = 'eh_x' // xidchar //'f_f'  
		  xvarname(15) = 'eh_x' // xidchar //'fgradpf'
		  xvarname(16) = 'fh_rxx' // xidchar		  
     !  endif
	 
	 
        ! get Reynolds composition fluctuations	 
		xnf_r(i,:,:) =  xnuc(i,:,:,n) - havg(2,iransname(n),i)
		! get Favre corrections
		xnf_c(i,n)   = (sum(density(i,:,:)*xnf_r(i,:,:)*fsterad(:,:)))/havg(2,iransname(n),i) 
        ! get Favre composition fluctuations		
		xnf_f(i,:,:) = xnf_r(i,:,:) - xnf_c(i,n)
		! get nuclear composition rate of change
		xd(i,:,:) = xdot(i,:,:,n)

	    ! get composition flux
        ifield = ifield+1		
		ransname(ifield) = xvarname(8)	
		fh_fxnx(2,i,n)  = sum(xnf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get flux of composition flux
		ifield = ifield+1		
		ransname(ifield) = xvarname(9)
		fh_fxnxx(2,i,n) = sum(xnf_f(i,:,:)*uxf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))	

        ! get composition variance
		ifield = ifield+1		
		ransname(ifield) = xvarname(10)
        fh_sxn(2,i,n) = sum(xnf_f(i,:,:)*xnf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get flux of variance
		ifield = ifield+1		
		ransname(ifield) = xvarname(11)
		fh_fsxnx(2,i,n) = sum(xnf_f(i,:,:)*xnf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))

		! get nuclear-variance coupling
		ifield = ifield+1		
		ransname(ifield) = xvarname(12)
		fh_sxndot(2,i,n) = sum(xnf_f(i,:,:)*xd(i,:,:)*density(i,:,:)*fsterad(:,:))

        ! get nuclear flux 
		ifield = ifield+1		
		ransname(ifield) = xvarname(13)
        fh_fxndotx(2,i,n) = sum(uxf_f(i,:,:)*xd(i,:,:)*density(i,:,:)*fsterad(:,:))

        ! get Reynolds average of X''
		ifield = ifield+1		
		ransname(ifield) = xvarname(14)		
        eh_xnf_f(2,i,n) = sum(xnf_f(i,:,:)*fsterad(:,:)) 	
		
        ! get pressure-force coupling
		ifield = ifield+1		
		ransname(ifield) = xvarname(15)			
        eh_xnfgradpf(2,i,n) = sum(xnf_f(i,:,:)*gradxppf_r(i,:,:)*fsterad(:,:)) 		
	 enddo
     ! get Reynolds stress (rr)
	 ifield = ifield+1		
	 ransname(ifield) = xvarname(16)		  
     fh_rxx(2,i) = sum(uxf_f(i,:,:)*uxf_f(i,:,:)*density(i,:,:)*fsterad(:,:))
  enddo 
  


  write(*,*) fh_sxn(2,1,:)
 
  ! update RANS running average: fh_fxxx(3,:)
  
  time = 1.
  rans_tavg = 0.
  rans_end = time
  dt = 0.1
  
  if(imode.eq.1) then 
     rans_tavg = rans_tavg + dt
!     rans_tend = time
     do i=1,qx 
!         fh_feix(3,i) = fh_feix(3,i) + &
!             (fh_feix(2,i) + fh_feix(1,n))*0.5d0*dt
		do n=1,rans_nnuc
!		  write(*,*) i,n
          fh_fxnx(3,i,n) = fh_fxnx(3,i,n) + &
             (fh_fxnx(2,i,n) + fh_fxnx(1,i,n))*0.5d0*dt
        enddo		
     enddo	 
  else if(imode.eq.0) then
     do i=1,qx
!        fh_feix(4,i) = fh_feix(2,i) !store first instance in this averaging interval
		do n=1,rans_nnuc
		  fh_fxnx(4,i,n) = fh_fxnx(2,i,n)
		enddo
     enddo
  endif 
 
   ! store 
 
!  return

end 
 
!end subroutine ransX_avg
  
  
  
  
  
  