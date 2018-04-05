! snippet for ransX_avg new approach


  ! Composition variable names
  INTEGER, PARAMETER :: ix = 18
  character(len=24) xvarname(ix)

  ! TEMP STORAGE
  real*8 enth(qx,qy,qz)
  real*8 xx(rans_nnuc)

  ! Index store for composition variables  
  ! integer*4 iransnuc(rans_nnuc)   
  
  ! Indices for specific thermodyanmic variables 
  
  integer*4 idd,iux,iuy,iuz,ienuc,ipp,ihh,itt,iss  
  
  ! Reynolds fluctuations 
  
  real*8 ddf_r(qx,qy,qz),eif_r(qx,qy,qz),enf_r(qx,qy,qz)
  real*8 uxf_r(qx,qy,qz),uyf_r(qx,qy,qz),uzf_r(qx,qy,qz)
  real*8 ttf_r(qx,qy,qz),ppf_r(qx,qy,qz),hhf_r(qx,qy,qz)
  real*8 ssf_r(qx,qy,qz)
  real*8 xnf_r(qx,qy,qz,rans_nnuc)

  real*8 ddfr(qx),hhfr(qx),ttfr(qx)
  real*8 uxfr(qx),uyfr(qx),uzfr(qx)
  real*8 ssfr(qx),ppfr(qx),gxppfr(qx)
  real*8 xnfr(qx,rans_nnuc)
  
  ! Favre corrections
  
  real*8 uxf_c(qx), uyf_c(qx),uzf_c(qx)
  real*8 hhf_c(qx),ttf_c(qx),ddf_c(qx)  
  real*8 xnf_c(qx,rans_nnuc)
  
  ! Favre fluctuations 
  
  real*8 ddf_f(qx,qy,qz),eif_f(qx,qy,qz),hhf_f(qx,qy,qz)
  real*8 uxf_f(qx,qy,qz),uyf_f(qx,qy,qz),uzf_f(qx,qy,qz)  
  real*8 xnf_f(qx,qy,qz,rans_nnuc)

  real*8 uxff(qx),uyff(qx),uzff(qx)
  real*8 hhff(qx),ddff(qx)
  real*8 xnff(qx,rans_nnuc)
  
  ! Pressure fluctuations for grad P' 
  real*8 pf(qx)

  ! For pressure-variance coupling terms
  real*8 gxppf_r(qx,qy,qz)

  ! Xdot terms  
  real*8 xdot(qx,qy,qz,rans_nnuc) 
  real*8 xdn(qx,rans_nnuc)
  
  ! For density flux equation
  real*8 divelx(qx,qy,qz)   
 
  ! Other required fields 
  integer*4 ifieldmark   

  real*8 cpd(qx),dvlx(qx)
   
  integer*4 nsubcycle
  real*8    eps,snu,deltat,dtnew

  ! for MPI_ALLREDUCE
  integer*4 root  
  real*8 havg_sum(4,nrans,qx )
  real*8 uxf_c_sum(qx),uyf_c_sum(qx),uzf_c_sum(qx)
  real*8 hhf_c_sum(qx),ddf_c_sum(qx)
  real*8 xnf_c_sum(qx,rans_nnuc)
  
  
  if(imode.eq.0) then
     do k=1,qz
        do j=1,qy  
           do i=1,qx
              ddf_r(i,j,k) = 0.d0
              ttf_r(i,j,k) = 0.d0
              uxf_r(i,j,k) = 0.d0
              uyf_r(i,j,k) = 0.d0
              uzf_r(i,j,k) = 0.d0
              hhf_r(i,j,k) = 0.d0
              ppf_r(i,j,k) = 0.d0
              ssf_r(i,j,k) = 0.d0
              uxf_f(i,j,k) = 0.d0
              uyf_f(i,j,k) = 0.d0
              uzf_f(i,j,k) = 0.d0
              ddf_f(i,j,k) = 0.d0
              hhf_f(i,j,k) = 0.d0
              do n=1,rans_nnuc
                 xdot(i,j,k,n) = 0.d0
                 xnf_r(i,j,k,n) = 0.d0
                 xnf_f(i,j,k,n) = 0.d0
              enddo
           enddo
        enddo
     enddo
  endif
  
  
    ! ransX get nuclear data dX/dt  
  deltat = 1. ! 1 second

  do k=1,qz
    do j=1,qy  

        do i=1,qx
           dd(i)   = densty(i,j,k)
           tt(i)   = temp(i,j,k)		   
        enddo	
	
        do i=1,qx
           do n=1,rans_nnuc
              xn(i,n) = xnuc(i,j,k,n)
           enddo
        enddo
		
        do i=1,qx
           xx(:)    = xn(i,:)
           call burnzone(tt(i),dd(i),xx,eps,snu,nsubcycle,deltat,dtnew)
           xdot(i,j,k,:) = (xx(:) - xn(i,:))/deltat 
        enddo
		
     enddo
  enddo
  
  
   !call MPI_BARRIER(commcart, ierr)

  !call MPI_ALLREDUCE(havg(1,1,1), havg_sum, 4*qx*nrans, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  ! ransX 
  ifieldmark = ifield	
  
  do k=1,qz
     do j=1,qy
        ! mesh descriptors
        fsteradjk = fsterad(j,k)
       
        do i=1,qx
           dd(i)   = densty(i,j,k)
           ux(i)   = velx(i,j,k)
           uy(i)   = vely(i,j,k)
           uz(i)   = velz(i,j,k)
           pp(i)   = press(i,j,k)
           tt(i)   = temp(i,j,k)
           ek(i)   = ekin(i,j,k)
           et(i)   = energy(i,j,k)
           ei(i)   = et(i) - ek(i)
           ss(i)   = entropy(i,j,k)
           hh(i)   = ei(i) + pp(i)/dd(i)
           do n=1,rans_nnuc
              xn(i,n) = xnuc(i,j,k,n)
           enddo
        enddo

        ! get Reynolds fluctuations
        do i=1,qx
           ddf_r(i,j,k) =  dd(i) - havg(2,idd,i)
           ttf_r(i,j,k) =  tt(i) - havg(2,itt,i)          
           uxf_r(i,j,k) =  ux(i) - havg(2,iux,i) 		  
           uyf_r(i,j,k) =  uy(i) - havg(2,iuy,i) 
           uzf_r(i,j,k) =  uz(i) - havg(2,iuz,i)
           ppf_r(i,j,k) =  pp(i) - havg(2,ipp,i) 	
           hhf_r(i,j,k) =  hh(i) - havg(2,ihh,i)
           ssf_r(i,j,k) =  ss(i) - havg(2,iss,i)
        enddo

        do n=1,rans_nnuc
           do i=1,qx
              xnf_r(i,j,k,n) =  xn(i,n) - havg(2,iransnuc(n),i)
           enddo
        enddo
        
        ! get Reynolds fluctuations
        !do i=1,qx
        !   ddf_r(i,j,k) =  dd(i) - havg_sum(2,idd,i)
        !   ttf_r(i,j,k) =  tt(i) - havg_sum(2,itt,i)           
        !   uxf_r(i,j,k) =  ux(i) - havg_sum(2,iux,i) 		  
        !   uyf_r(i,j,k) =  uy(i) - havg_sum(2,iuy,i) 
        !   uzf_r(i,j,k) =  uz(i) - havg_sum(2,iuz,i)
        !   ppf_r(i,j,k) =  pp(i) - havg_sum(2,ipp,i) 	
        !   hhf_r(i,j,k) =  hh(i) - havg_sum(2,ihh,i)
        !   ssf_r(i,j,k) =  ss(i) - havg_sum(2,iss,i)
        !enddo

        !do n=1,rans_nnuc
        !   do i=1,qx
        !      xnf_r(i,j,k,n) =  xn(i,n) - havg_sum(2,iransnuc(n),i)
        !   enddo
        !enddo
        
     enddo
  enddo
  
  do i=1,qx
     uxf_c(i) = 0.d0
     uyf_c(i) = 0.d0
     uzf_c(i) = 0.d0
     hhf_c(i) = 0.d0
     ddf_c(i) = 0.d0
     do n=1,rans_nnuc
        xnf_c(i,n) = 0.d0
     enddo
  enddo
  
  do k=1,qz
     do j=1,qy
        ! mesh descriptors
        fsteradjk = fsterad(j,k)
        
        do i=1,qx
           dd(i)   = densty(i,j,k)
           ddfr(i) = ddf_r(i,j,k)
           uxfr(i) = uxf_r(i,j,k)
           uyfr(i) = uyf_r(i,j,k)
           uzfr(i) = uzf_r(i,j,k)          
           hhfr(i) = hhf_r(i,j,k)
           do n=1,rans_nnuc
              xnfr(i,n) = xnf_r(i,j,k,n)
           enddo
        enddo

        ! get Favre corrections
        do i=1,qx	
           uxf_c(i) = uxf_c(i) + (dd(i)*uxfr(i)*fsteradjk)/havg(2,idd,i)
           uyf_c(i) = uyf_c(i) + (dd(i)*uyfr(i)*fsteradjk)/havg(2,idd,i)
           uzf_c(i) = uzf_c(i) + (dd(i)*uzfr(i)*fsteradjk)/havg(2,idd,i)
           hhf_c(i) = hhf_c(i) + (dd(i)*hhfr(i)*fsteradjk)/havg(2,idd,i)
           ddf_c(i) = ddf_c(i) + (dd(i)*ddfr(i)*fsteradjk)/havg(2,idd,i)
        enddo

        do n=1,rans_nnuc
           do i=1,qx
              xnf_c(i,n) = xnf_c(i,n) + &
                   (dd(i)*xnfr(i,n)*fsteradjk)/havg(2,idd,i)
           enddo
        enddo

        ! get Favre corrections
        !do i=1,qx	
        !   uxf_c(i) = uxf_c(i) + &
        !        (dd(i)*uxfr(i)*fsteradjk)/havg_sum(2,idd,i)
        !   uyf_c(i) = uyf_c(i) + &
        !        (dd(i)*uyfr(i)*fsteradjk)/havg_sum(2,idd,i)
        !   uzf_c(i) = uzf_c(i) + &
        !        (dd(i)*uzfr(i)*fsteradjk)/havg_sum(2,idd,i)
        !   hhf_c(i) = hhf_c(i) + &
        !        (dd(i)*hhfr(i)*fsteradjk)/havg_sum(2,idd,i)
        !   ddf_c(i) = ddf_c(i) + &
        !        (dd(i)*ddfr(i)*fsteradjk)/havg_sum(2,idd,i)
        !enddo

        !do n=1,rans_nnuc
        !   do i=1,qx
        !      xnf_c(i,n) = xnf_c(i,n) + &
        !           (dd(i)*xnfr(i,n)*fsteradjk)/havg_sum(2,idd,i)
        !   enddo
        !enddo
        
     enddo
  enddo

  !call MPI_BARRIER(commcart, ierr)

  !call MPI_ALLREDUCE(uxf_c(1), uxf_c_sum, qx, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  !call MPI_ALLREDUCE(uyf_c(1), uyf_c_sum, qx, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  !call MPI_ALLREDUCE(uzf_c(1), uzf_c_sum, qx, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  !call MPI_ALLREDUCE(hhf_c(1), hhf_c_sum, qx, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  !call MPI_ALLREDUCE(ddf_c(1), ddf_c_sum, qx, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)

  !call MPI_ALLREDUCE(xnf_c(1,1), xnf_c_sum, qx*rans_nnuc, &
  !     MPI_DOUBLE_PRECISION, MPI_SUM, commyz, ierr)
  
  do k=1,qz
     do j=1,qy
        
        do i=1,qx
           ddfr(i) = ddf_r(i,j,k)
           uxfr(i) = uxf_r(i,j,k)
           uyfr(i) = uyf_r(i,j,k)
           uzfr(i) = uzf_r(i,j,k)          
           hhfr(i) = hhf_r(i,j,k)
           ttfr(i) = ttf_r(i,j,k)
           do n=1,rans_nnuc
              xnfr(i,n) = xnf_r(i,j,k,n)
           enddo
        enddo

        ! get Favre fluctuations
        do i=1,qx
           uxf_f(i,j,k) =  uxfr(i) - uxf_c(i)
           uyf_f(i,j,k) =  uyfr(i) - uyf_c(i)
           uzf_f(i,j,k) =  uzfr(i) - uzf_c(i)	 
           hhf_f(i,j,k) =  hhfr(i) - hhf_c(i)
           ddf_f(i,j,k) =  ddfr(i) - ddf_c(i)
        enddo

        do n=1,rans_nnuc
           do i=1,qx
              xnf_f(i,j,k,n) = xnfr(i,n) - xnf_c(i,n)
           enddo
        enddo

        ! get Favre fluctuations
        !do i=1,qx
        !   uxf_f(i,j,k) =  uxfr(i) - uxf_c_sum(i)
        !   uyf_f(i,j,k) =  uyfr(i) - uyf_c_sum(i)
        !   uzf_f(i,j,k) =  uzfr(i) - uzf_c_sum(i)	 
        !   hhf_f(i,j,k) =  hhfr(i) - hhf_c_sum(i)
        !   ddf_f(i,j,k) =  ddfr(i) - ddf_c_sum(i)	       
        !enddo

        !do n=1,rans_nnuc
        !   do i=1,qx
        !      xnf_f(i,j,k,n) = xnfr(i,n) - xnf_c_sum(i,n)
        !   enddo
        !enddo
        
     enddo
  enddo
  
  ! get Grad P'  
  
  do k=1,qz
     do j=1,qy
        
        do i=1,qx
           pf(i) = ppf_r(i,j,k)
        enddo
	
        do i=2,qx-1 
           ii         = coords(1)*qx + i
           dx         = (gxznr(ii) - gxznl(ii))
           gxppf_r(i,j,k) = 0.5d0*(pf(i+1) - pf(i-1))/dx
        enddo
		
     enddo
  enddo    
     
  do k=1,qz
     do j=1,qy
        fsteradjk = fsterad(j,k)
        ifield = ifieldmark

        do i=1,qx
           dd(i) = densty(i,j,k)
           uy(i)   = vely(i,j,k)
           uz(i)   = velz(i,j,k)
           uxff(i) = uxf_f(i,j,k)
           uyff(i) = uyf_f(i,j,k)
           uzff(i) = uzf_f(i,j,k)
           hhff(i) = hhf_f(i,j,k)
           ddff(i) = ddf_f(i,j,k)
           ddfr(i) = ddf_r(i,j,k)
           ttfr(i) = ttf_r(i,j,k)
           ppfr(i) = ppf_r(i,j,k)
           ssfr(i) = ssf_r(i,j,k)
           gxppfr(i) = gxppf_r(i,j,k)
           dvlx(i) = divelx(i,j,k)
           cpd(i) = cpdat(i,j,k)
           do n=1,rans_nnuc
              xnfr(i,n) = xnf_r(i,j,k,n)
              xnff(i,n) = xnf_f(i,j,k,n)
              xdn(i,n)  = xdot(i,j,k,n)
           enddo
        enddo
        
  
