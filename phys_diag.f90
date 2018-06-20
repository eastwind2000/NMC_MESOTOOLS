
subroutine  phys_diag

        use diag
!============================================================================
       implicit none
       character*120 ::  workdir
       character*60  ::  case_name       
       integer       ::  nt    
       character*120 ::  gribprep_fname, diag_fname
!============================================================================
       
       integer,parameter:: nx=221, ny=121, nz=21
       
       real, dimension(nx,ny)   :: cape,cin,u10m, v10m,pwatm, psfc, slp, rh2m, tmp2m, hgtsfc
 	   real, dimension(nx,ny,nz)::h, u,v, rh,  tmp, w, vor, theta,tse,q
       
       real, dimension(nx,ny)::qs, pwatc

       real, dimension(nx,ny,nz)::pv1,pv2, mpv1,mpv2,                & 
                                hp, h1, h2, dens, omega, h1w, h2w, mf
       
	 real, dimension(nx,ny,nz):: prs, qvp, the, eth, sateth, twb, tdp
     real, dimension(nx,ny,nz):: uz, vz 
     
  !================adding vars by dynamic.gs Geostrophic Diag=======================
     
     real,dimension(nx,ny,nz):: vadv, tadv, qvec1,qvec2, FG_func, ug, vg
     real,dimension(nx,ny,nz):: sm_h, sm_tmp, sm_u, sm_v
     real,dimension(nx,ny)  :: lat, lon, fe
     
 !==================================================================================
	 
	 real lat0,lon0, dtheta, dz, dp, dx, dy

	 real plev(nz)
	
	 data plev/100000, 97500, 95000, 92500, 90000, 85000, 80000, 75000, 70000, 65000, 60000,  &
                  55000, 50000, 45000, 40000, 35000, 30000, 25000, 20000, 15000, 10000/

	 integer  i, j, k, it , icerel
	 real    Re, omega0,  y, f0 , dtr
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!================================== QG_DIAG==================================
     
     
     real,dimension(nx,ny,nz) ::  qg_qvec1, qg_qvec2, qg_FG_func, qg_tadv, qg_vadv, qg_ug, qg_vg
     integer                  ::  smth_method
     
     
!============================================================================
      smth_method = -1 
     
      f0 = 1e-4

      Re     =  6370.949e3
	  omega0 =  2*pi/(24*3600)
      
      lat0   =   10
      lon0   =   50 
      dtheta =   0.5
       
      dtr = pi/180.0
      
      do j=1, ny
      do i=1, nx
         lat(i,j)=lat0+(j-1)*dtheta
         lon(i,j)=lon0+(i-1)*dtheta     
      enddo
      enddo
      
      do j=1, ny
      do i=1, nx
         
         fe(i,j)= 2 * omega0 * sin( lat(i,j)*dtr )
          
	  enddo
      enddo
      
!	plev   =  plev*100                  ! pay attention
      
!==================================================================================   

       
!    gribprep_fname=trim(workdir)//trim(gribprep_fname)
       
!    diag_fname=trim(workdir)//trim(diag_fname)
       
	 open(12, file=trim(gribprep_fname),   form='binary', status='old')
	 open(22, file=trim(diag_fname),       form='binary')
       
       print*
       print*, trim(diag_fname)

       
	 do it=1,nt
     
	    print*, trim(case_name), ' :  ',  'phys_diag ---- :', it
        
        
        read(12)vadv, tadv, qvec1, qvec2, FG_func,   ug, vg
        read(12)cape,  cin,  u10m,  v10m,   pwatm, psfc,  slp, rh2m, tmp2m, hgtsfc
	    read(12)h,u,v,rh, tmp, omega, vor, theta, tse, q


!cccccccccccc  From RIPv4 source code cccccccccccccccccccccccccc
      
	do   k=1, nz
	do   j=1, ny
	do   i=1, nx
           
	     prs(i,j,k) = plev(k)            

	enddo
	enddo
    enddo
      

   
      icerel=1
	
      call  qvpcalc(rh2m, tmp2m, psfc/100, icerel, qs, nx, ny, 1)
       
	  call  qvpcalc(rh, tmp, prs/100, icerel, qvp, nx, ny, nz)       ! prs is in hPa!!  important!

      call  thecalc(prs/100, tmp, qvp, the, nx, ny, nz)
      
      call  eqthecalc(qvp, tmp, prs/100, eth, nx, ny, nz)
  
      call  sateqthecalc(tmp, prs/100, sateth, nx, ny, nz)

!	call  wetbulbcalc(prs, tmp, qvp, twb, nx, ny, nz)
      
	  call  tdpcalc(qvp, prs/100, tdp, nx, ny, nz)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
	  do i=1, nx
      do j=1, ny
         call pwatcalc(psfc(i,j), qs(i,j), qvp(i,j,1:nz), plev(1:nz), pwatc(i,j), nz)
      enddo
      enddo
    
        
        do k=1,nz
	    do j=1,ny
		do i=1,nx
             dens(i,j,k)=prs(i,j,k)/( rgas* tmp(i,j,k) )
	    enddo
	    enddo
	    enddo
  
		do k=1,nz
	    do j=1,ny
		do i=1,nx
             w(i,j,k)= -omega(i,j,k)/( grav * dens(i,j,k) )
	    enddo
	    enddo
	    enddo
    
 
        do k=2,nz-1
        do i=1,nx
	    do j=1,ny
           
           dp =  prs(i,j,k+1)-prs(i,j,k-1)
        
	       pv1(i,j,k) = -grav*( vor(i,j,k) + fe(i,j) ) *                         &
                              ( theta(i,j,k+1) - theta(i,j,k-1) )/dp
           
	      mpv1(i,j,k) = -grav*( vor(i,j,k) + fe(i,j) )   *						&
                          ( tse(i,j,k+1)- tse(i,j,k-1) ) /dp

        enddo
	    enddo
        enddo
        

        
	  do k=2,nz-1
	  do j=2,ny-1
	  do i=2,nx-1
		
	    
        dx =  Re*cos(lat(i,j)*dtr)*dtheta*dtr
	    dy =  Re*dtheta*dtr    
        dp =  prs(i,j,k+1)-prs(i,j,k-1)

        !c	mpv2(i,j,k)= ( (w(i,j+1,k)-w(i,j-1,k))/dy
!c    $                 -(v(i,j,k+1)-v(i,j,k))/dz )*(tse(i+1,j,k)-tse(i-1,j,k))/dx
!c    $               +( (u(i,j,k+1)-u(i,j,k))/dz
!c    $                 -(w(i+1,j,k)-w(i-1,j,k))/dx )*(tse(i,j+1,k)-tse(i,j-1,k))/dy		 	

	    pv2(i,j,k)=  grav * ( (v(i,j,k+1)-v(i,j,k-1))/dp ) *                       &
					        ( theta(i+1,j,k)  - theta(i-1,j,k))/(2*dx)             &							     &
                    -grav * ( (u(i,j,k+1)-u(i,j,k-1))/dp ) *                       &
                            (theta(i,j+1,k) - theta(i,j-1,k))/(2*dy)		
  

	   mpv2(i,j,k)=  grav * ( (v(i,j,k+1)-v(i,j,k-1))/dp ) *                &
                           (  tse(i+1,j,k)-tse(i-1,j,k) )/(2*dx)           &
				    -grav * ( (u(i,j,k+1)-u(i,j,k-1))/dp ) *                &
                           ( tse(i,j+1,k)-tse(i,j-1,k)  )/(2*dy)		
	    
	  
	  enddo
	  enddo
	  enddo

        
	do   k=2, nz-1
	do   j=2, ny-1
	do   i=2, nx-1
         
	    
        dx =  Re*cos(lat(i,j)*dtr)*dtheta*dtr
	    dy =  Re*dtheta*dtr   
        
        dp =  prs(i,j,k+1)-prs(i,j,k-1)
           
		 h1w(i,j,k)= u(i,j,k)*(  (w(i,j+1,k)-w(i,j-1,k))/(2*dy) +                       &
                                 dens(i,j,k)*grav*( v(i,j,k+1)-v(i,j,k-1))/dp   )

		 h2w(i,j,k)=-v(i,j,k)*(  (w(i+1,j,k)-w(i-1,j,k))/(2*dx) +                       &
                                dens(i,j,k)*grav*( u(i,j,k+1)-u(i,j,k-1))/dp   )
           
		 h1(i,j,k) =  u(i,j,k)*( dens(i,j,k)*grav*( v(i,j,k+1)-v(i,j,k-1))/dp   )
 
		 h2(i,j,k) = -v(i,j,k)*( dens(i,j,k)*grav*( u(i,j,k+1)-u(i,j,k-1))/dp   )
	
	     hp(i,j,k)= w(i,j,k)*vor(i,j,k)	
         
	enddo
	enddo
	enddo 

	do   k=1, nz
	do   j=1, ny
	do   i=1, nx
    
	         y= re* ( lat(i,j)-lat0)*pi/180
	        mf(i,j,k)= u(i,j,k) - f0 * y 
	enddo
	enddo
	enddo 
  
      do k=2, nz-1
      do i=2, nx-1
      do j=2, ny-1
                
        dx =  Re*cos(lat(i,j)*dtr)*dtheta*dtr
	    dy =  Re*dtheta*dtr   
	      
          
          vz(i,j,k)= - ( grav/( fe(i,j)*theta(i,j,k) ) )*( theta(i+1,j,k)-theta(i-1,j,k)  )/(2*dx)
 
          uz(i,j,k)= - ( grav/( fe(i,j)*theta(i,j,k) ) )*( theta(i,j+1,k)-theta(i,j-1,k)  )/(2*dy)
 
      enddo
      enddo
      enddo

!===================================================================================================
      
      smth_method = 306
      
      print*, ' =============== Getting into QG_dynamics, with smth_method = ', smth_method
     
      do  k = 1, nz
          
         
          call QG_dynamics( h(:,:,k), tmp(:,:,k), u(:,:,k), v(:,:,k), prs(:,:,k), lat, lon,  nx, ny,      &             
                            qg_qvec1(:,:,k), qg_qvec2(:,:,k), qg_FG_func(:,:,k), qg_tadv(:,:,k),          &
                            qg_vadv(:,:,k),     qg_ug(:,:,k),      qg_vg(:,:,k), smth_method,            &
                            sm_h(:,:,k),    sm_tmp(:,:,k), sm_u(:,:,k), sm_v(:,:,k)   )    
      enddo
      
!=================================Output Result ====================================================
      
        write(22) vadv, tadv, Qvec1,Qvec2, FG_func, ug, vg
        write(22) qg_vadv, qg_tadv, qg_qvec1,qg_qvec2, qg_FG_func, qg_ug, qg_vg
        write(22) sm_h, sm_tmp, sm_u, sm_v      
	    write(22) cape, cin, u10m, v10m, qvp, the, eth, sateth, tdp, rh,dens,  h, u, v, w,               &
                  omega, vor,tse,theta,                                                                  &
                  pv1,pv2,mpv1,mpv2, tmp, h1, h2, hp, h1w, h2w, mf,q,                                    &
                  pwatc, pwatm, psfc,slp, hgtsfc, qs, uz, vz
        
!===================================================================================================	 
	 enddo


	 close(12)

     close(22)

    end subroutine



!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
 

subroutine QG_dynamics(hgtprs, tmpprs, ugrdprs, vgrdprs, prs, lat, lon, nx, ny,     &
                       qvec1,  qvec2,  FG_func, tadv, vadv, ug, vg, smth_method,    &
                       sm_hgt, sm_tmp, sm_u, sm_v )

implicit none

integer  :: nx, ny, nz

!=======================Input Vars ================================
real,dimension(nx,ny)   :: hgtprs, tmpprs, ugrdprs, vgrdprs
real,dimension(nx,ny)   :: prs                                ! in Pa
real,dimension(nx,ny)   :: lat, lon
integer                 :: smth_method
!====================Output========================================
real,dimension(nx,ny):: ug, vg, qvec1, qvec2, FG_func 
real,dimension(nx,ny):: divg, vort, vadv, tadv

!=================== Work Vars==============================

real, dimension(nx,ny) :: sm_hgt, sm_tmp, sm_u, sm_v, work

real                   :: pi, dtr, Re, g, R, omega
real				   :: dlat, dlon, dx, dy
real, dimension(nx,ny) :: f

real, dimension(nx,ny) :: dugdx, dugdy, dvgdx, dvgdy, dtdx, dtdy, dudx, dudy, dvdx, dvdy,       &
                          grad_tmp, f1, f2, fnx, fny, fn, fs, def1, def2
!================================================================
integer                :: i, j


!==================================================================
 


  pi    = 3.14159
  dtr   = pi/180
  Re     = 6.37122e6
  omega = 7.2921e-5
  g     = 9.8
  R     = 287
  
  dlon = lon(2,1)-lon(1,1)
  dlat = lat(1,2)-lat(1,1) 
  
  do i=1, nx
  do j=1, ny
  
    f(i,j) = 2 * omega*sin( lat(i,j) *dtr )
    
  enddo
  enddo
  
      sm_hgt = hgtprs
      sm_tmp = tmpprs
      sm_u   = ugrdprs
      sm_v   = vgrdprs
      work   = 0.0 
  
  if ( smth_method > 0 ) then
       

      
      call smooth( sm_hgt, work, smth_method, nx, nx, ny)
      
      call smooth( sm_tmp, work, smth_method, nx, nx, ny)
      
      call smooth( sm_u,   work, smth_method, nx, nx, ny)
      
      call smooth( sm_v,   work, smth_method, nx, nx, ny)
      
      !hgtprs = sm_hgt
      !
      !tmpprs = sm_tmp
      !
      !ugrdprs = sm_u
      !
      !vgrdprs = sm_v
      
  endif
  

  
  do i=2, nx-1
  do j=2, ny-1

     dx =  Re * cos( lat(i,j) * dtr ) * dlon * dtr     
     dy =  Re * dlat * dtr   

     ug(i,j) = -1 *( g/f(i,j) ) * ( sm_hgt(i,j+1) - sm_hgt(i,j-1) )/(2*dy)
     vg(i,j) =     ( g/f(i,j) ) * ( sm_hgt(i+1,j) - sm_hgt(i-1,j) )/(2*dx) 
     
  enddo
  enddo

  
!'dugx=cdiff(ug,x)'
!'dugy=cdiff(ug,y)'
!'dvgx=cdiff(vg,x)'
!'dvgy=cdiff(vg,y)'
!
!'dtdx=cdiff(tmp,x)/dx'
!'dtdy=cdiff(tmp,y)/dy'  
  
  do i=2, nx-1
  do j=2, ny-1
  	
     dx =  Re * cos( lat(i,j) * dtr ) * dlon * dtr     
     dy =  Re * dlat * dtr   
	
     dugdx(i,j) = ( ug(i+1,j) - ug(i-1,j) ) / ( 2.0*dx )
     
     dugdy(i,j) = ( ug(i,j+1) - ug(i,j-1) ) / ( 2.0*dy )
     
     dvgdx(i,j) = ( vg(i+1,j) - vg(i-1,j) ) / ( 2.0*dx ) 
     
     dvgdy(i,j) = ( vg(i,j+1) - vg(i,j-1) ) / ( 2.0*dy )
     
     dudx(i,j) = ( sm_u(i+1,j) - sm_u(i-1,j) ) / ( 2.0*dx )
     
     dudy(i,j) = ( sm_u(i,j+1) - sm_u(i,j-1) ) / ( 2.0*dy )
     
     dvdx(i,j) = ( sm_v(i+1,j) - sm_v(i-1,j) ) / ( 2.0*dx ) 
     
     dvdy(i,j) = ( sm_v(i,j+1) - sm_v(i,j-1) ) / ( 2.0*dy )  
     
     dtdx(i,j) = ( sm_tmp(i+1,j) - sm_tmp(i-1,j) ) /( 2.0*dx  )
     
     dtdy(i,j) = ( sm_tmp(i,j+1) - sm_tmp(i,j-1) ) /( 2.0*dy  )
  
  enddo
  enddo
  
  do i=2, nx-1
  do j=2, ny-1
     
     vort(i,j) =    dvdx(i,j) - dudy(i,j)
     
	 divg(i,j) = -( dudx(i,j) + dvdy(i,j)  )
  enddo
  enddo
  
  
!'define Q1=-1*(R/p)*(dugx/dx*dtdx + dvgx/dx*dtdy)'
!'define Q2=-1*(R/p)*(dugy/dy*dtdx + dvgy/dy*dtdy)'
  
  do i=2, nx-1
  do j=2, ny-1
  
	 dx =  Re * cos( lat(i,j) * dtr ) * dlon * dtr     
     dy =  Re * dlat * dtr   
     
     Qvec1(i,j) = -1* (R/prs(i,j) ) * ( dugdx(i,j)*dtdx(i,j) + dvgdx(i,j)*dtdy(i,j) )
	 
     Qvec2(i,j) = -1* (R/prs(i,j) ) * ( dugdy(i,j)*dtdx(i,j) + dvgdy(i,j)*dtdy(i,j) )
     
     tadv(i,j) = -ugrdprs(i,j)*dtdx(i,j) - vgrdprs(i,j)*dtdy(i,j) 
     
     vadv(i,j) = -sm_u(i,j)* ( vort(i+1,j) - vort(i-1,j) )/(2.0*dx) - &
                  sm_v(i,j)* ( vort(i,j+1) - vort(i,j-1) )/(2.0*dy)   
     
  enddo
  enddo
  
!'def1=cdiff(uwnd,x)/dx-cdiff(vwnd,y)/dy-vwnd*tan(dtr*lat)/a'
!'def2=cdiff(vwnd,x)/dx+cdiff(uwnd,y)/dy+uwnd*tan(dtr*lat)/a'
!
  do i=2, nx-1
  do j=2, ny-1
	    
    def1(i,j)=  dudx(i,j)  -  dvdy(i,j)   -   sm_v(i,j)*tan(dtr*lat(i,j))/Re
    def2(i,j)=  dvdx(i,j)  +  dudy(i,j)   +   sm_u(i,j)*tan(dtr*lat(i,j))/Re
  enddo
  enddo
  

!f1=-( dtdx*(divg+def1)+dtdy*(vort+def2) ) / 2'
!'f2=(dtdx*(vort-def2)-dtdy*(divg-def1))/2'
!
!
!'fn=(dtdx*f1+dtdy*f2)/mag(dtdx,dtdy)'
!'fs=(dtdx*f2-dtdy*f1)/mag(dtdx,dtdy)'
!
!'fnx=-1*((dtdx*f1)/mag(dtdx,dtdy))*10e9'
!'fny=-1*((dtdy*f2)/mag(dtdx,dtdy))*10e9'  
  
  
  do i=2, nx-1
  do j=2, ny-1
     
 
    f1(i,j) = -( dtdx(i,j) * ( divg(i,j) + def1(i,j)) +        &
				 dtdy(i,j) * ( vort(i,j) + def2(i,j)) ) / 2.0
                 
    f2(i,j) =  ( dtdx(i,j) * ( vort(i,j) - def2(i,j)) -        &
                 dtdy(i,j) * ( divg(i,j) - def1(i,j)) ) / 2.0
  enddo
  enddo
  
  do i=2, nx-1
  do j=2, ny-1
  
    grad_tmp(i,j)=sqrt( dtdx(i,j)*dtdx(i,j) + dtdy(i,j)*dtdy(i,j)  )
    
  enddo
  enddo
  
  
  do i=2, nx-1
  do j=2, ny-1
  
  fn(i,j)=( dtdx(i,j)*f1(i,j)+dtdy(i,j)*f2(i,j) )/grad_tmp(i,j)
  fs(i,j)=( dtdx(i,j)*f2(i,j)-dtdy(i,j)*f1(i,j) )/grad_tmp(i,j)
!
  fnx(i,j)= -1*( (dtdx(i,j)*f1(i,j))/grad_tmp(i,j) )*10e9
  fny(i,j)= -1*( (dtdy(i,j)*f2(i,j))/grad_tmp(i,j) )*10e9
  
  FG_func(i,j) = ( fn(i,j) + fs(i,j) ) * 10e9 
  
  enddo
  enddo
  
    !
    !if ( smth_method > 0 ) then
    !  
    !  sm_hgt = hgtprs
    !  sm_tmp = tmpprs
    !  sm_u   = ugrdprs
    !  sm_v   = vgrdprs
    !  
    !endif
    !
  
  
!
!say '........'
!
!'define F=(fn+fs)*10e9'
  
  
endsubroutine
    
    
    
!
!
!function dynamic(args)
!
!'q dims'
!xline=sublin(result,2)
!yline=sublin(result,3)
!zline=sublin(result,4)
!tline=sublin(result,5)
!
!lons=subwrd(xline,6)' 'subwrd(xline,8)
!lats=subwrd(yline,6)' 'subwrd(yline,8)
!if(subwrd(zline,7)='Z');levs=subwrd(zline,6);else;levs=subwrd(zline,6)' 'subwrd(zline,8);endif
!time=subwrd(tline,6)
!
!check=1
!a=1
!v=1
!while(check=1)
!  line=subwrd(args,a)
!  if(line='');break;endif
!
!  if(line='-help')
!    help()
!    return
!  endif
!
!  if(line='-var');hgt=subwrd(args,a+1);tmp=subwrd(args,a+2);uwnd=subwrd(args,a+3);vwnd=subwrd(args,a+4);v=0;endif
!
!  a=a+1
!endwhile
!
!  if(v=1)
!    say ''
!    say 'You have Chosen the Default variables for Height, Temperature, and Wind'
!    say 'To specify your variables, type "-var" then your variable names in the order of Height, Temperature, Uwind, Vwind'
!    say ''
!    'hgt=hgtprs'
!    'tmp=tmpprs'
!    'vwnd=vgrdprs'
!    'uwnd=ugrdprs'
!  endif
!
!  if(v!=1)
!    say ''
!    say 'You have Chosen to specify the variables for Height, Temperature, and Wind'
!    say 'You have specified the variables as:'
!    say 'Height: 'hgt
!    say 'Temp:   'tmp
!    say 'U-wind: 'uwnd
!    say 'V-wind: 'vwnd
!    'hgt='hgt
!    'tmp='tmp
!    'vwnd='vwnd
!    'uwnd='uwnd
!  endif
!say 'Please wait while I calculate your variables'
!
!'pi=3.14159'
!'dtr=pi/180'
!'a=6.37122e6'
!'omega=7.2921e-5'
!'g=9.8'
!'R=287'
!
!'define f=2*omega*sin(lat*dtr)'
!'define p=lev*100'
!
!say '...'
!
!'dy=cdiff(lat,y)*dtr*a'
!'dx=cdiff(lon,x)*dtr*a*cos(lat*dtr)'
!
!'dhgtx=cdiff(hgt,x)'
!'dhgty=cdiff(hgt,y)'
!
!'define ug=-1*(g/f)*(dhgty/dy)'
!'define vg=(g/f)*(dhgtx/dx)'
!
!'define ua=uwnd-ug'
!'define va=vwnd-vg'
!
!say '....'
!
!'dugx=cdiff(ug,x)'
!'dugy=cdiff(ug,y)'
!'dvgx=cdiff(vg,x)'
!'dvgy=cdiff(vg,y)'
!
!'dtdx=cdiff(tmp,x)/dx'
!'dtdy=cdiff(tmp,y)/dy'
!
!say '.....'
!
!'define Q1=-1*(R/p)*(dugx/dx*dtdx + dvgx/dx*dtdy)'
!'define Q2=-1*(R/p)*(dugy/dy*dtdx + dvgy/dy*dtdy)'
!'define divq=hdivg(Q1,Q2)'
!'define tadv=(-uwnd*dtdx-vwnd*dtdy)'
!
!say '......'
!
!'divg=hdivg(uwnd,vwnd)'
!'vort=hcurl(uwnd,vwnd)'
!
!'dvdx=cdiff(vort,x)/dx'
!'dvdy=cdiff(vort,y)/dy'
!
!'define vadv=(-uwnd*dvdx-vwnd*dvdy)'
!
!say '.......'
!
!'def1=cdiff(uwnd,x)/dx-cdiff(vwnd,y)/dy-vwnd*tan(dtr*lat)/a'
!'def2=cdiff(vwnd,x)/dx+cdiff(uwnd,y)/dy+uwnd*tan(dtr*lat)/a'
!
!'f1=-(dtdx*(divg+def1)+dtdy*(vort+def2))/2'
!'f2=(dtdx*(vort-def2)-dtdy*(divg-def1))/2'
!
!
!'fn=(dtdx*f1+dtdy*f2)/mag(dtdx,dtdy)'
!'fs=(dtdx*f2-dtdy*f1)/mag(dtdx,dtdy)'
!
!'fnx=-1*((dtdx*f1)/mag(dtdx,dtdy))*10e9'
!'fny=-1*((dtdy*f2)/mag(dtdx,dtdy))*10e9'
!
!say '........'
!
!'define F=(fn+fs)*10e9'
!
!say 'DONE!'
!say ''
!
!say 'The following variables have been defined for the dimensions:'
!say 'Longitude: 'lons
!say 'Latitude: 'lats
!say 'Pressure Levels: 'levs
!say 'Time: 'time
!say '--------------------------------------------------------------------'
!say 'Defined Variables: Variable              Name          Units        '                 
!say '                  -Geostrophic Wind :    ug,vg         [m/s]        '
!say '                  -Ageostrophic Wind:    ua,va         [m/s]        '
!say '                  -Q-Vectors        :    Q1,Q1         [pa/m2/s]    '
!say '                  -Temp Advection   :    tadv          [K/s]        '  
!say '                  -Vort Advection   :    vadv          [-]          '
!say '                  -Frontogenesis    :    F             [K/m/s]x10^9 '
!say '                  -Fn Vector        :    fnx,fny       [K/m/s]x10^9 '
!say '                  -Deformation      :    def1,def2     [m]          '
!say '--------------------------------------------------------------------'
!say 'To plot a variable simply type "d " then pick a Name from the list  '
!
!
!return
!
!
!
  
  
      
      
      subroutine generate_diag_ctl(case_name, nt, timestr)
      
      use data_dir
      implicit none
      integer       :: nt
      character*60  :: case_name
      character*120 :: timestr
      character*120 :: diag_ctl_fname
    
    
    
      diag_ctl_fname=trim(home_dir)//trim(case_name)//'\diag_phys.ctl'
    
      print*, trim(diag_ctl_fname)
    
      open(22,file=trim(diag_ctl_fname))
      write(22,'(a)')'dset '//trim(home_dir)//trim(case_name)//'.\diag_phys.grd'
      write(22,'(a)')'title diag_phys '//trim(case_name)
      write(22,'(a)')'undef -999.0'
      write(22,'(a)')'xdef 221  linear  50 0.5'
      write(22,'(a)')'ydef 121  linear  10 0.5'
      write(22,'(a)')'zdef 21  levels  1000 975 950 925 900 850 800 750 700 650 600 550 500 450 400 350 300 250 200 150 100'
      write(22,'(a)')timestr
      write(22,'(a)')'vars   39'
      write(22,'(a)')'cape      00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'cin       00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'u10m      00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'v10m      00 0 GEOPHYSICAL HEIGHT'      
      write(22,'(a)')'qvp       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'the       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'eqthe     21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'sateqthe  21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'tdp       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'rh        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'dens      21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'h         21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'u         21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'v         21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'w         21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'omega     21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'vor       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'tse       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'theta     21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'pv1       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'pv2       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'mpv1      21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'mpv2      21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'tmp       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'h1        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'h2        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'hp        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'h1w       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'h2w       21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'mf        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'q         21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'pwatc     00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'pwatm     00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'psfc      00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'slp       00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'hgtsfc    00 0 GEOPHYSICAL HEIGHT'      
      write(22,'(a)')'qs        00 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'uz        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'vz        21 0 GEOPHYSICAL HEIGHT'
      write(22,'(a)')'endvars'    
    
      close(22)
   
    
    endsubroutine
     

    
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
!================================================================================
    
  
    
      subroutine draw_phys_diag( case_name )
    
      use data_dir
      implicit none

    
      character*60  case_name
!==============================================================================    
    
      open(22,file='draw_phys_diag.bat')
      write(22,*)'@echo on'
      write(22,*)'cd '//trim(home_dir)//trim(case_name)//'\'
      write(22,*)'copy /y '//'D:\WarmRain_GYHY2014\data_process\*.gs '//' .\'    
      write(22,*)'opengrads -blc draw_500h500w.gs'
      write(22,*)'opengrads -blc draw_500h700w.gs'
      write(22,*)'opengrads -blc draw_500h850w.gs'
      write(22,*)'opengrads -blc draw_500h925w.gs' 
      write(22,*)'opengrads -blc draw_cape.gs'
      write(22,*)'opengrads -blc draw_h200.gs'
      write(22,*)'opengrads -blc draw_pwat.gs'
      write(22,*)'opengrads -blc draw_hw58_tse.gs'
      write(22,*)'opengrads -blc draw_h500+slp.gs'  

      close(12)
    
      call system( 'draw_phys_diag.bat' ) 


      close(22)
    
    
    endsubroutine
    
    
!=========================================================    
