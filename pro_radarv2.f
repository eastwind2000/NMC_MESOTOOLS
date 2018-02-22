      program read_radar
      implicit none    
	include 'define_DataType.f'

	integer*4 i
	integer*4 year,mon,day
	integer*4 hour,minute,seconds
	integer*4 test
      
	type(radar_record) :: radar1
	type(radar_grads ) :: radar1ch


	real x0,y0,xc,yc
	integer nr,nElevationAngle
	character*40 fname, basedata_name

	integer flag, lev
	real    tim
	character*8 sid
      


      open(32,file='f:\gx07\radar\files.list')
	
c	do while(.not.eof(32))
  
      read(32,*)basedata_name
	open(12,file='f:\gx07\radar\200706\'//basedata_name,form='binary')

      print*, basedata_name
c      print*,'The recold length of SA Radar data: ', sizeof(radar1)
      
	nr=0
      nElevationAngle	= 1
	write(fname,'(a23,i2.2,a4)')'radar_'//trim(basedata_name) 
     $                             //'_ne',nElevationAngle,'.grd'
	open(22,file=fname,form='binary')

	print*, fname

	tim=0
      flag=1
	lev=1

	do while(.not.eof(12))
	   read(12)radar1
        
c	  if(nr==0)then
c	   write(21,*),nr,(dble(radar1%ElevationAngle)/8)*(180./4096.), 
c     $           radar1%ReflectivityGateSize, radar1%DopplerGateSize, 
c     $           radar1%ReflectivityGates,    radar1%DopplerGates
c        endif
	   

c	   write(21,*)  i,radar1%DataNumber,(dble(radar1%ElevationAngle)/8)*(180./4096.), 
c     $              (dble(radar1% AzimuthAngle)/8)*(180./4096.)
	
	   call radar_transform(radar1, radar1ch)

	   if (nr>1.and.radar1%DataNumber==1)then
	      
		  lev=0
	      write(22)sid, radar1ch%y_dbz(460), radar1ch%x_dbz(460),
     $               tim, lev, flag
	      close(22)

   	      nElevationAngle=nElevationAngle+1
	      write(fname,'(a23,i2.2,a4)')'radar_'//trim(basedata_name) 
     $                                  //'_ne',nElevationAngle,'.grd' 	      
	      open(22,file=fname,form='binary')

	      tim=0
	      flag=1
	      lev=1

	      do i=1,460
	         write(22)sid, radar1ch%y_dbz(i), radar1ch%x_dbz(i),
     $                 tim, lev, flag, radar1ch%dbz(i)
            enddo

         else

	      do i=1,460
	         write(22)sid, radar1ch%y_dbz(i), radar1ch%x_dbz(i),
     $                 tim, lev, flag, radar1ch%dbz(i)
c	         print*, radar1ch%y_dbz(i), radar1ch%x_dbz(i),
c     $                 radar1ch%dbz(i)
            enddo
            

	   endif
	   
	nr=nr+1
	enddo

      lev=0
	write(22)sid, radar1ch%y_dbz(460), radar1ch%x_dbz(460),
     $               tim, lev, flag
	close(22)


	print*, 'There are ',nr,' records in this Archive File'


	call get_date(radar1%radical_collect_date,year,mon,day)

	call GetRadialTime(radar1%radical_collect_time,hour,minute,seconds)
      

c	enddo

	end program 



	subroutine radar_transform(ra1, ra1ch)
	implicit none
	include 'define_DataType.f'
	type(radar_record) :: ra1
	type(radar_grads ) :: ra1ch
	integer*4 i,j,k
	real  alpha, De, DL

	real distance, dx,dy,dlat,dlon
	real refc,velc,swc
	real undef

	undef = -999.0

!	ra1ch%x0=117.717 ; 	ra1ch%y0=39.0439  ! Tianjing
      ra1ch%x0=109.4561 ; ra1ch%y0=24.3569 ! Liuzhou

      
!	ra1ch%x0=116.4719 ;  ra1ch%y0=39.8089  ! Beijing



	ra1ch % AzimuthAngle   = (dble(ra1%AzimuthAngle   )/8)*(180./4096.)
      ra1ch % ElevationAngle = (dble(ra1%ElevationAngle )/8)*(180./4096.)
      ra1ch % FirstGateRangeOfRef = dble(ra1%FirstGateRangeOfRef)
	ra1ch % FirstGateRangeOfDoppler=dble(ra1%FirstGateRangeOfDoppler)
	do k=1,11
	   if ( ra1ch%ElevationAngle >= (vcp21(k)-0.2).and. 
     $        ra1ch%ElevationAngle <= (vcp21(k)+0.2) )then
              ra1ch%ElevationAngle=vcp21(k)
	   endif
      enddo
      
	do i=1,460
	 

	   distance= ra1ch%FirstGateRangeOfRef + (i-1)*ra1%ReflectivityGateSize 

         DL=distance*cos(ra1ch%ElevationAngle*pi/180. )
	   
	   alpha=atan( DL/Re )

	   De=alpha*Re
         
	   dx=De*cos( (90-ra1ch%AzimuthAngle)*pi/180  )

	   dy=De*sin( (90-ra1ch%AzimuthAngle)*pi/180  )
         
	   dlon= (dx/(6371e3*cos( ra1ch%y0*pi/180 )))*180/pi


         dlat= (dy/6371e3)*180/pi
	       
 
	   ra1ch%x_dbz(i)=ra1ch%x0 + dlon 
         ra1ch%y_dbz(i)=ra1ch%y0 + dlat
         ra1ch%z_dbz(i)=distance*sin(ra1ch%ElevationAngle*pi/180. )
         
	   refc=real( ichar( ra1%dbz(i) )  )
	  
	   if(refc==0)then
	      ra1ch%dbz(i)=-999.0
	   else if (refc==1)then
	      ra1ch%dbz(i)=(refc-2.)/2.0-32.0
	   else
	      ra1ch%dbz(i)=(refc-2.)/2.0-32.0
	   endif
   	   
         
	enddo
      
	
	do i=1,920
	
         distance = ra1ch%FirstGateRangeOfDoppler +(i-1)*ra1%DopplerGateSize
         ra1ch%x_dpl(i)=distance*cos(ra1ch%ElevationAngle*pi/180. )
     $                          *cos( (90+ra1ch%AzimuthAngle)*pi/180  )
         ra1ch%y_dpl(i)=distance*cos(ra1ch%ElevationAngle*pi/180. )
     $                          *sin( (90+ra1ch%AzimuthAngle)*pi/180  )
         ra1ch%z_dpl(i)=distance*sin(ra1ch%ElevationAngle*pi/180. )	

	   velc=real( ichar( ra1%vel(i) )  )
	   swc =real( ichar( ra1%sw(i)  )  )
	   
	   if(velc==0)then
	      ra1ch%vel(i)=-999.0
	   else if (velc==1)then
	      if(ra1%VelResolution==2 ) ra1ch%vel(i)=(velc-2.)/2.0-63.5
            if(ra1%VelResolution==4 ) ra1ch%vel(i)=(velc-2.)/2.0-127.
   
	   else

	      if(ra1%VelResolution==2 ) ra1ch%vel(i)=(velc-2.)/2.0-63.5
            if(ra1%VelResolution==4 ) ra1ch%vel(i)=(velc-2.)/2.0-127.

	   endif
         

	   if(swc==0)then
	        ra1ch%sw(i)=-999.0
	   else if (swc==1)then
              ra1ch%sw(i)=(swc-2.)/2.0-63.5
  	   else
              ra1ch%sw(i)=(swc-2.)/2.0-63.5
	   endif
       
	
	enddo	



	end subroutine


      
	subroutine get_date(JulianDate,year, month,day)
	integer*2 JulianDate
	integer*4   year, month, day
	integer*4 JLDAYN,L,N,I,J
      JLDAYN = JulianDate+2440587
	  L  = JLDAYN + 68569 
	  N  = 4 * L / 146097 
        L  = L - (146097 * N + 3) / 4 
	  I  = 4000 * (L + 1) / 1461001 
	  L  = L - 1461 * I / 4 + 31 
	  J  = 80 * L / 2447 
	  day= L - 2447 * J / 80 
        L  = J / 11 
	  month  = J + 2 - 12 * L 
	  year  = 100 * (N - 49) + I + L 
	
         print*,year,month,day
	end subroutine




	subroutine GetRadialTime(liMilliSeconds, hour,minute,seconds )
	integer*4 liMilliSeconds
	integer*4 hour,minute,seconds

     	Seconds=liMilliSeconds/1000  
	hour = Seconds/3600 
	minute= (Seconds-hour*3600)/60 
	seconds = Seconds -(60*hour+minute)*60 

      print*,hour,minute,seconds

	end subroutine






	
	
	
	
	
	
	





