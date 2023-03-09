c************************************************************************
c       Program for the computation of the b-value, the b-positive, 
c     	b-more positive and b-more incomplete.
c       Format input file (2 columns) occurrence time,magnitude
c	Output file bmth.dat contains the estmate of b(m_{th}) as function of time
c       Output file b+dmth.dat contains the estmate of b(m+dmth) (b-positive) as function of time	
c       Output file b++dmth.dat contains the estmate of b(m++dmth) (b-more-positive) as function of time
c	Output file b+fdmth.dat contains the estmate of b(m+fdmth) (b-more-incomplete) as function of time	
	
c*************************************************************************
c**************** VARIABLE DECLARATION ***********************************
c*************************************************************************
	
	real*8 itime(0:10000000),q(0:10000000)
	real*8 itime2(0:10000000),q2(0:10000000)
	character*40 file_name

c*************************************************************************





	

c*************************************************************************
c**************** OUTPUT FILE ********************************************
c*************************************************************************

	open(80,file='bmth.dat',status='unknown')
	open(81,file='b+dmth.dat',status='unknown')	
	open(82,file='b++dmth.dat',status='unknown')
	open(86,file='b+fdmth.dat',status='unknown')

c*************************************************************************	








	
c*************************************************************************
c******	User Input Parameters ********************************************
c*************************************************************************
	
	print*,'Please type the catalogue file name'
	read(*,*)file_name
	open(70,file=file_name,status='unknown')
	print*,'Please type the initial time'
	read(*,*)tmin                                 
	print*,'Please type the final time in day'
	read(*,*)tmax
	tmax=tmax*3600*24                             

c************************************************************************





	
c*************************************************************************
c******	Internal Parameters **********************************************
c*************************************************************************
	
	kstep=400		!# of events in each temporal bin
	kstep2=230		!# of events in each temporal bin in the filtered catalog with PhiA 
	qmin=3			!#minimum magnitude mth for b(mth)
	dqmin=0.2		!#minimum magnitude difference dmth for b+(dmth)
	tau=100			!blind time tau in b+f(tau) !TAKE CARE OF TIME UNITS

c*************************************************************************
	
	



c*************************************************************************
c******	Reading catalogue ************************************************
c*************************************************************************
	
	jj=0
	nz=0
	do i=1,100000000
	   read(70,*,end=99)tz,qz                !reading occurrence time and magnitude  
	   if(tz.ge.tmin.and.tz.le.tmax)then     !filtering catalogue in the chosen time window
	      jj=jj+1
	      itime(jj)=tz
	      q(jj)=qz+0.1*(rand()-0.5)      	!roundings to 0.1 precision according to Godano et al.(2014) 
	   endif
	enddo
 99	continue
	close(70)
	nevent=jj



	
c************************************************************************
c********************Filtering with PhiA ********************************
c************************************************************************

	nevent2=1
	q2(nevent2)=q(1)
	itime2(nevent2)=itime(1)
	do jj=2,nevent
	   do ii=jj-1,1,-1
	      dt=itime(jj)-itime(ii)
	      if(dt.gt.tau)exit          	 !applying blind time 
      	      if(q(jj).lt.q(ii))goto 56          !considering only \delta m greater than 0
	   enddo
	   nevent2=nevent2+1
	   q2(nevent2)=q(jj)
	   itime2(nevent2)=itime(jj)
 56	   continue
	enddo

c************************************************************************



c************************************************************************
c************Computing bmth, b+dmth and b++dmth *************************
c************************************************************************
	
	do ii=3,nevent-kstep,20
	   
	   nt=0
	   ntp=0
	   ntpp=0
	   ntap=0
	   aveq=0.
	   aveqp=0.
	   aveqpp=0.
	   aveqap=0.
		
	   do k=1,kstep
	      jj=ii+kstep-k

	      if(q(jj).ge.qmin)then
		 aveq=aveq+q(jj)            !computing the average of q for bmth
		 nt=nt+1
	      endif
	   
	      dq=q(jj)-q(jj-1)              !\delta m for b-positive 

	      if(dq.ge.dqmin)then     
		 ntp=ntp+1
		 aveqp=aveqp+dq             !computing the average of q for b+dmth
		 if(q(jj-1).ge.q(jj-2))then
		    ntpp=ntpp+1
		    aveqpp=aveqpp+dq        !computing the average of q for b++dmth
		 endif
	      endif
	      
	   enddo

	   bbb=1d0/(aveq/nt-qmin)/log(10.)            !b
	   bbbp=1d0/(aveqp/ntp-dqmin)/log(10.)        !b-postive
	   bbbpp=1d0/(aveqpp/ntpp-dqmin)/log(10.)     !b-more

	  
	   if(nt.ge.10)write(80,*)itime(ii),bbb*(nt-1.)/nt,
     &          bbb/sqrt(1.*nt),nt	      
	   if(ntp.ge.10)write(81,*)itime(ii),bbbp*(ntp-1.)/ntp,
     &   	bbbp/sqrt(1.*ntp),ntp
	   if(ntpp.ge.10)write(82,*)itime(ii),bbbpp*(ntpp-1.)/ntpp,
     &	        ntpp

	   open(86,file='bftau.dat',status='unknown')


c********* Computing b-more-incomplete **********************************
	   
	   if(ii.le.nevent2-kstep2)then
	      do k=1,kstep2
		 jj=ii+kstep2-k
		 dq=q2(jj)-q2(jj-1)

		 if(dq.ge.0)then
		    ntap=ntap+1
		    aveqap=aveqap+dq
		 endif
	      enddo
	      bbbap=1d0/(aveqap/ntap-0)/log(10.)
	      
	      if(ntap.ge.10)write(86,*)itime2(ii),bbbap*(ntap-1.)/ntap,bbbap/sqrt(1.*ntap),ntap
	   endif
	   call flush(80)
	   call flush(81)
	   call flush(83)
	   call flush(86)
	   
	enddo	
	stop
	end
	
c******************************************************************************************



