cc Silhouette score
        parameter(maxd=11000,maxfea=4000,maxc=90)
	real*8  ddd(maxd,maxfea),pairdis(maxd,maxd)
     &   ,xa(maxd),xb(maxd),sl(maxd)
     & ,cent(maxc,maxfea),avedis(maxc),avedis2(maxc)
     &   ,avedis3(maxc),clsdis(maxc,maxc)
        integer icls(maxd),ires(maxc,maxd),sz(maxc)
        character fcls*100,fcoor*100
        
        read(*,*) idim,ncc
        read(*,*) fcls
        read(*,*) fcoor
        do i=1,ncc
        sz(i)=0
        enddo
c         write(91,*) 'start'
c         close(91)
      
        open(unit=10,file=fcls,status='old') 
        nd=1
5        read(10,*,end=100) ii,icls(nd)
        sz(icls(nd))=sz(icls(nd))+1
        ires(icls(nd),sz(icls(nd)))=nd 
        nd=nd+1
        goto 5
100     close(10)
        nd=nd-1

c         write(90,*) nd,idim,ncc
        open(unit=11,file=fcoor,status='old') 
        do i=1,nd
        read(11,*) (ddd(i,l),l=1,idim)
        enddo
        close(11)

        do i=1,ncc
           do l=1,idim
           cent(i,l)=0
           enddo
           do j=1,sz(i)
            do l=1,idim
            cent(i,l)=cent(i,l)+ddd(ires(i,j),l)
            enddo
           enddo
          
           do l=1,idim
           cent(i,l)=cent(i,l)/sz(i)
c           write(*,*) i,l,cent(i,l)
           enddo

          avedis(i)=0
          avedis2(i)=0
          avedis3(i)=0   
           do j=1,sz(i)
            dd=0
            do l=1,idim
            dd=dd+(cent(i,l)-ddd(ires(i,j),l))**2
            enddo
            dd1=sqrt(dd)
            avedis2(i)=avedis2(i)+dd
            avedis(i)=avedis(i)+dd1
            if(avedis3(i).lt.dd1) avedis3(i)=dd1
           enddo
         enddo 
          do i=1,ncc-1
          do k=i+1,ncc
             dd=0
            do l=1,idim
            dd=dd+(cent(i,l)-cent(k,l))**2
            enddo
           clsdis(i,k)=sqrt(dd)
           clsdis(k,i)=sqrt(dd)
          enddo
          enddo




       write(*,*) 'ave dis ',(avedis(l)/sz(l),l=1,ncc)
       write(*,*) 'max dis ',(avedis3(l),l=1,ncc)
       write(*,*) 'cluster dis ',((l,k,clsdis(l,k),k=l+1,ncc),l=1,ncc-1)
c       write(*,*) 'within cluster SS ',(avedis2(l),l=1,ncc)

c        
c         write(*,*) i,sz(i), (ires(i,l),l=1,sz(i))
c        enddo
c        write(*,*) (ddd(1,l),l=1,idim)
         do i=1,nd-1
         do j=i+1,nd
          ss=0.0
          do k=1,idim
          ss=ss+(ddd(i,k)-ddd(j,k))**2
          enddo
          rd=sqrt(ss)
          pairdis(i,j)=rd
          pairdis(j,i)=rd
          enddo
         enddo
       

          sil=0.0
         do ii=1,nd
            ic=icls(ii)
            nc=sz(ic)
            sl(ii)=0
            xa(ii)=0
            xb(ii)=0
            if(nc.eq.1) goto 45            
            tmp=0
            do l=1,nc
            j=ires(ic,l)
            if(j.ne.ii) then
c            write(90,*) ii,j,l,nc
            tmp=tmp+pairdis(ii,j)
            endif
            enddo
            xa(ii)=tmp/(nc-1)
            tmp0=10000000.0
            do k=1,ncc
            if(k.ne.ic) then
            nc2=sz(k)
            tmp=0.0
            do l=1,nc2
            j=ires(k,l)
            tmp=tmp+pairdis(ii,j)
            enddo
            tmp=tmp/nc2
            if(tmp.lt.tmp0) tmp0=tmp
            endif           
            enddo
            xb(ii)=tmp0
             
            if(xa(ii).lt.xb(ii)) then
            sl(ii)=1-xa(ii)/xb(ii)
            else if ( xa(ii).eq.xb(ii)) then
            sl(ii)=0.0
            else
            sl(ii)=xb(ii)/xa(ii)-1
            endif 
45	  continue

c           write(90,*) ii,sl(ii),xa(ii),xb(ii)
           sil=sil+sl(ii)
         enddo



         write(*,*) 'SC ',sil/nd
c         write(92,*) 'stop ',sil/nd
c         close(92)
          




       stop
       end

