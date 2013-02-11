           program anneal
          real xs(10000),ys(10000)
          real xt(20000),yt(20000)
          integer itype(20000),itypet(20000)
          integer iweight(20000),iweightt(20000)
          real xtt(20000),ytt(20000)
          integer indexyt(20000)
          integer stp(10000,100),tsp(20000,10)
          integer sta(10000,10)
          integer tsa(20000,2)
          integer ttex(20000,100),tsex(20000,100,3)
          common/nblk/nx,ny
          integer ntt(3)
          logical lcross

          write(*,*)'number of passes?'
          read(*,*)np

          write(*,*)'nx,ny=?'
          read(*,*)nx,ny

          write(*,*)'log Kfact=?'
          read(*,*)lkfact
          kfact=10**lkfact

          write(*,*)'patrol radius/pitch?'
          read(*,*)patr
          patr2=patr**2

          write(*,*)'exclusion radius/pitch?'
          read(*,*)excl
          dist2ex=excl**2


          ns=nx*ny
          sqrt3by2=sqrt(3.)/2
          do i=1,ns
            xs(i) = mod(i-1,nx) - mod(i,2)/2
            ys(i) = sqrt3by2 * int((i-1)/nx)
            write(32,*)i,xs(i),ys(i)
          enddo

          ntt(1)=0
          ntt(2)=0
          ntt(3)=0

          write(*,*)'random (0) or stephs test target set (1)?'
          read(*,*)irand
          if(irand.eq.0)then
            nt=np*ns
            do j=1,nt
              xtt(j)=ran1(iseed) * nx
              ytt(j)=ran1(iseed) * ny * sqrt3by2
              itypet(j)=1
              iweightt(j)=1
            enddo
            ntt(1)=nt

          elseif(irand.eq.1)then
            open(10,file='elg_test_sample.dat',status='old')
            do j=1,20000
              read(10,*,end=11)ai,aimag,aimage,z,atype,size,itarget,x,y
              xtt(j)=(y+0.6)*nx/1.2     ! so rotate and rescale to fit nx pitches
              ytt(j)=(x+0.55)*nx/1.2   
              itypet(j)=2
              iweightt(j)=1
            enddo
 11         close(10)
            ntt(2)=j-1
            write(*,*)ntt(2),' elg targets'

            write(*,*)'mag lim for LRGs?'
            read(*,*)aimaglim
            open(10,file='lrg_test_sample.dat',status='old')
            do jt=i,20000
              read(10,*,end=12)ai,aimag,aimage,z,it,size,itarget,x,y
              if(aimag.lt.aimaglim)then
                xtt(j)=(y+0.6)*nx/1.2     ! so rotate and rescale to fit 40 pitches
                ytt(j)=(x+0.55)*nx/1.2  
                itypet(j)=3
                iweightt(j)=1
                j=j+1
              endif
            enddo
 12         close(10)
            ntt(3)=j-1-ntt(2)
            write(*,*)ntt(3),' lrg targets'
            nt=j-1
           endif
          write(*,*)ns,nt,np
          write(*,*)'iseed?'
          read(*,*)iseed
c          iseed=7481885

 
          call indexx(nt,ytt,indexyt)
          do j=1,nt
            xt(j)=xtt(indexyt(j))
            yt(j)=ytt(indexyt(j))
            itype(j)=itypet(indexyt(j))
            iweight(j)=iweightt(indexyt(j))
            write(33,*)j,xt(j),yt(j),itype(j),iweight(j)
          enddo

          write(*,*)'field set up'

C  set up close pairs
c          dist2ex=(0.7/6)**2    !700um, with pitch=6mm
          jtmax=0
          npp=0
          do j1=1,nt
           jt=0
           call hunt(yt,nt,yt(j1)-excl,j2start)
           do j2=max(j2start-1,1),nt
            if(yt(j2).gt.yt(j1)+excl)goto 2   ! done with j1
            if((j1.ne.j2).and.
     !       (dist2(xt(j1),yt(j1),xt(j2),yt(j2)).lt.dist2ex))then
              npp=npp+1
              jt=jt+1
              ttex(j1,jt)=j2
              jtmax=max(jtmax,jt)
            endif
           enddo
2          if(mod(j1,1000).eq.0)write(*,*)j1,jtmax,npp/2
          enddo
          write(*,*)jtmax,npp/2,npp/(2.*nt)

C    set up list stp of all target possibilities for each spine and vice versa
          ntot=0
          itmax=0
          do i=1,ns
           it=0
           call hunt(yt,nt,ys(i)-patr,jstart)
           do j=max(jstart-1,1),nt
            if(yt(j).gt.ys(i)+patr)goto 3   ! done with i
            if(dist2(xs(i),ys(i),xt(j),yt(j)).lt.patr2)then   ! configurable
               it=it+1
               stp(i,it)=j
               ntot=ntot+1
               itmax=max(itmax,it)
               jt=1
               do while(tsp(j,jt).ne.0)
                 jt=jt+1
               enddo
               tsp(j,jt)=i
             endif
           enddo
3          if(mod(i,1000).eq.0)write(*,*)i,itmax,ntot
          enddo
          write(*,*)itmax,ntot

C    set up list tsex of all possibile spine-target crossings
         nx=0
         nxtmax=0
         do j1=1,nt
           call hunt(yt,nt,yt(j1)-2*patr,j2start)
           do j2=max(j2start-1,1),nt
            if(yt(j2).gt.yt(j1)+2*patr)goto 4   ! done with j1
            if((j1.ne.j2).and.
     !       (dist2(xt(j1),yt(j1),xt(j2),yt(j2)).lt.4*patr2))then  !candidates for crossings
              nxt=0
              jt1=1
              do while(tsp(j1,jt1).ne.0)
                jt2=1
                do while(tsp(j2,jt2).ne.0)
                  if(dist2(xs(tsp(j1,jt1)),ys(tsp(j1,jt1)),
     !              xs(tsp(j2,jt2)),ys(tsp(j2,jt2))).lt.4*patr2)then 
                     if(lcross(xs(tsp(j1,jt1)),ys(tsp(j1,jt1)),
     !               xs(tsp(j2,jt2)),ys(tsp(j2,jt2)),
     !               xt(j1),yt(j1),xt(j2),yt(j2)))then

                      nxt=nxt+1
                      nxtmax=max(nxtmax,nxt)
                      nx=nx+1
                      write(34,*)xs(tsp(j1,jt1)),ys(tsp(j1,jt1)),
     !                xt(j1),yt(j1),
     !                xs(tsp(j2,jt2)),ys(tsp(j2,jt2)),
     !                xt(j2),yt(j2),ak1,ak2
                      tsex(j1,nxt,1)=tsp(j1,jt1)
                      tsex(j1,nxt,2)=j2
                      tsex(j1,nxt,3)=tsp(j2,jt2)
                    endif
                  endif
                  jt2=jt2+1
                enddo
                jt1=jt1+1
               enddo
             endif
           enddo
4          if(mod(j1,1000).eq.0)write(*,*)j1,nx,nxtmax
          enddo
          write(*,*)nx,nxtmax

C do the annealing

c          kmax=10000000
           kmax=np*ns*kfact
          dak=10./kmax
          nall=0
          do k=1,kmax
           i=1+int(ran1(iseed)*np*ns)
           ii=1+i/ns        ! ii is pass number
           i=mod(i,ns)       ! i is fiber number
           if(sta(i,ii).eq.0)then         !fiber not allocated, put on a target if possible
            jj=0
            ranmax=0
            jt=1
            do while (stp(i,jt).ne.0)
             if(tsa(stp(i,jt),1).eq.0)then !available unfibered target
              jtt=1
              do while(ttex(stp(i,jt),jtt).ne.0) !check not too near an already configured target
               if((tsa(ttex(stp(i,jt),jtt),1).ne.0).and.
     !          (tsa(ttex(stp(i,jt),jtt),2).eq.ii))goto 1! reject, but only if in same pass
               jtt=jtt+1
              enddo                 !not a close pair
              itt=1
              do while(tsex(stp(i,jt),itt,1).ne.0)   !check for xings
               if((tsex(stp(i,jt),itt,1).eq.i).and.   !this target/spine combo
     !         (tsa(tsex(stp(i,jt),itt,2),2).eq.ii).and.   !other target configured in this pass
     !        (tsa(tsex(stp(i,jt),itt,2),1).eq.tsex(stp(i,jt),itt,3)))  !and the other pair is also in use
     !          goto 1 
                itt=itt+1
               enddo

               rant=ran1(iseed)      ! choose one of ok targets at random
               if(rant.gt.ranmax)then
                ranmax=rant
                jj=stp(i,jt)
               endif
             endif
  1          jt=jt+1
            enddo
            if(jj.ne.0)then   !update the allocation
             sta(i,ii)=jj
             tsa(jj,1)=i
             tsa(jj,2)=ii   
             nall=nall+1
            endif
           elseif(iweight(sta(i,ii))*ran1(iseed).lt.exp(-k*dak))then     !deallocate fiber 
            tsa(sta(i,ii),1)=0
            tsa(sta(i,ii),2)=0
            sta(i,ii)=0
            nall=nall-1
           endif
           if(mod(k,kmax/100).eq.0)
     !      write(*,*)k,nall,nall/(np*ns*1.),nall/real(nt)
          enddo 
          
          do it=1,3
            nallt=0
            do j=1,nt
              if((itype(j).eq.it).and.(tsa(j,1).ne.0))then
                nallt=nallt+1
                write(10+it,*)xt(j),yt(j),itype(j),iweight(j)
              elseif((itype(j).eq.it).and.(tsa(j,1).eq.0))then
                write(20+it,*)xt(j),yt(j),itype(j),iweight(j)
              endif
            enddo
            write(*,*)it,ntt(it),nallt,nallt/max(1.,real(ntt(it)))
           enddo
          write(*,*)nall,nall/real(nt),nall/real(np*ns), 
     !      nall/sqrt(1.*nt*np*ns) 
          end
C---------------------------------------------
          function dist2(x1,y1,x2,y2)
          common/nblk/nx,ny

          dist2=min((x1-x2)**2,(abs(x1-x2)-nx)**2)+ 
     !    min((y1-y2)**2,((abs(y1-y2)-ny)*sqrt3by2)**2)

          return
          end
C-------------------------------------------------
          function lcross(xs1,ys1,xs2,ys2,xt1,yt1,xt2,yt2)
          logical lcross
          common/nblk/nx,ny
           lcross=.false.
           ax=xs2-xs1
           if(ax.gt.nx/2)ax=ax-nx
           if(ax.lt.-nx/2)ax=ax+nx
           ay=ys2-ys1
           if(ay.gt.ny/2)ay=ay-ny
           if(ay.lt.-ny/2)ay=ay+ny
           b1x=xt1-xs1
           if(b1x.gt.ny/2)b1x=b1x-nx
           if(b1x.lt.-ny/2)b1x=b1x+nx
           b1y=yt1-ys1
           if(b1y.gt.ny/2)b1y=b21-ny
           if(b1y.lt.-ny/2)b1y=b1y+ny
           b2x=xt2-xs2
           if(b2x.gt.ny/2)b2x=b2x-nx
           if(b2x.lt.-ny/2)b2x=b2x+nx
           b2y=yt2-ys2
           if(b2y.gt.ny/2)b2y=b2y-ny
           if(b2y.lt.-ny/2)b2y=b2y+ny
           det=(b1x*b2y) - (b1y*b2x)
           ak1= ((b2y*ax) - (b2x*ay))/det
           ak2= ((b1y*ax) - (b1x*ay))/det
           if ((ak1.lt.1).and.(ak1.gt.0).and.
     !      (ak2.lt.1).and.(ak2.gt.0))lcross=.true.  !crossing
          return
          end
C-----------------------------------------------
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C------------------------------------------------ 
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C------------------------------------------------------------
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END

C---------------------------------
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C------------------------------------------------------
