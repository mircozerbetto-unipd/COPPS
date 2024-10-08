c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c       POWELL
c
c       Minimizes a function func of n variables, using Powell
c       method, starting from initial direction xi(n,n)
c
c       p(np)     : starting point (I) / minimum (O)
c       xi(np,np) : initial direction (I) / current direction (O)
c       n         : logical dimensions (I)
c       np        : physical dimensions (I)
c       ftol      : fractional tolerance (I)
c       iter      : current number of iterations (O)
c       fret      : current value of func at p (O) 
c
c       Uses: linmin / f1dim / mnbrak / brentc
c
c       Numerical Recipes, AA VV (Cambridge, 1990), pp. 294-300
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        subroutine powell(p,xi,n,np,ftol,atol,iter,fret,iflag)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        implicit real*8 (a-h,o-z)

        include 'iounits.inc'

        external func
        parameter (nmax=20,itmax=200)
        dimension p(np),xi(np,np)
        dimension pt(nmax),ptt(nmax),xit(nmax)

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        fret=func(p)
c-------MIRCO ZERBETTO---02/04/2012|
        if (fret < atol) then
          iflag = 0
          return
        end if
c----------------------------------|

        do j=1,n
           pt(j)=p(j)
        end do

        iter=0
1       iter=iter+1
        fp=fret
        ibig=0
        del=0.d0
        do i=1,n
           do j=1,n
              xit(j)=xi(j,i)
           end do
           fptt=fret
           call linmin(p,xit,n,fret)
           if(dabs(fptt-fret).gt.del) then
              del=dabs(fptt-fret)
              ibig=i
           endif
        end do

        ratio=2.d0*dabs(fp-fret)/(dabs(fp)+dabs(fret))
        diffr=dabs(fp-fret)
        if (ratio.le.ftol) then
           iflag=0
           return
        else if (diffr.lt.atol) then
           iflag=1
           return
        else
           write (consunit,"('+ iter n. ',i4,'/ rel = ',
     *          g12.6,'/ abs =',g12.6)") 
     *          iter,ratio,diffr
           write (logunit,"('+ iter n. ',i4,'/ rel = ',
     *          g12.6,'/ abs =',g12.6)") 
     *          iter,ratio,diffr

        end if
           
        if (iter.eq.itmax) then
c-------MIRCO ZERBETTO---02/04/2012|
           write(*,*)
           write(*,*) "Max iterations (", itmax,") reached in powell
     *                 routine."
           write(*,*)
c----------------------------------|
           iflag=-1
           return
        end if

        do j=1,n
          ptt(j)=2.d0*p(j)-pt(j)
          xit(j)=p(j)-pt(j)
          pt(j)=p(j)
        end do
        fptt=func(ptt)
        if(fptt.ge.fp) go to 1
        t=2.d0*(fp-2.d0*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
        if(t.ge.0.d0) go to 1
        call linmin(p,xit,n,fret)
        do j=1,n
          xi(j,ibig)=xi(j,n)
          xi(j,n)=xit(j)
        end do

        go to 1

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        subroutine linmin(p,xi,n,fret)

        implicit real*8 (a-h,o-z)
        parameter(nmax=200,tol=1.d-04)
        dimension p(n),xi(n)
        dimension pcom(nmax),xicom(nmax)
        common/f1com/pcom,xicom,ncom
        external f1dim

        ncom=n
        do j=1,n
          pcom(j)=p(j)
          xicom(j)=xi(j)
        end do
        ax=0.d0
        xx=1.d0
        call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
        fret=brent(ax,xx,bx,f1dim,tol,xmin)
        do j=1,n
          xi(j)=xmin*xi(j)
          p(j)=p(j)+xi(j)
        end do

        return
        end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        function f1dim(x)
        implicit real*8 (a-h,o-z)
        parameter (nmax=200)
        dimension pcom(nmax),xicom(nmax),xt(nmax)
        common/f1com/pcom,xicom,ncom

        do j=1,ncom
          xt(j)=pcom(j)+x*xicom(j)
        end do
        f1dim=func(xt)

        return
        end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        subroutine mnbrak(ax,bx,cx,fa,fb,fc,func)
        implicit real*8 (a-h,o-z)
        external func
        parameter(gold=1.618034d0,glimit=100.d0,tiny=1.d-20)
        
        fa=func(ax)
        fb=func(bx)
        if(fb.gt.fa) then
           dum=ax
           ax=bx
           bx=dum
           dum=fb
           fb=fa
           fa=dum
        endif
        cx=bx+gold*(bx-ax)
        fc=func(cx)
 1      if(fb.ge.fc) then
           r=(bx-ax)*(fb-fc)
           q=(bx-cx)*(fb-fa)
           u=bx-((bx-cx)*q-(bx-ax)*r)/
     *          (2.d0*dsign(dmax1(dabs(q-r),tiny),q-r))
           ulim=bx+glimit*(cx-bx)
           if((bx-u)*(u-cx).gt.0.d0) then
              fu=func(u)
              if(fu.lt.fc) then
                 ax=bx
                 fa=fb
                 bx=u
                 fb=fu
                 return
              else if(fu.gt.fb) then
                 cx=u
                 fc=fu
                 return
              endif
              u=cx+gold*(cx-bx)
              fu=func(u)
           else if((cx-u)*(u-ulim).gt.0.d0) then
              fu=func(u)
              if(fu.lt.fc) then
                  bx=cx
                  cx=u
                  u=cx+gold*(cx-bx)
                  fb=fc
                  fc=fu
                  fu=func(u)
              endif
           else if((u-ulim)*(ulim-cx).ge.0.d0) then
              u=ulim
              fu=func(u)
           else
              u=cx+gold*(cx-bx)
              fu=func(u)
           endif
           ax=bx
           bx=cx
           cx=u
           fa=fb
           fb=fc
           fc=fu
           go to 1
        endif
        
        return
        end

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        function brent(ax,bx,cx,func,tol,xmin)
        implicit real*8 (a-h,o-z)
        external func
        parameter(itmax=2000,cgold=.381966d0,zeps=1.d-10)
        
        a=dmin1(ax,cx)
        b=dmax1(ax,cx)
        v=bx
        w=v
        x=v
        e=0.d0
        fx=func(x)
        fv=fx
        fw=fx

        do iter=1,itmax
           xm=0.5d0*(a+b)
           tol1=tol*dabs(x)+zeps
           tol2=2.d0*tol1
           if(dabs(x-xm).le.(tol2-0.5d0*(b-a))) go to 3
           if(dabs(e).gt.tol1) then
              r=(x-w)*(fx-fv)
              q=(x-v)*(fx-fw)
              p=(x-v)*q-(x-w)*r
              q=2.d0*(q-r)
              if(q.gt.0.d0) p=-p
              q=dabs(q)
              etemp=e
              e=d
              if(dabs(p).ge.dabs(0.5d0*q*etemp).or.p.le.
     *        q*(a-x).or.p.ge.q*(b-x)) go to 1
              d=p/q
              u=x+d
              if(u-a.lt.tol2.or.b-u.lt.tol2) d=dsign(tol1,xm-x)
              go to 2
           endif
 1         if(x.ge.xm) then
              e=a-x
           else
              e=b-x
           endif
           d=cgold*e
 2         if(dabs(d).ge.tol1) then
              u=x+d
           else
              u=x+dsign(tol1,d)
           endif
           fu=func(u)
           if(fu.le.fx) then
              if(u.ge.x) then
                 a=x
              else
                 b=x
              endif
              v=w
              fv=fw
              w=x
              fw=fx
              x=u
              fx=fu
           else
              if(u.lt.x) then
                 a=u
              else
                 b=u
              endif
              if(fu.lt.fw.or.w.eq.x) then
                 v=w
                 fv=fw
                 w=u
                 fw=fu
              else if(fu.le.fv.or.v.eq.x.or.v.eq.w) then
                 v=u
                 fv=fu
              endif
           endif
        end do
c       pause ' brent exceed maximum iterations'
        write (*,*) ' brent exceed maximum iterations'
 3      xmin=x
        brent=fx

        return
        end


c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
