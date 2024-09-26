c
c      TESTING SUITE FOR STOKES KERNELS:  June 18th, 2020
c
c
      implicit none
      real *8 ztrg(3),source(3,20),targ(3,20),targ2(3,20)
      real *8 errgrad(3,20), pert(3), dmu(3), dnu(3)
      real *8 rlet(3), rvec(3)

      real *8, allocatable :: pot(:,:), pre(:), grad(:,:,:)
      real *8, allocatable :: pot2(:,:), pre2(:), grad2(:,:,:)

      real *8, allocatable :: potl(:,:), gradl(:,:,:)

      real *8, allocatable :: stoklet(:,:), strslet(:,:)
      real *8, allocatable :: strsvec(:,:)
      real *8, allocatable :: charge(:,:), dipvec(:,:,:)
      real *8, allocatable :: rotstr(:,:), rotvec(:,:)
      real *8, allocatable :: doubstr(:,:), doubvec(:,:)
      real *8 df,done,pi,h,thresh,errmax,gave1,gave2,gave3


      integer ipass(5)
      integer istress, irotlet, idoublet
      integer i,nd,ns,nt,ntest,nd1
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
      write(*,*) "=========================================="
      write(*,*) "Testing suite for stokkernels rotlet"

      open(unit=33,file='print_testres.txt',access='append')

      ntest = 4
      do i=1,ntest
        ipass(i) = 0
      enddo


      ns = 10
      nt = 12

      allocate(stoklet(3,ns),strslet(3,ns),strsvec(3,ns))
      allocate(rotstr(3,ns),rotvec(3,ns))
      allocate(doubstr(3,ns),doubvec(3,ns))

      do i = 1,ns
         source(1,i) = cos(i*done)
         source(2,i) = cos(10*i*done)
         source(3,i) = cos(100*i*done)
         stoklet(1,i) = 0.0d0
         stoklet(2,i) = 0.0d0
         stoklet(3,i) = 0.0d0
         strslet(1,i) = 0.0d0
         strslet(2,i) = 0.0d0
         strslet(3,i) = 0.0d0
         strsvec(1,i) = 0.0d0
         strsvec(2,i) = 0.0d0
         strsvec(3,i) = 0.0d0
         rotstr(1,i) = sin(4*i*done)
         rotstr(2,i) = sin(44*i*done)
         rotstr(3,i) = sin(444*i*done)
         rotvec(1,i) = sin(5*i*done)
         rotvec(2,i) = sin(55*i*done)
         rotvec(3,i) = sin(555*i*done)
         doubstr(1,i) =sin(4*i*done)
         doubstr(2,i) =sin(44*i*done)
         doubstr(3,i) =sin(444*i*done)
         doubvec(1,i) =sin(5*i*done)
         doubvec(2,i) =sin(55*i*done)
         doubvec(3,i) =sin(555*i*done)
      enddo

c      call prin2('src *',source,3*ns)
c      call prin2('stoklet *',stoklet,3*ns)
c      call prin2('strslet *',strslet,3*ns)
c      call prin2('strsvec *',strsvec,3*ns)
c

      h = 1d-5
      pert(1) = cos(1000*done)
      pert(2) = cos(2000*done)
      pert(3) = cos(3000*done)

      do i = 1,nt
         targ(1,i) = cos(12*i*done)
         targ(2,i) = cos(123*i*done)
         targ(3,i) = cos(1234*i*done)
         targ2(1,i) = targ(1,i) + h*pert(1)
         targ2(2,i) = targ(2,i) + h*pert(2)
         targ2(3,i) = targ(3,i) + h*pert(3)
      enddo

      allocate(pot(3,nt),pre(nt),grad(3,3,nt))
      allocate(pot2(3,nt),pre2(nt),grad2(3,3,nt))

c     test gradient of rotlet formula

      thresh = 1d-15
      nd1 = 1

      do i = 1,nt
         pre(i) = 0
         pot(1,i) = 0
         pot(2,i) = 0
         pot(3,i) = 0
         grad(1,1,i) = 0
         grad(2,1,i) = 0
         grad(3,1,i) = 0
         grad(1,2,i) = 0
         grad(2,2,i) = 0
         grad(3,2,i) = 0
         grad(1,3,i) = 0
         grad(2,3,i) = 0
         grad(3,3,i) = 0
         pre2(i) = 0
         pot2(1,i) = 0
         pot2(2,i) = 0
         pot2(3,i) = 0
         grad2(1,1,i) = 0
         grad2(2,1,i) = 0
         grad2(3,1,i) = 0
         grad2(1,2,i) = 0
         grad2(2,2,i) = 0
         grad2(3,2,i) = 0
         grad2(1,3,i) = 0
         grad2(2,3,i) = 0
         grad2(3,3,i) = 0
      enddo

      errmax = 0

      istress = 0
      irotlet = 1
      idoublet = 1
      call st3ddirectstokstrsrotdoubg(nd1,source,stoklet,istress,
     1     strslet,strsvec,irotlet,rotstr,rotvec,
     2     idoublet,doubstr,doubvec,
     3     ns,targ,nt,pot,pre,grad,thresh)
      call st3ddirectstokstrsrotdoubg(nd1,source,stoklet,istress,
     1     strslet,strsvec,irotlet,rotstr,rotvec,
     2     idoublet,doubstr,doubvec,
     3     ns,targ2,nt,pot2,pre2,grad2,thresh)


      do i = 1,nt
         df = pot2(1,i) - pot(1,i)
         gave1 = (grad(1,1,i) + grad2(1,1,i))/2.0d0
         gave2 = (grad(2,1,i) + grad2(2,1,i))/2.0d0
         gave3 = (grad(3,1,i) + grad2(3,1,i))/2.0d0         
         errgrad(1,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot2(2,i) - pot(2,i)
         gave1 = (grad(1,2,i) + grad2(1,2,i))/2.0d0
         gave2 = (grad(2,2,i) + grad2(2,2,i))/2.0d0
         gave3 = (grad(3,2,i) + grad2(3,2,i))/2.0d0         
         errgrad(2,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         df = pot2(3,i) - pot(3,i)
         gave1 = (grad(1,3,i) + grad2(1,3,i))/2.0d0
         gave2 = (grad(2,3,i) + grad2(2,3,i))/2.0d0
         gave3 = (grad(3,3,i) + grad2(3,3,i))/2.0d0         
         errgrad(3,i) = abs(1.0d0 - h*(gave1*pert(1) + gave2*pert(2)
     1        + gave3*pert(3))/df)

         errmax = max(errmax,errgrad(1,i))
         errmax = max(errmax,errgrad(2,i))
         errmax = max(errmax,errgrad(3,i))                  
         
      enddo

      call prin2('err grad (just rotlet)*',errgrad,3*nt)

      if (errmax .lt. 1d-6) ipass(1) = 1

      stop
      end
c
c
c
c
c
c
      subroutine errprint(pot,opot,fld,ofld,errs)
      implicit real *8 (a-h,o-z)
      real *8 pot,opot,fld(3),ofld(3)
      real *8 errs(2)
 1000  format(4D15.5)
      err = 0
      ddd = 0
      err = err + abs(fld(1)-ofld(1))**2
      err = err + abs(fld(2)-ofld(2))**2
      err = err + abs(fld(3)-ofld(3))**2
      ddd = ddd + abs(ofld(1))**2
      ddd = ddd + abs(ofld(2))**2
      ddd = ddd + abs(ofld(3))**2
      err = sqrt(err)
      ddd = sqrt(ddd)

      err1 = abs(pot-opot)/abs(opot)
      write(*,'(a,e11.4,a,e11.4)')
     1     'pot error=',err1,'   grad error=',err/ddd

      errs(1) = err1
      errs(2) = err/ddd

      return
      end
c
