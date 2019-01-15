c
c
c

      subroutine hfmm3dpartstoscp_vec(nd,eps,zk,nsource,source,
     1    charge,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstoscg_vec(nd,eps,zk,nsource,source,
     1    charge,pot,grad)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine hfmm3dpartstosdp_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstosdg_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,pot,grad)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstoscdp_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,1)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c
c
c
c
c

      subroutine hfmm3dpartstoscdg_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,grad)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,1)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,1),gradtarg(nd,3,1)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      ntarg = 0

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine hfmm3dpartstotcp_vec(nd,eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine hfmm3dpartstotcg_vec(nd,eps,zk,nsource,source,
     1    charge,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdp_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotdg_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstotcdp_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine hfmm3dpartstotcdg_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,1),grad(nd,3,1)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcp_vec(nd,eps,zk,nsource,source,
     1    charge,pot,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcg_vec(nd,eps,zk,nsource,source,
     1    charge,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,1)
      double complex dipvec(nd,3,1)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
c

      subroutine hfmm3dpartstostdp_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,pot,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostdg_vec(nd,eps,zk,nsource,source,
     1    dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,1)
      double complex dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdp_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,ntarg,targ,pottarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource)
      double complex pottarg(nd,ntarg)

      double complex grad(nd,3),gradtarg(nd,3)
      double complex hess(nd,6),hesstarg(nd,6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c

      subroutine hfmm3dpartstostcdg_vec(nd,eps,zk,nsource,source,
     1    charge,dipstr,dipvec,pot,grad,ntarg,targ,pottarg,
     2    gradtarg)
      implicit none
      double precision eps
      double complex zk

      integer nsource,ntarg,ifcharge,ifdipole,ifpgh,ifpghtarg
      integer nd
      
      double precision source(3,nsource),targ(3,ntarg)
      double complex charge(nd,nsource),dipstr(nd,nsource)
      double complex dipvec(nd,3,nsource)

      double complex pot(nd,nsource),grad(nd,3,nsource)
      double complex pottarg(nd,ntarg),gradtarg(nd,3,ntarg)

      double complex hess(6),hesstarg(6)

      
      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm3dpart(nd,eps,zk,nsource,source,ifcharge,charge,
     1      ifdipole,dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,
     2      ifpghtarg,pottarg,gradtarg,hesstarg)

      return
      end
c
c
c
c
c
