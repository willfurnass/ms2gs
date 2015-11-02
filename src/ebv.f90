! 2015-01-28 11:26:19 : include q as input for vanradenB
! 2015-01-26 11:32:07 : so far I include m=2000 for vanradenB
! 2015-01-21 14:39:14 : incorporates Vanraden's approx non linear approach B (AL)
! 2014-08-11 15:17:41 : new conv criterion
! 2014-08-07 14:33:30 : error with miss in blup discovered by Andres L.
! 2014-08-04 14:39:14 : incorporates Vanraden's approx non linear approach A (AL)
! 15.04.2014 13:51:17 : accounts for missing
! 10/04/2014 13:17:04 : error in computing freqs after standardizing X
! This module should contain subroutines associated with EBV estimation: GS or blup

! this v of solveSNP works as a subroutine, returns EBVs & vars
! input:  id(i), y(i), X(i,1:nsnp+1), va_prior, ve_prior
! output: ebv(i)=Xb, b(), va, ve
! USAGE
!   call solveSNP (ebv, sol, va, ve, X, y, vce, maxiter, f)
! WARNING
!   to minimize memory usage, first column in X must be '1', next are genotypes


!-----------------
 module solveByPCG
!-----------------
 use kinds

 contains

  subroutine pcgru(X,y,D,vare,sol,maxiter)
  ! preconditioned conjugated gradients

  implicit none
   real(r8), parameter :: minConv=1e-14
   real(r8)::X(:,:),D(:),sol(:)
   ! solve A sol = rhs by preconditioned conjugate gradient for densem A
   real (r8)::y(:)

   real (r8),allocatable::m(:),r(:),p(:),z(:),w(:),rhs(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val,vare
   real(r8):: errold(size(y)) !!AL
   integer::i,j,k,n,ndata,top
   logical:: verbose
   integer,optional:: maxiter

   top=100000
   if(present(maxiter)) top=maxiter

   n=size(sol);ndata=size(y)
   allocate(m(size(sol)))

   allocate (r(n),p(n),z(n),w(n),rhs(n))
   ! find the preconditioner
   !diagonal preconditioner
      m=0
      do i=1,n
         val=dot_product(x(:,i),x(:,i))/vare+d(i)
         if (val /= 0) m(i)=1d0/val
      enddo

   ! build rhs
   rhs=matmul(transpose(X),y)/vare
   sol=0d0
   r=rhs

   do k=1,top
      errold=y-matmul(X,sol) !!AL
      z=m*r
      tau=dot_product(z,r)
      if (k == 1) then
         beta=0; p=z
       else
         beta=tau/oldtau
         p=z+beta*p
      end if
      !w=matmul(a,p)
      w=matmul(transpose(X),matmul(X,p))/vare+D*p ! if id vare
      alpha=tau/dot_product(p,w)
      sol=sol+alpha*p
      if (mod(k,100) /= 0) then
           r=r-alpha*w
         else
           r=rhs - matmul(transpose(X),matmul(X,sol))/vare - D*sol
      endif
      !conv=(alpha**2)*dot_product(p,p)/dot_product(sol,sol)
      !conv=dot_product(r,r)/dot_product(rhs,rhs)
      !if((mod(k,1)==0))  print*,'round ',k,'   convergence=',conv
      conv=sum( (errold- (y-matmul(X,sol)))**2 )/dot_product(errold,errold) !!AL
      if (conv < minConv) exit
      oldtau=tau
   enddo

   deallocate (m,r,p,z,w)
   !print*,k,' iterations,   convergence criterion=',conv
   if(conv > minConv) print*, 'ULL: convergence not attained pcgru!'

  end subroutine

!----------
 end module
!----------

!--------------------
 MODULE EBV
!--------------------
! This module should contain subroutines associated with EBV estimation: GS or blup
 use kinds
 use random
 use stat
 use solveByPCG
 implicit none

 CONTAINS

!------------------
 subroutine do_Ainv (Ainv, ped)
!------------------
! old bluppers never die
! Ignacy's A inverse implementation
 integer   :: ped(:,:), p(3), par_stat, i, k, l
 real(r8)  :: Ainv(:,:)
 real      :: w(3)    = (/1., -.5, -.5/), &
              res(4)  = (/2., 1.333, 1., 0./)

 Ainv(:,:) = 0
 do i=1, size(ped,1)
    !--> ind + parents' ids
    p(1:3)= ped(i,:)
    par_stat = 1
    if(p(2)==0) par_stat = par_stat+1
    if(p(3)==0) par_stat = par_stat+1
    do k=1,3
       do l=1,3
          if (p(k) /= 0 .and. p(l) /=0) then
             Ainv(p(k),p(l)) = Ainv(p(k),p(l)) + w(k) * w(l) * res(par_stat)
          endif
       enddo
    enddo
 enddo
!-----------------------
 end subroutine
!-----------------------

!-----------------
 subroutine pcgrua (sol,X,y0,Ai,varg,vare,miss,maxiter)
!-----------------
! andres legarra
! makeup by mpe
! BLUP by preconditioned conjugated gradients
! Ai is A-inverse, possibly obtained using Henderson's rules
! solve (X'X/vare + Ai/varg) sol = rhs by preconditioned conjugate gradient for densem LHS
 real(r8), parameter :: minConv=1e-14

 real(r8) :: X(:,:),Ai(:,:),sol(:),y0(:), y(size(y0))
 real (r8),allocatable::m(:),r(:),p(:),z(:),w(:),rhs(:),Sigma(:,:)

 real (r8)::alpha,beta,tau,oldtau,conv,val,vare,varg
 integer:: miss, i,j,k,n,ndata,top=100000, offset
 logical:: verbose
 integer,optional :: maxiter

 if(present(maxiter)) top=maxiter

 n=size(sol)
 ndata=size(y)
 y = y0
 where(int(y).eq.miss) y=0.

 offset=n-size(Ai,1)
 ! this thing assumes that polygenic solutions are the last ones
 !---------
 allocate(Sigma(n,n))
 ! inverse covariance matrix of ALL unknowns in the system, this is easier to fit into the PCG
 ! typically [0 0]
 !           [0 Ai]
 Sigma=0
 Sigma(offset+1:n,offset+1:n)=Ai

 allocate(m(n),r(n),p(n),z(n),w(n),rhs(n))
 ! find the preconditioner
 !diagonal preconditioner
 m=0
 do i=1,n
    val=dot_product(x(:,i),x(:,i))/vare+Sigma(i,i)/varg
    if (val /= 0) m(i)=1d0/val
 enddo

 ! build rhs
 rhs=matmul(transpose(X),y)/vare
 sol=0d0
 r=rhs

 do k=1,top
    z=m*r
    tau=dot_product(z,r)
    if (k == 1) then
       beta=0; p=z
    else
       beta=tau/oldtau
       p=z+beta*p
    end if
    !w=matmul(a,p)
    w=matmul(transpose(X),matmul(X,p))/vare+matmul(Sigma,p)/varg ! if id vare
    alpha=tau/dot_product(p,w)
    sol=sol+alpha*p
    if (mod(k,100) /= 0) then
       r=r-alpha*w
    else
       r=rhs - matmul(transpose(X),matmul(X,sol))/vare - matmul(Sigma,sol)/varg
    endif
    !conv=(alpha**2)*dot_product(p,p)/dot_product(sol,sol)
    conv=dot_product(r,r)/dot_product(rhs,rhs)
    !if((mod(k,1)==0))  print*,'round ',k,'   convergence=',conv
    if (conv < minConv) exit
    oldtau=tau
 enddo

 deallocate (m,r,p,z,w,rhs,sigma)
 if(conv > minConv) print*, 'ULL: convergence not attained in pcgrua!'

!--------------
 end subroutine
!--------------

!----------------------
 subroutine solveSNPsub0 (ebv, sol, varg, vare, X, y, f, dovce, dopcg, mx_iter)
!----------------------
! SNP-RR by Andres Legarra (Andres.Legarra@toulouse.inra.fr)
! reshaped to subroutine by M. Perez-Enciso (miguel.perez@uab.es)
! allows for mssing values by removing those in X and y

 real(r8) :: ebv(:), sol(:), X(:,:), y(:), varg, vare, f(:)
 logical  :: dovce, dopcg
 integer  :: mx_iter
 optional :: dovce, mx_iter, dopcg, f

 integer  :: ndata, nsnp, neq, maxiter, nburning
 integer  :: i, j, io, iter, n, ioin=1
 real(r8) :: a, vara, mu, eps, val, lhs, rhs, &
             vare_ap, vara_ap, df_vare, df_vara, twosumpq
 real(r8),allocatable:: e(:), D(:), freq(:), ssol(:), xpx(:)
 character:: ofile*20='solutions'
 logical  :: VCE, wantPCG

 !--> inits
 neq      = size(X,dim=2)
 nsnp     = neq - 1
 ndata    = size(y)
 nburning = 10000
 maxiter  = 100000
 vce      = .false. !--> VarianceComponentEstimation
 wantPCG  = .true.  !--> conjugate gradient
 if(present(mx_iter)) maxiter = mx_iter
 if(present(dopcg))   wantPCG = dopcg
 if(present(dovce))   vce     = dovce

 ! Prior information for variances
 ! NOTE the 'chin' invertedchi square subroutine is parameterized in terms of
 ! degrees of freedom and expectation of the random variable and NOT
 ! in terms of d.f. and scale.
 ! That is: E(chin(S,v))=S and not to S*v
 ! Hence the way of writing the d.f. and expectations in a certain manner
 vare_ap = 9d0 !--> simulated was 10
 vara_ap = 10/(2*nsnp*0.5*0.15) !--> for 10 of genetic variance and 50000 SNPs with p=0.15
 df_vare = 4 ! low informative
 df_vara = 4

 allocate (xpx(neq),e(ndata),D(neq),freq(nsnp),ssol(neq))
 !vare=10.0
 !vara=10/(2*nsnp*.25)      !2sumpq for 0.5

 !--> computes freq if not provided, could be important for seq data
 if(.not.present(f)) then
    do i=1,nsnp
       freq(i)=real(sum(X(:,i+1)))/(2.*ndata)
    enddo
 else
    freq = f
 endif
 twosumpq = 2*sum(freq*(1-freq))
 !--> var per marker
 vara = varg/twosumpq
 e    = y
 sol  = 0d0
 ssol = 0d0
 D    = 1/vara
 D(1) = 0 ! in fact D=G^{-1}
 !print *,'vara varap vare varep:', vara, vara_ap, vare, vare_ap

 !--> genotype centering
 X(:,2:neq)=X(:,2:neq)-1d0

 if (wantPCG) then
    call pcgru(X,y,D,vare,sol,maxiter)
 else
    do i=1,neq
       xpx(i)=dot_product(X(:,i),X(:,i)) !form diagonal of X′X
    enddo
    !--> Gibbs sampling to get var estimates, slooow
    do iter=1,maxiter
       eps=0
       do i=1,neq
          !form lhs
          lhs=xpx(i)/vare+D(i)
          ! form rhs with y corrected by other effects (formula 1)
          rhs=dot_product(X(:,i),e)/vare +xpx(i)*sol(i)/vare
          ! do Gauss Seidel
          val=rhs/lhs
          ! MCMC sample solution from its conditional (commented out here)
          if(VCE) val=normal(val,1d0/lhs)
          eps=eps+(val-sol(i))**2
          ! update e with current estimate (formula 2)
          e=e - X(:,i)*(val-sol(i))
          ! update sol
          sol(i)=val
       enddo
       if (.not.VCE)then
          !print *,iter,'eps: ',eps
          if(eps<10d-12) EXIT
       else !--> VCE
          vara = (dot_product(sol(2:neq),sol(2:neq))+ df_vara*vara_ap) / (nsnp+df_vara)
          vara = chin(vara, nsnp+df_vara)
          vare = (dot_product(e,e) + df_vare*vare_ap) / (ndata+df_vare)
          vare = chin(vare, ndata+df_vare)
          !print*,'vara,vare',vara,vare
       endif
       if(iter>nburning) ssol=ssol+sol
    enddo !--> end iter
    if(VCE) sol=ssol/(maxiter-nburning)
 endif

 ebv(:) = matmul(X,sol)

 deallocate(xpx,e,D,freq,ssol)

!--------------
 end subroutine
!--------------

!-----------------------
 subroutine solveSNPsub (ebv, sol, varg, vare, X0, y0, f, dovce, dopcg, mx_iter, miss, nonl, nB)
!-----------------------
! AL, 21/01/2015 : nonl is character, if 0 nothing is done, 'A' = nonlinearA, 'B'=
! nonlinear B
! I monitor best missing values and fixed snps ....
! SNP-RR by Andres Legarra (Andres.Legarra@toulouse.inra.fr)
! Reshaped to subroutine by M. Perez-Enciso (miguel.perez@uab.es)
! Allows for mssing values by removing rows in X and y (mpe)
! incorporates Vanraden's approx non linearA and nonlinearB approaches
! the constant q is a putative number of QTLs, could be fixed at the actual
! number of QTLs or yet something different
 real(r8) :: ebv(:), sol(:), X0(:,:), y0(:), varg, vare, f(:)
 logical  :: dovce, dopcg, nonlA,nonlB
 integer  :: mx_iter, miss, nB
 optional :: dovce, mx_iter, dopcg, f, miss, nB

 integer  :: ndata, nsnp, neq, maxiter, nburning
 integer  :: i, j, io, iter, n, ioin=1
 real(r8) :: a, vara, mu, eps, val, lhs, rhs, &
             vare_ap, vara_ap, df_vare, df_vara, twosumpq
 real(r8),allocatable:: e(:), D(:), freq(:), ssol(:), xpx(:), X(:,:), y(:)
 character:: ofile*20='solutions'
 logical  :: VCE, wantPCG
 real(r8) :: curvature=1.12d0,asd(size(sol)-1),sda
 logical  :: nonlinearA=.false.,nonlinearB=.false.
 real(r8) :: solold(size(sol)),conv1,conv2
 ! nonlinearB stuff (VanRaden 2008)
 real(r8) :: q,m,ferr,fQTLerr,p1,q1,err2QTL
 character,optional:: nonl*1

 !--> inits
 neq      = size(X0,dim=2)
 nsnp     = neq - 1
 ndata    = size(pack(y0,y0.ne.miss))
 nburning = 10000
 maxiter  = 100000
 vce      = .false. !--> VarianceComponentEstimation
 wantPCG  = .true.  !--> conjugate gradient

 if(present(mx_iter)) maxiter = mx_iter
 if(present(dopcg))   wantPCG = dopcg
 if(present(dovce))   vce     = dovce
 if(present(nonl))    then
    if(nonl=='A') nonlinearA = .true.  !--> Van Raden's non linearA method
    if(nonl=='B') then                 !--> Van Raden's non linearB method
       nonlinearB = .true.
       m=nsnp; q=2000
       if(present(nB)) q=nB
       q=min(int(q),nsnp/2)
       err2QTL=m/q 
    endif
 endif
 ! Prior information for variances
 ! NOTE the 'chin' invertedchi square subroutine is parameterized in terms of
 ! degrees of freedom and expectation of the random variable and NOT
 ! in terms of d.f. and scale.
 ! That is: E(chin(S,v))=S and not to S*v
 ! Hence the way of writing the d.f. and expectations in a certain manner
 vare_ap = 9d0 !--> simulated was 10
 vara_ap = 10/(2*nsnp*0.5*0.15) !--> for 10 of genetic variance and 50000 SNPs with p=0.15
 df_vare = 4 ! low informative
 df_vara = 4

 allocate (xpx(neq),e(ndata),D(neq),freq(nsnp),ssol(neq),X(ndata,neq),y(ndata))

 
 ! Remove missing
 y = pack(y0,y0.ne.miss)
 do i=1, neq
    X(:,i) = pack(X0(:,i),y0.ne.miss)
 enddo

!--> computes freq if not provided, could be important for seq data
 if(.not.present(f)) then
    do i=1,nsnp
       freq(i)=real(sum(X(:,i+1)))/(2.*ndata)
    enddo
 else
    freq = f
 endif

!print*, count(freq==0.)+count(freq==1.), count(freq<=0.001), count(freq>=.99), size(freq), 2*sum(freq*(1-freq)), &
!       varg/(2*sum(freq*(1-freq))), vare
 

 !vare=10.0
 !vara=10/(2*nsnp*.25)      !2sumpq for 0.5

 twosumpq = 2*sum(freq*(1-freq))

 !--> var per marker
 vara = varg/twosumpq
 e    = y
 sol  = 0d0
 ssol = 0d0
 D    = 1/vara
 D(1) = 0 ! in fact D=G^{-1}
 !print *,'vara varap vare varep:', vara, vara_ap, vare, vare_ap

 !--> genotype centering
 X(:,2:neq)=X(:,2:neq)-1d0

 if (wantPCG) then
    call pcgru(X,y,D,vare,sol,maxiter)
    sda=sqrt(compute_var(sol(2:neq)))
    !print *,'crude var ahat',sda**2
    if ((sda/=0) .and. (nonlinearA .or. nonlinearB)) then
        do iter=1,100
            solold=sol
            ! iterate until convergence
            ! update following PVR 2008 nonlinearA/B
            ! sd of marker solutions
            sda=sqrt(compute_var(sol(2:neq)))
            !print *,'sd of marker solutions',sda
            do i=1,nsnp
                ! get standardised value of the marker
                asd(i)=(sol(i+1)/sda)
                !D(i+1)=1/(vara*curvature**abs(asd-2))
                !if(mod(i,100)==0) print *,curvature**(asd-2)
            enddo
            !print *,'meanasd=',mean(asd),'max=',maxval(asd),'min=',minval(asd)
            if(nonlinearA) then
                    asd=abs(asd)
                    where(asd>5d0)
                        asd=5d0
                    endwhere
                    !print *,'meanasd=',mean(asd),'max=',maxval(asd),'min=',minval(asd)
                    D(2:neq)=1d0/(vara*curvature**(asd-2))
            elseif(nonlinearB) then
                    if(iter==1) print *,'nonlinearB, number of QTL: ',int(q),'number of markers', nsnp
                    do i=1,nsnp
                        asd(i)=abs(asd(i))
                        if(asd(i)>15d0) asd(i)=15d0
                        ! density under the null (noQTL)
                        ferr=dnormal(asd(i),0d0,1d0)
                        ! density under the alternative (QTL)
                        fQTLerr=dnormal(asd(i),0d0,err2QTL)
                        D(1+i)=(1d0/(vara))*(q/m+(1-q/m)*(ferr/fQTLerr))
                    enddo
            endif
            !    asd=minval((/asd,5d0/))
            !print *,'vara0=',vara,'maxvara=',maxval(1/D(2:neq)),'minvara=',minval(1/D(2:neq))
            !print *,'changes in shrinkage by nonlinear:',nonl,'nsnp',int(nsnp),'nqtl assumed',int(q)
            !print *,'relatives',1,'maxvara=',maxval(1/D(2:neq))/vara,'minvara=',minval(1/D(2:neq))/vara
            call pcgru(X,y,D,vare,sol,maxiter)
            ! convergence on all solutions
            conv1=sum( (sol-solold)**2/sum(sol**2) )
            ! convergence on markers
            conv2=sum( (sol(2:neq)-solold(2:neq))**2/sum(sol(2:neq)**2) )
            !print *,iter,'convv=',conv1,conv2
            if(conv1<1e-8) then
                print *,'converged',iter,conv1
                print *,'relative changes in shrinkage by nonlinear:',nonl,'nsnp',int(nsnp),'nqtl assumed',int(q)
                print *,'maxvara/vara=',maxval(1/D(2:neq))/vara,'minvara/vara=',minval(1/D(2:neq))/vara
                exit
            endif
            if (iter==100) then
                print *,'nonlinear failed convergence'
                !stop
            endif
        enddo
    endif
 else
    do i=1,neq
       xpx(i)=dot_product(X(:,i),X(:,i)) !form diagonal of X′X
    enddo
    !--> Gibbs sampling to get var estimates, slooow
    do iter=1,maxiter
       eps=0
       do i=1,neq
          !form lhs
          lhs=xpx(i)/vare+D(i)
          ! form rhs with y corrected by other effects (formula 1)
          rhs=dot_product(X(:,i),e)/vare +xpx(i)*sol(i)/vare
          ! do Gauss Seidel
          val=rhs/lhs
          ! MCMC sample solution from its conditional (commented out here)
          if(VCE) val=normal(val,1d0/lhs)
          eps=eps+(val-sol(i))**2
          ! update e with current estimate (formula 2)
          e = e - X(:,i)*(val-sol(i))
          ! update sol
          sol(i)=val
       enddo
       if (.not.VCE)then
          eps=eps/sum(sol**2)
          print*,iter,'eps: ',eps
          if(eps<10d-12) EXIT
       else !--> VCE
          vara = (dot_product(sol(2:neq),sol(2:neq))+ df_vara*vara_ap) / (nsnp+df_vara)
          vara = chin(vara, nsnp+df_vara)
          vare = (dot_product(e,e) + df_vare*vare_ap) / (ndata+df_vare)
          vare = chin(vare, ndata+df_vare)
          print*,'vara,vare',vara,vare
       endif
       if(iter>nburning) ssol=ssol+sol
    enddo !--> end iter
    if(VCE) sol=ssol/(maxiter-nburning)
 endif

 !--> X0 contains genotypes of missing y individuals as well
 ebv(:) = matmul(X0,sol)

 deallocate(xpx,e,D,freq,ssol,X,y)

!--------------
 end subroutine
!--------------


!-----------------------
 subroutine solveSNPsub1 (ebv, sol, varg, vare, X0, y0, f, dovce, dopcg, mx_iter, miss, nonl)
!-----------------------
! SNP-RR by Andres Legarra (Andres.Legarra@toulouse.inra.fr)
! Reshaped to subroutine by M. Perez-Enciso (miguel.perez@uab.es)
! Allows for mssing values by removing rows in X and y (mpe)
! incorporates Vanraden's approx non linear approach
 real(r8) :: ebv(:), sol(:), X0(:,:), y0(:), varg, vare, f(:)
 logical  :: dovce, dopcg, nonl
 integer  :: mx_iter, miss
 optional :: dovce, mx_iter, dopcg, f, miss, nonl

 integer  :: ndata, nsnp, neq, maxiter, nburning
 integer  :: i, j, io, iter, n, ioin=1
 real(r8) :: a, vara, mu, eps, val, lhs, rhs, &
             vare_ap, vara_ap, df_vare, df_vara, twosumpq
 real(r8),allocatable:: e(:), D(:), freq(:), ssol(:), xpx(:), X(:,:), y(:)
 character:: ofile*20='solutions'
 logical  :: VCE, wantPCG
 real(r8) :: curvature=1.12d0,asd(size(sol)-1),sda
 logical  :: nonlinearA=.false.
 real(r8) :: solold(size(sol)),conv1,conv2

 !--> inits
 neq      = size(X0,dim=2)
 nsnp     = neq - 1
 ndata    = size(pack(y0,y0.ne.miss))
 nburning = 10000
 maxiter  = 100000
 vce      = .false. !--> VarianceComponentEstimation
 wantPCG  = .true.  !--> conjugate gradient
 if(present(mx_iter)) maxiter = mx_iter
 if(present(dopcg))   wantPCG = dopcg
 if(present(dovce))   vce     = dovce
 if(present(nonl))    nonlinearA = nonl  !--> Van Raden's non linear method

 ! Prior information for variances
 ! NOTE the 'chin' invertedchi square subroutine is parameterized in terms of
 ! degrees of freedom and expectation of the random variable and NOT
 ! in terms of d.f. and scale.
 ! That is: E(chin(S,v))=S and not to S*v
 ! Hence the way of writing the d.f. and expectations in a certain manner
 vare_ap = 9d0 !--> simulated was 10
 vara_ap = 10/(2*nsnp*0.5*0.15) !--> for 10 of genetic variance and 50000 SNPs with p=0.15
 df_vare = 4 ! low informative
 df_vara = 4

 allocate (xpx(neq),e(ndata),D(neq),freq(nsnp),ssol(neq),X(ndata,neq),y(ndata))

 ! Remove missing
 y = pack(y0,y0.ne.miss)
 do i=1, neq
    X(:,i) = pack(X0(:,i),y0.ne.miss)
 enddo

 !vare=10.0
 !vara=10/(2*nsnp*.25)      !2sumpq for 0.5

 !--> computes freq if not provided, could be important for seq data
 if(.not.present(f)) then
    do i=1,nsnp
       freq(i)=real(sum(X(:,i+1)))/(2.*ndata)
    enddo
 else
    freq = f
 endif
 twosumpq = 2*sum(freq*(1-freq))
 !--> var per marker
 vara = varg/twosumpq
 e    = y
 sol  = 0d0
 ssol = 0d0
 D    = 1/vara
 D(1) = 0 ! in fact D=G^{-1}
 !print *,'vara varap vare varep:', vara, vara_ap, vare, vare_ap

 !--> genotype centering
 X(:,2:neq)=X(:,2:neq)-1d0

 if (wantPCG) then
    call pcgru(X,y,D,vare,sol,maxiter)
    sda=sqrt(compute_var(sol(2:neq)))
    print *,'crude var ahat',sda**2
    if ((sda/=0) .and. nonlinearA) then
        do iter=1,100
            solold=sol
            ! iterate until convergence
            ! update following PVR 2008 nonlinearA
            ! sd of marker solutions
            sda=sqrt(compute_var(sol(2:neq)))
            !print *,'sd of marker solutions',sda
            do i=1,nsnp
                ! get standardised value of the marker
                asd(i)=(sol(i+1)/sda)
                !D(i+1)=1/(vara*curvature**abs(asd-2))
                !if(mod(i,100)==0) print *,curvature**(asd-2)
            enddo
            !print *,'meanasd=',mean(asd),'max=',maxval(asd),'min=',minval(asd)
            asd=abs(asd)
            where(asd>5d0)
                asd=5d0
            endwhere
            !print *,'meanasd=',mean(asd),'max=',maxval(asd),'min=',minval(asd)
            D(2:neq)=1d0/(vara*curvature**(asd-2))

            !    asd=minval((/asd,5d0/))
            !print *,'changes in shrinkage by nonlinearA'
            !print *,'vara0=',vara,'maxvara=',maxval(1/D(2:neq)),'minvara=',minval(1/D(2:neq))
            !print *,'relatives',1,'maxvara=',maxval(1/D(2:neq))/vara,'minvara=',minval(1/D(2:neq))/vara
            call pcgru(X,y,D,vare,sol,maxiter)
            ! convergence on all solutions
            conv1=sum( (sol-solold)**2/sum(sol**2) )
            ! convergence on markers
            conv2=sum( (sol(2:neq)-solold(2:neq))**2/sum(sol(2:neq)**2) )
            !print *,iter,'convv=',conv1,conv2
            if(conv1<1e-8) then
                print *,'converged',iter,conv1
                exit
            endif
            if (iter==100) then
                print *,'nonlinearA failed convergence'
                !stop
            endif
        enddo
    endif
 else
    do i=1,neq
       xpx(i)=dot_product(X(:,i),X(:,i)) !form diagonal of X′X
    enddo
    !--> Gibbs sampling to get var estimates, slooow
    do iter=1,maxiter
       eps=0
       do i=1,neq
          !form lhs
          lhs=xpx(i)/vare+D(i)
          ! form rhs with y corrected by other effects (formula 1)
          rhs=dot_product(X(:,i),e)/vare +xpx(i)*sol(i)/vare
          ! do Gauss Seidel
          val=rhs/lhs
          ! MCMC sample solution from its conditional (commented out here)
          if(VCE) val=normal(val,1d0/lhs)
          eps=eps+(val-sol(i))**2
          ! update e with current estimate (formula 2)
          e=e - X(:,i)*(val-sol(i))
          ! update sol
          sol(i)=val
       enddo
       if (.not.VCE)then
          eps=eps/sum(sol**2)
          print *,iter,'eps: ',eps
          if(eps<10d-12) EXIT
       else !--> VCE
          vara = (dot_product(sol(2:neq),sol(2:neq))+ df_vara*vara_ap) / (nsnp+df_vara)
          vara = chin(vara, nsnp+df_vara)
          vare = (dot_product(e,e) + df_vare*vare_ap) / (ndata+df_vare)
          vare = chin(vare, ndata+df_vare)
          print*,'vara,vare',vara,vare
       endif
       if(iter>nburning) ssol=ssol+sol
    enddo !--> end iter
    if(VCE) sol=ssol/(maxiter-nburning)
 endif

 ebv(:) = matmul(X0,sol)

 deallocate(xpx,e,D,freq,ssol,X,y)

!--------------
 end subroutine
!--------------
!----------
 END MODULE
!----------
