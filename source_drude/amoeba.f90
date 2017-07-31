module amoeba_routines
implicit none

contains

	SUBROUTINE amoeba(p,y,ftol,iter)
	IMPLICIT NONE
	INTEGER, INTENT(OUT) :: iter
	REAL*8, INTENT(IN) :: ftol
	REAL*8, DIMENSION(:), INTENT(INOUT) :: y
	REAL*8, DIMENSION(:,:), INTENT(INOUT) :: p
	INTEGER, PARAMETER :: ITMAX=200
	REAL*8, PARAMETER :: TINY=1.0e-10
	INTEGER :: ihi,ndim
	REAL*8, DIMENSION(size(p,2)) :: psum
	call amoeba_private
	CONTAINS
!BL
	SUBROUTINE amoeba_private
	IMPLICIT NONE
	INTEGER :: i,ilo,inhi
	REAL*8 :: rtol,ysave,ytry,ytmp
	ndim=size(p,2)
        if ( (size(p,2) .ne. size(p,1)-1 ) .or. (size(p,2) .ne. size(y)-1) ) then
           stop "error in input array sizes in amoeba"
        endif


iter=0
	psum(:)=sum(p(:,:),dim=1)
	do
       
           ilo=iminloc(y(:))
		ihi=imaxloc(y(:))
		ytmp=y(ihi)
		y(ihi)=y(ilo)
		inhi=imaxloc(y(:))
		y(ihi)=ytmp
		rtol=2.0d0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
		if (rtol < ftol) then
			call swap(y(1),y(ilo))
			call swap2(p(1,:),p(ilo,:))
			RETURN
		end if
		if (iter >= ITMAX) then
                   stop "ITMAX exceeded in amoeba"
                endif
		ytry=amotry(-1.0d0)
		iter=iter+1
		if (ytry <= y(ilo)) then
			ytry=amotry(2.0d0)
			iter=iter+1
		else if (ytry >= y(inhi)) then
			ysave=y(ihi)
			ytry=amotry(0.5d0)
			iter=iter+1
			if (ytry >= ysave) then
				p(:,:)=0.5d0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
				do i=1,ndim+1
					if (i /= ilo) y(i)=func(p(i,:))
				end do
				iter=iter+ndim
				psum(:)=sum(p(:,:),dim=1)
			end if
		end if
	end do
	END SUBROUTINE amoeba_private
!BL
	FUNCTION amotry(fac)
	IMPLICIT NONE
	REAL*8, INTENT(IN) :: fac
	REAL*8 :: amotry
	REAL*8 :: fac1,fac2,ytry
	REAL*8, DIMENSION(size(p,2)) :: ptry
	fac1=(1.0d0-fac)/ndim
	fac2=fac1-fac
	ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
	ytry=func(ptry)
        	if (ytry < y(ihi)) then
		y(ihi)=ytry
		psum(:)=psum(:)-p(ihi,:)+ptry(:)
		p(ihi,:)=ptry(:)
	end if
	amotry=ytry
	END FUNCTION amotry
	END SUBROUTINE amoeba

	FUNCTION iminloc(arr)
	REAL*8, DIMENSION(:), INTENT(IN) :: arr
	INTEGER, DIMENSION(1) :: imin
	INTEGER :: iminloc
	imin=minloc(arr(:))
	iminloc=imin(1)
	END FUNCTION iminloc

	FUNCTION imaxloc(arr)
	REAL*8, DIMENSION(:), INTENT(IN) :: arr
	INTEGER :: imaxloc
	INTEGER, DIMENSION(1) :: imax
	imax=maxloc(arr(:))
	imaxloc=imax(1)
	END FUNCTION imaxloc

	SUBROUTINE swap(a,b)
	REAL*8, INTENT(INOUT) :: a,b
	REAL*8 :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap

	SUBROUTINE swap2(a,b)
	REAL*8,dimension(:), INTENT(INOUT) :: a,b
	REAL*8,dimension(size(a)) :: dum
	dum=a
	a=b
	b=dum
	END SUBROUTINE swap2


real*8 function func(param)
  use variables
  use routines
  use eq_drude
  real*8,dimension(:),intent(in) :: param

  real*8 :: rms,max_d,min_d,temp_charge
  integer:: i,j,atms,atoms

  atms=size(atmtype)
  atoms=size(atomtype)

  call fill_charges(param)

  write(*,*) shellcharge(:)
  call equilibrate_drude_positions
  call rms_p2p(rms,max_d,min_d)
  write(*,*) rms

 ! add penalty if charges are less than a certain values
  do i=1,atms
     do j=1,atoms
        if (atmtype(i) .eq. atomtype(j) ) then
           temp_charge=shellcharge(j)
        endif
     enddo
     if (abs (temp_charge) < shell_min ) then
        rms = rms * shell_min/( abs (temp_charge) ) 
     endif
  enddo
           

  func= rms

end function func

end module amoeba_routines
