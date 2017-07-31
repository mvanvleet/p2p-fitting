module eq_drude

contains

  subroutine rms_p2p(rms,max_d,min_d)
    use variables
    real*8,intent(out)::rms,max_d,min_d

    real*8 ::sum1,potential
    integer:: i,j,n_points,count

    rms=0d0
    n_points = size(xyz_pts(:,1))
    max_d=0d0
    min_d=0d0
    sum1=0d0
  count=0

    ! loop over all point to point responses
    do i=1,n_points
       do j=1,i

          potential = electrostatic_potential(i,j)
          ! max and min
          if (  ( potential - p2p_response(i_freq,i,j) ) > max_d ) then
             max_d = potential - p2p_response(i_freq,i,j)
          elseif (  (potential - p2p_response(i_freq,i,j) ) < min_d ) then
             min_d = potential - p2p_response(i_freq,i,j)
          endif
          
       sum1= sum1 + ( potential - p2p_response(i_freq,i,j))**2
       count=count+1
    enddo
 enddo

 rms = sqrt(sum1/dble(count))

end subroutine rms_p2p


subroutine initialize_drude_oscillators
use variables
integer :: i,atoms,points

atoms=size(xyz_molecule(:,1))
points=size(xyz_pts(:,1))

allocate(xyz_drude_eq(points,atoms,3),shellcharge(atoms) )

! static dipole-dipole polarizabilitiess stored in dist_pol_tensor(:,n_freq+1,1)
do i=1,atoms
   shellcharge(i) = -sqrt( dist_pol_tensor(i,n_freq+1,1) * springconstant )
enddo

do i=1,points
   ! put oscillators on atoms as intitial guess
   do j=1,atoms
      xyz_drude_eq(i,j,:)=xyz_molecule(j,:)
   enddo
enddo

call equilibrate_drude_positions

end subroutine initialize_drude_oscillators


subroutine equilibrate_drude_positions
use variables

integer :: i,j,points,atoms

atoms=size(xyz_molecule(:,1))
points=size(xyz_pts(:,1))

! get positions of drude oscillators for these polarizabilities

do i=1,points
   call findpositionshell(i)
enddo

end subroutine equilibrate_drude_positions


!*****************************
! this subroutine evaluates the electrostatic potential at grid point P
! due to the molecular linear response to a point charge at point Q
!*****************************
real*8 function  electrostatic_potential(P,Q)
use variables
integer,intent(in) :: P,Q

integer :: i,atoms
real*8,dimension(3) :: xyz
real*8 :: dist,sum

  atoms=size(xyz_molecule(:,1))
  sum=0d0

  do i=1,atoms
     ! shells
     xyz(:)= xyz_pts(P,:)-xyz_drude_eq(Q,i,:)
     dist=sqrt(dot_product(xyz,xyz))
     sum = sum + shellcharge(i)/dist
     ! atom centers
     xyz(:)= xyz_pts(P,:)-xyz_molecule(i,:)     
     dist=sqrt(dot_product(xyz,xyz))
     sum = sum - shellcharge(i)/dist
  enddo


electrostatic_potential=sum

end function electrostatic_potential


 recursive subroutine findpositionshell(n)
    use variables
    integer,intent(in)::n
    real*8,dimension(:,:),allocatable::forces
    integer::atoms,i,j,k,l,iter,converged,anharm
    real*8::delta,r,sum
    real*8,dimension(3)::xxout,del,dxyz
    real*8,parameter::step=0.1,step1=.005,small=1D-7


    atoms=size(xyz_molecule(:,1))

    allocate(forces(atoms,3))
    converged=1

    do i=1,atoms
       if ( abs(shellcharge(i)) > 0d0 ) then
       dxyz=xyz_drude_eq(n,i,:)-xyz_molecule(i,:)
       r=sqrt(dot_product(dxyz,dxyz))
       forces(i,:)=shellcharge(i)*field(i,n)


       do j=1,3
          delta=abs(forces(i,j)-springconstant*dxyz(j))
          if(delta > deltashell) converged=0
       enddo
       endif
    enddo
    if (converged.eq.0) then
       do i=1,atoms
       if ( abs(shellcharge(i)) > 0d0 ) then
          xxout(:)=forces(i,:)/springconstant
          del(:)=xyz_molecule(i,:)+xxout(:)-xyz_drude_eq(n,i,:)

          xyz_drude_eq(n,i,:)=xyz_drude_eq(n,i,:)+step*del(:)

          else
             ! if zero charge, put on atom
           xyz_drude_eq(n,i,:)= xyz_molecule(i,:)
        endif
       enddo

       deallocate(forces)
       call findpositionshell(n)

     endif

  end subroutine findpositionshell


function field(shell,n)
  use variables

  real*8,dimension(3)::field
  integer,intent(in)::shell,n
  real*8::dist
  real*8,dimension(3)::xyz,sum_field
  integer::i,atoms
  real*8 :: sign

  ! this determines whether response is to positive or negative point charge
  sign=-1d0

  atoms=size(xyz_molecule(:,1))


!!!!!!!!!!!!!!!field at shell due to point charge on grid

     xyz(:)= xyz_drude_eq(n,shell,:)-xyz_pts(n,:)
     dist=sqrt(dot_product(xyz,xyz))
     sum_field= sign * xyz(:)/dist**3

!!!!!!!!!!!!! field due to intra molecular contributions
  do i=1,atoms
     if ( abs(shellcharge(i)) > 0d0 ) then
!!!!!!!!!!!!! field due to other shells
     if(i .ne.shell) then
        xyz(:)=xyz_drude_eq(n,shell,:)-xyz_drude_eq(n,i,:)
        dist=sqrt(dot_product(xyz,xyz))
        sum_field(:)=sum_field(:)+screen(shell,i,dist)*shellcharge(i)*xyz(:)/dist**3
        sum_field(:)=sum_field(:)+(-1)*dscreen(shell,i,xyz)*shellcharge(i)/dist


!!!!!!!!!!!!! field from shell charge on atom center
        xyz(:)=xyz_drude_eq(n,shell,:)-xyz_molecule(i,:)
        dist=sqrt(dot_product(xyz,xyz))
        sum_field(:)=sum_field(:)+screen(shell,i,dist)*(-shellcharge(i))*xyz(:)/dist**3
        sum_field(:)=sum_field(:)+(-1)*dscreen(shell,i,xyz)*(-shellcharge(i))/dist
     endif
     endif
  enddo

  field=sum_field

end function field



function screen(atom1,atom2,dist)
  use variables
  real*8::screen
  integer,intent(in)::atom1,atom2
  real*8,intent(in)::dist
  integer::i,j
  real*8::pol1,pol2,a,u
  real*8::au3, expau3, approx_gamma
  real*8,parameter::small=1D-6
  
  i= atom1
  j= atom2

  ! set screen = zero for zero charge otherwise blows up

     if((abs(shellcharge(i)).lt. small).or.( abs(shellcharge(j)).lt. small)) then
        screen=0d0
     else
      pol1=(shellcharge(i)**2)/springconstant
      pol2=(shellcharge(j)**2)/springconstant
      !! a=thole/(pol1*pol2)**(1.0/6.0)
      !! screen=1.0-(1.0+(a*dist)/(2.*(pol1*pol2)**(1./6.)))*exp(-a*dist/(pol1*pol2)**(1./6.))

      a=thole
      u=dist/((pol1*pol2)**(1.0/6.0))

      au3 = a*u**(3.0)
      expau3 = exp(-au3)
      approx_gamma = expau3*((4.0/9.0) + 2*(1+au3)**2) / (2*(1+ au3)**(7.0/3.0))
      screen = 1.0-expau3 + u*a**(1.0/3.0)*approx_gamma

                  
                  !! *dist*(4/9.0 + 2*(1 + a*dist**3)**2))/(6.*exp(a*dist**3)*(1 + a*dist**3)**(7.0/3.0))

      !! write(*,*) 'screen = ', screen
      endif

end function screen

function  dscreen(atom1,atom2,xyz)
  use variables
  real*8,dimension(3)::dscreen
  integer,intent(in)::atom1,atom2
  !! real*8,dimension(3),intent(in)::xyz
  real*8,dimension(3) :: xyz
  integer::i,j
  real*8::pol1,pol2,r,fac,Ex,a, fac1
  real*8:: au3, expau3
  real*8,parameter::small=1D-6

  !! a=thole
  i= atom1
  j= atom2

  r=(xyz(1)**2+xyz(2)**2+xyz(3)**2)**.5

     ! set dscreen = zero for zero charge otherwise blows up
     if((abs(shellcharge(i)).lt. small).or.( abs(shellcharge(j)).lt. small)) then
        dscreen =0d0
     else
        pol1=(shellcharge(i)**2)/springconstant
        pol2=(shellcharge(j)**2)/springconstant
        ! a=thole/(pol1*pol2)**(1.0/6.0)

        a=thole
        u=r/(pol1*pol2)**(1.0/6.0)
        !! pol=1.0/(pol1*pol2)**(1.0/6.0)
        !! Ex=exp(-fac*dist)
        !! dscreen(:)=(xyz(:)/dist)*(fac*(1.0+fac*dist/2.)*Ex-(fac/2.)*Ex)
        !! dscreen(:) = (11 + 27*a**4*xyz(:)*(r**2)**5*(-xyz(:) + (1 + a*(r**2)**1.5)**(1/3.0)) +&
        !!                 3*a**2*(r**2)**2*(-26*xyz(:)**2 + 9*(r**2 - xyz(:)**2) +&
        !!                 27*xyz(:)*(1 + a*(r**2)**1.5)**(1/3.0)) +&
        !!                 9*a**3*(r**2)**3.5*(r**2 - xyz(:)**2 + 9*xyz(:)*&
        !!                 (-xyz(:) + (1 + a*(r**2)**1.5)**(1/3.0)))&
        !!                 + a*r*(29*(r**2 - xyz(:)**2) +& 
        !!                 27*xyz(:)*(-xyz(:) + (1 + a*(r**2)**1.5)**(1/3.0))))/&
        !!                 (9.*exp(a*(r**2)**1.5)* (1 + a*(r**2)**1.5)**(10/3.0))

        !! a = 1
        !! xyz(:) = (/ 2,0,0 /)
        !! r = 2
        !! pol = 1
        !! u=r*pol

        au3 = a*u**3
        expau3 = exp(-au3)
        !! approx_gamma = expau3*((4.0/9.0) + 2*(1+au3)**2) / (2*(1+ au3)**(7.0/3.0))

        !! dscreen(:) = xyz(:)*(a**(1.0/3.0) * expau3 &
        !!                 * ( -23*au3 - 18*(au3)**2 - 9*(au3)**3 + 11*(1 +au3)&
        !!                     - 15*au3*(1+au3) - (45*(au3)**2)*(1 + au3) &
        !!                     - (27 *(au3)**3)*(1+au3) &
        !!                     + 27*a**(2.0/3.0)*u**2*(1+au3)**(10.0/3.0))) &
        !!                  / (9.0*pol*r*(1+au3)**(10.0/3.0) )

        fac = 11 - 27*au3 - 78*au3**2 - 81*au3**3 - 27*au3**4 &
                               + 27*a**(2.0/3.0)*u**2*(1+ au3)**(10.0/3.0)

        !! write(*,*)  'au3', au3
        !! write(*,*)  'expau3', expau3
        !! write(*,*)  'fac1', - 27*au3
        !! write(*,*)  'fac2', - 78*au3**2
        !! write(*,*)  'fac3', - 81*au3**3
        !! write(*,*)  'fac4', - 27*au3**4
        !! write(*,*)  'fac5', + 27*a**(2.0/3.0)*u**2*(1+ au3)**(10.0/3.0)
        !! write(*,*) 'fac', fac
        dscreen(:) = xyz(:)*(a**(1.0/3.0) * expau3 * fac ) &
                         / (9.0*u*(1+au3)**(10.0/3.0) )

        fac1=1d0/(pol1*pol2)**(1./6.)
        dscreen(:) = dscreen(:) * fac1**3

        !! write(*,*) dscreen(:)
        !! stop
     endif

end function dscreen


end module eq_drude
