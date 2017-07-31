module electrostatic_interaction_tensor

contains

!********************************************************
! this function evaluates the electrostatic interaction functions
! Ttu that arise in a spherical-tensor multipole expansion of the 
! coulomb operator. "t" represents the rank of site a, u represents
! the rank of site b.  these have values 00,20,21c,21s, etc
! from The Theory of Intermolecular Forces, AJ Stone
!
! this function returns tensor components Tab in the 1st element of the 
! 3rd index, and tensor components Tba in the 2nd element of the 3rd index
! the first two indices run over all components of the given rank for the
! first site, with the second site defined as 0,0,0
! so for Tab these are la,ma,csa,0,0,0
! and for Tba these are lb,mb,csb,0,0,0
! where la, lb are the rank of sites a and b
! ma, csa, mb,csb determine tensor component
! cs=0 means cosine, cs=1 means sine
!
! for example (la,ma,csa)=(2,1,0) = 21c
! (la,ma,csa)=(2,2,1) = 22s
! 
! in practice, site A is an atom, and site B is a grid point
! remember, distance is always defined from a to b, so 
! if we are generating Tba tensors, need to change sign on direction
! cosine
!*********************************************************

function elec_int_tensor(Rab,ra_cos,rb_cos,rank,inverse)
  real*8,dimension(3),intent(in)::Rab,ra_cos,rb_cos
  integer,intent(in) :: rank,inverse
  real*8,dimension(3,0:3,0:1,2) :: elec_int_tensor

  real*8 :: R_inv2,R_inv3,R_inv4

! get inverse distances
  R_inv2=Rab(1);  R_inv3=Rab(2);  R_inv4=Rab(3);

! unit vector is from site a to b, so make sure signs are correct


! don't need c matrix for these tensors
!!$  do j=1,3
!!$     c_a_b(i,j) = dot_product(locA(i,:),locB(j,:))
!!$  enddo

! if statements will slow the code down, so evaluate all components up to rank 3

 elec_int_tensor=0d0

! if inverse is 0, calculate Tab, otherwise if inverse is 1, calculate Tba
Select Case(inverse)
   Case(0)

!**************** 1xx000 tensors
  elec_int_tensor(1,0,0,1)= t_1_0_0_0_0_0(R_inv2, ra_cos)
  elec_int_tensor(1,1,0,1)= t_1_1_0_0_0_0(R_inv2, ra_cos)
  elec_int_tensor(1,1,1,1)= t_1_1_1_0_0_0(R_inv2, ra_cos)

!****************2xx000 tensors
  elec_int_tensor(2,0,0,1)= t_2_0_0_0_0_0(R_inv3, ra_cos)
  elec_int_tensor(2,1,0,1)= t_2_1_0_0_0_0(R_inv3, ra_cos)
  elec_int_tensor(2,1,1,1)= t_2_1_1_0_0_0(R_inv3, ra_cos)
  elec_int_tensor(2,2,0,1)= t_2_2_0_0_0_0(R_inv3, ra_cos)
  elec_int_tensor(2,2,1,1)= t_2_2_1_0_0_0(R_inv3, ra_cos)

!****************3xx000 tensors
  elec_int_tensor(3,0,0,1)= t_3_0_0_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,1,0,1)= t_3_1_0_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,1,1,1)= t_3_1_1_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,2,0,1)= t_3_2_0_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,2,1,1)= t_3_2_1_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,3,0,1)= t_3_3_0_0_0_0(R_inv4, ra_cos)
  elec_int_tensor(3,3,1,1)= t_3_3_1_0_0_0(R_inv4, ra_cos)

  Case(1)

!*****************0001xx tensors
  elec_int_tensor(1,0,0,2)= t_1_0_0_0_0_0(R_inv2, -rb_cos)
  elec_int_tensor(1,1,0,2)= t_1_1_0_0_0_0(R_inv2, -rb_cos)
  elec_int_tensor(1,1,1,2)= t_1_1_1_0_0_0(R_inv2, -rb_cos)

!****************0002xx tensors
  elec_int_tensor(2,0,0,2)= t_2_0_0_0_0_0(R_inv3, -rb_cos)
  elec_int_tensor(2,1,0,2)= t_2_1_0_0_0_0(R_inv3, -rb_cos)
  elec_int_tensor(2,1,1,2)= t_2_1_1_0_0_0(R_inv3, -rb_cos)
  elec_int_tensor(2,2,0,2)= t_2_2_0_0_0_0(R_inv3, -rb_cos)
  elec_int_tensor(2,2,1,2)= t_2_2_1_0_0_0(R_inv3, -rb_cos)

!****************0003xx tensors
  elec_int_tensor(3,0,0,2)= t_3_0_0_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,1,0,2)= t_3_1_0_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,1,1,2)= t_3_1_1_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,2,0,2)= t_3_2_0_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,2,1,2)= t_3_2_1_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,3,0,2)= t_3_3_0_0_0_0(R_inv4, -rb_cos)
  elec_int_tensor(3,3,1,2)= t_3_3_1_0_0_0(R_inv4, -rb_cos)

End Select


end function elec_int_tensor



!**********************************HERE ARE THE TENSOR FUNCTIONS**********************************
!*************************************************************************************************


!***************************** rank 1xx000 tensors*************************************

function t_1_0_0_0_0_0(R_inv2, r_cos )
  real*8,intent(in) :: R_inv2
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_1_0_0_0_0_0
t_1_0_0_0_0_0 = R_inv2 * r_cos(3)
end function t_1_0_0_0_0_0

function t_1_1_0_0_0_0(R_inv2, r_cos )
  real*8,intent(in) :: R_inv2
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_1_1_0_0_0_0
t_1_1_0_0_0_0 = R_inv2 * r_cos(1)
end function t_1_1_0_0_0_0

function t_1_1_1_0_0_0(R_inv2, r_cos )
  real*8,intent(in) :: R_inv2
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_1_1_1_0_0_0
t_1_1_1_0_0_0 = R_inv2 * r_cos(2)
end function t_1_1_1_0_0_0

!**************************** rank 2xx000 tensors****************************************

function t_2_0_0_0_0_0(R_inv3, r_cos )
  real*8,intent(in) :: R_inv3
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_2_0_0_0_0_0
t_2_0_0_0_0_0 = R_inv3 * .5d0*(3d0*r_cos(3)**2-1d0)
end function t_2_0_0_0_0_0

function t_2_1_0_0_0_0(R_inv3, r_cos)
  real*8,intent(in) :: R_inv3
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_2_1_0_0_0_0
t_2_1_0_0_0_0 = R_inv3 * sqrt(3d0)*r_cos(1)*r_cos(3)
end function t_2_1_0_0_0_0

function t_2_1_1_0_0_0(R_inv3, r_cos)
  real*8,intent(in) :: R_inv3
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_2_1_1_0_0_0
t_2_1_1_0_0_0 = R_inv3 * sqrt(3d0)*r_cos(2)*r_cos(3)
end function t_2_1_1_0_0_0

function t_2_2_0_0_0_0(R_inv3, r_cos)
  real*8,intent(in) :: R_inv3
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_2_2_0_0_0_0
t_2_2_0_0_0_0 = R_inv3 * .5d0*sqrt(3d0)*(r_cos(1)**2 - r_cos(2)**2)
end function t_2_2_0_0_0_0

function t_2_2_1_0_0_0(R_inv3, r_cos)
  real*8,intent(in) :: R_inv3
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_2_2_1_0_0_0
t_2_2_1_0_0_0 = R_inv3 * sqrt(3d0)*r_cos(1)*r_cos(2)
end function t_2_2_1_0_0_0

!*********************** rank 3xx000 tensors**********************************************

function t_3_0_0_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_0_0_0_0_0
t_3_0_0_0_0_0 = R_inv4 * .5d0*(5d0*r_cos(3)**3-3d0*r_cos(3))
end function t_3_0_0_0_0_0

function t_3_1_0_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_1_0_0_0_0
t_3_1_0_0_0_0 = R_inv4 * .25d0*sqrt(6d0)*r_cos(1)*(5d0*r_cos(3)**2-1d0)
end function t_3_1_0_0_0_0

function t_3_1_1_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_1_1_0_0_0
t_3_1_1_0_0_0 = R_inv4 * .25d0*sqrt(6d0)*r_cos(2)*(5d0*r_cos(3)**2-1d0)
end function t_3_1_1_0_0_0

function t_3_2_0_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_2_0_0_0_0
t_3_2_0_0_0_0 = R_inv4 * .5d0*sqrt(15d0)*r_cos(3)*(r_cos(1)**2-r_cos(2)**2)
end function t_3_2_0_0_0_0

function t_3_2_1_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_2_1_0_0_0
t_3_2_1_0_0_0 = R_inv4 * sqrt(15d0)*r_cos(3)*r_cos(1)*r_cos(2)
end function t_3_2_1_0_0_0

function t_3_3_0_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_3_0_0_0_0
t_3_3_0_0_0_0 = R_inv4 * .25d0*sqrt(10d0)*r_cos(1)*(r_cos(1)**2-3d0*r_cos(2)**2)
end function t_3_3_0_0_0_0

function t_3_3_1_0_0_0(R_inv4, r_cos )
  real*8,intent(in) :: R_inv4
  real*8,dimension(3),intent(in) ::r_cos
  real*8  :: t_3_3_1_0_0_0
t_3_3_1_0_0_0 = R_inv4 * .25d0*sqrt(10d0)*r_cos(2)*(3d0*r_cos(1)**2-r_cos(2)**2)
end function t_3_3_1_0_0_0

end module electrostatic_interaction_tensor
