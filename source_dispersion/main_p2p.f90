!*****************************************************************
! this program is for fitting polarizability tensors (freq dependent and static) to the
! point to point response calculations performed by camcasp
! this only fits isotropic diagonal tensor components, as this is needed for 
! a purely isotropic C6,C8,C10 treatment
!*****************************************************************

program p2p
  use variables
  use routines
  interface
     subroutine construct_G_matrix(G_matrix)
       use variables
       real*8,dimension(:,:), intent(out) :: G_matrix
     end subroutine construct_G_matrix
     subroutine construct_f_vector(f_vector)
       use variables
       real*8,dimension(:),intent(out) :: f_vector
     end subroutine construct_f_vector
     SUBROUTINE svdcmp_sp(a,w,v)
       USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
       USE nr, ONLY : pythag
       IMPLICIT NONE
       REAL*8, DIMENSION(:,:), INTENT(INOUT) :: a
       REAL*8, DIMENSION(:), INTENT(OUT) :: w
       REAL*8, DIMENSION(:,:), INTENT(OUT) :: v
     END SUBROUTINE svdcmp_sp
     SUBROUTINE svbksb_sp(u,w,v,b,x)
       USE nrtype; USE nrutil, ONLY : assert_eq
       REAL*8, DIMENSION(:,:), INTENT(IN) :: u,v
       REAL*8, DIMENSION(:), INTENT(IN) :: w,b
       REAL*8, DIMENSION(:), INTENT(OUT) :: x
     END SUBROUTINE svbksb_sp
  end interface

  character(30) :: base
  character(80) :: path,consfile
  character(8)::date
  character(10)::time
  integer:: i,j,k,n_param,atoms,atms,h_atms
  real*8,dimension(:),allocatable::p0,p1,df0,delta,f_vector,w_svd
  real*8,dimension(:,:),allocatable::G_matrix,G_matrix_store,v_svd
  real*8 :: f0,f1,fdel,rms,max_d,min_d
  real*8,parameter :: small=1d-10

  Select Case(hard_constraints)
  Case("no")
     call getarg(1,base); call getarg(2,path); 
  Case("yes")
     call getarg(1,base); call getarg(2,path); call getarg(3,consfile);
  End Select

  !***********path************!
  write(*,*) "path input should be path relative to ", absolute_path

  ! get hard_constraints
  if ( hard_constraints .eq. "yes" ) then
     call get_hard_constraints(consfile)
  endif

  call get_molecule_data(path,base)
  call numindatoms(atms)
  allocate(atmtype(atms),link_list(atms))
  call collapseatomtype(h_atms)


  ! get point to point responses for all frequencies
  call get_p2p(path,base)

  ! set local_axis equal to global axis since this is isotropic model
  local_axis=0d0
  do j=1,3
     local_axis(:,j,j)=1d0
  enddo


  ! for isotropic model, all dipole-quadrupole, quadrupole-octapole polarizabilities are zero
  ! also, all dipole-dipole, quadrupole-quadrupole tensors are diagonal with one eigenvalue
  ! therefore, instead of storing a two component tensor with different elements for each rank
  ! we just need one value per rank

  atoms=size(atomtype)
  allocate(dist_pol_tensor(atoms,n_freq+1,rank))

  ! need to fill in hard constraints
  if ( hard_constraints .eq. "yes" ) then
     call fill_in_hard_constraints
  endif

write(*,*) atms
write(*,*) atmtype

! number of parameters depends on rank of hydrogen atoms
  n_param= ( atms - h_atms) * rank + h_atms * rankh

  ! allocate vectors and matrix for linear fitting
  allocate(p0(n_param),G_matrix(n_param,n_param),G_matrix_store(n_param,n_param),f_vector(n_param),w_svd(n_param),v_svd(n_param,n_param))
  p0=1d0


  !***************** G_matrix and static tensor (located in n_freq+1) ********************
  i_freq= n_freq+1
  call link_param(p0)

  write(*,*) "***** initializing distance and direction cosine data****"
  ! initialize data for interaction tensor calculations
  call initialize_elec_int_data

  ! construct G_matrix
  write(*,*) "***** started constructing G_matrix***********"
  call construct_G_matrix(G_matrix)
  ! store G_matrix as svd will overwrite it
  G_matrix_store=G_matrix

  call svdcmp_sp(G_matrix,w_svd,v_svd)
     DO j=1, SIZE(w_svd(:))
        write(*,*) w_svd(j)
!!$        IF(w_svd(j) < small) THEN
!!$           w_svd(j)=0.0
!!$        ENDIF
     ENDDO

  ! only fit static frequency if we are not using hard constraints, since we don't have 
  ! constraints for this value
  Select Case(hard_constraints)
     Case("yes")
  write(*,*) "using hard constraints, therefore skipping frequency 0 "
     Case("no")
  write(*,*) "**** fitting frequency 0 *********"
  ! construct f_vector
  call construct_f_vector(f_vector)

     CALL svbksb_sp(G_matrix,w_svd,v_svd,f_vector,p0)

! fill in dist pol tensor for this frequency
     call link_param(p0)

     write(*,*) p0

  End Select


 !**************************** now frequency dependent polarizabilities *********************
 ! loop over all frequencies, no need to reconstruct G_matrix, but we need to reconstruct f_vector
  do i_freq=1,n_freq
    ! G_matrix=G_matrix_store
  write(*,*) "**** fitting frequency ",i_freq," ********"
    
 call construct_f_vector(f_vector)
    !call svdcmp_sp(G_matrix,w_svd,v_svd)

     CALL svbksb_sp(G_matrix,w_svd,v_svd,f_vector,p0)
! fill in dist pol tensor for this frequency
     write(*,*) p0
     call link_param(p0)
enddo   

!************************************** done *************************************************


! now print
call print_casimir(base)

! print rms for all frequency
do i_freq=1,n_freq
  call rms_p2p(rms,max_d,min_d)
  write(*,*) "rms for frequency ",i_freq, " is ",rms
  write(*,*) "max and min deviations are ",max_d,"   ",min_d, " respectively"
enddo

end program p2p
