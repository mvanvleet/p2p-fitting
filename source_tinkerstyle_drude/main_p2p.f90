!*****************************************************************
! this program is for fitting polarizability tensors (freq dependent and static) to the
! point to point response calculations performed by camcasp
! this only fits isotropic diagonal tensor components, as this is needed for 
! a purely isotropic C6,C8,C10 treatment
!*****************************************************************

program p2p
  use variables
  use routines
  use eq_drude
  use amoeba_routines
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
  character(80) :: path,consfile,thole_string
  character(8)::date
  character(10)::time
  integer:: i,j,k,n_param,atoms,atms,h_atms,iter
  real*8,dimension(:),allocatable::p0,p1,df0,delta,f_vector,w_svd
  real*8,dimension(:,:),allocatable::G_matrix,G_matrix_store,v_svd
  real*8 :: f0,f1,fdel,rms,max_d,min_d
  real*8,parameter :: small=1d-10
  real*8:: r_tol=1D-4
  real*8,dimension(:),allocatable:: y_amoeba
  real*8,dimension(:,:),allocatable::p_amoeba

  Select Case(hard_constraints)
  Case("no")
     call getarg(1,base); call getarg(2,path); call getarg(3,thole_string);
  Case("yes")
     call getarg(1,base); call getarg(2,path); call getarg(3,consfile); call getarg(4,thole_string);
  End Select

     read (thole_string, *) thole
     write(*,*) 'thole param = ',thole

  !***********path************!
  write(*,*) "path input should be path relative to ", absolute_path

  ! get hard_constraints
  if ( hard_constraints .eq. "yes" ) then
     call get_hard_constraints(consfile)
  endif

  call get_molecule_data(path,base)
  ! atom types will be fit if initial values are non zero, and will not be fit if initial values are zero

  atoms=size(atomtype)
  allocate(dist_pol_tensor(atoms,n_freq+1,rank))

  ! need to fill in hard constraints
  if ( hard_constraints .eq. "yes" ) then
     call fill_in_hard_constraints
  endif

  call numindatoms(atms)
  allocate(atmtype(atms),link_list(atms))
  call collapseatomtype(h_atms)


  ! get point to point responses for all frequencies
  call get_p2p(path,base)

  write(*,*) "finished reading p2p"

  ! set local_axis equal to global axis since this is isotropic model
  local_axis=0d0
  do j=1,3
     local_axis(:,j,j)=1d0
  enddo



  i_freq= n_freq+1
  ! initialize drude oscillators
  call initialize_drude_oscillators
  call rms_p2p(rms,max_d,min_d)

write(*,*) shellcharge
write(*,*) rms


 ! initialize amoeba data
  allocate( p_amoeba(atms+1,atms), y_amoeba(atms+1) )
  do i=1,atms
     do j=1,atoms
        if ( atomtype(j) .eq. atmtype(i) ) then
  p_amoeba(:,i) = shellcharge(j)
      endif
     enddo
  enddo
 y_amoeba(1) = rms


!! write(*,*) 'stopping'
!! stop

  do i= 2, atms+1
       p_amoeba(i,i-1) = 1.5d0 * p_amoeba(i,i-1)
       call fill_charges(p_amoeba(i,:))
       call equilibrate_drude_positions
       call rms_p2p(rms,max_d,min_d)
       y_amoeba(i)=rms
  
  !! write(*,*) 'amoeba params:', p_amoeba
  !! write(*,*) rms
  !! stop
  
  enddo

  !! Main function call
  call amoeba(p_amoeba,y_amoeba,r_tol,iter) 


  call rms_p2p(rms,max_d,min_d)

write(*,*) "parameters"
do i=1,atms
   do j=1,atoms
      if ( atomtype(j) .eq. atmtype(i) ) then  
     write(*,*) atmtype(i), shellcharge(j)
     endif
   enddo
enddo
write(*,*) "rms, max_d,min_d"
write(*,*) rms, max_d,min_d
stop

write(*,*) atms
write(*,*) atmtype
  !number of parameters is number of non hydrogen atoms we are fitting times the rank, plus one parameter for each hydrogen, since we only include one rank for hydrogen

  n_param= ( atms - h_atms) * rank + h_atms

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


  write(*,*) "**** fitting frequency 0 *********"
  ! construct f_vector
  call construct_f_vector(f_vector)

     CALL svbksb_sp(G_matrix,w_svd,v_svd,f_vector,p0)

! fill in dist pol tensor for this frequency
     call link_param(p0)

     write(*,*) p0



!!$ !**************************** now frequency dependent polarizabilities *********************
!!$ ! loop over all frequencies, no need to reconstruct G_matrix, but we need to reconstruct f_vector
!!$  do i_freq=1,n_freq
!!$    ! G_matrix=G_matrix_store
!!$  write(*,*) "**** fitting frequency ",i_freq," ********"
!!$    
!!$ call construct_f_vector(f_vector)
!!$    !call svdcmp_sp(G_matrix,w_svd,v_svd)
!!$
!!$     CALL svbksb_sp(G_matrix,w_svd,v_svd,f_vector,p0)
!!$! fill in dist pol tensor for this frequency
!!$     write(*,*) p0
!!$     call link_param(p0)
!!$enddo   

!************************************** done *************************************************


! now print
!!$call print_casimir(base)

! print rms for all frequency
!!$do i_freq=1,n_freq
  call rms_p2p(rms,max_d,min_d)
  write(*,*) "rms for frequency ",i_freq, " is ",rms
  write(*,*) "max and min deviations are ",max_d,"   ",min_d, " respectively"
!!$enddo

end program p2p
