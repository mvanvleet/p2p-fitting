module variables


  !********************* ABSOLUTE PATH --might need to change this
  character(80) :: absolute_path="/home/mvanvleet/"
  !***************************************************************
  

  integer,parameter :: MAX_CONS_TYPE=30
  real*8,dimension(:,:),allocatable:: xyz_pts,xyz_molecule
  real*8,dimension(:,:,:),allocatable:: p2p_response,local_axis
  character(3),dimension(:),allocatable::atomtype,atom_label,atmtype
  integer,dimension(:),allocatable :: link_list
  integer,parameter:: n_freq=10        ! number of imag freq.  10 is default for camcasp
  integer::i_freq ! frequency that we are currently fitting (0 means static)
  integer,parameter::rank=3, rankh=3
  character(3) :: hard_constraints="yes"  ! if we are plugging in definite values for pol tensors 
  real*8,dimension(:,:,:),allocatable :: dist_pol_tensor
  character(3),dimension(MAX_CONS_TYPE) :: hard_cons_type
  real*8,dimension(MAX_CONS_TYPE,n_freq,rank) :: hard_cons_values
  
  !*************** stored data for calculating interaction tensor values **************
  real*8,dimension(:,:,:),allocatable :: R_inv_n, ra_cosine , rp_cosine

end module variables
