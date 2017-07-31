module routines

contains

!***************************************************
! this subroutine reads a .clt file and obtains the molecule coordinates and atom types
!**********************************************
subroutine get_molecule_data(path,base)
    use variables
    character(*),intent(in)::base,path
    character(160) :: molecule_file
    character(50) :: line
    integer::inputstatus,ind,count,i,j
    real*8 :: junkr
    character(4) :: junkc

  molecule_file=trim(absolute_path)//trim(path)//trim(base)//".clt"

 open(unit=7,file=molecule_file,status='old')

    ! get number of atom types

    do
       read(7,'(A)',Iostat=inputstatus) line
       if(inputstatus < 0 ) exit
       ind=index(line,'Charge')
       if(ind .ne.0) exit
    enddo

    count=0
    do
       read(7,'(A)',Iostat=inputstatus) line
       if(inputstatus < 0 ) exit
       ind=index(line,'End')
       if(ind .ne.0) exit
       count=count+1
    enddo

    allocate(atomtype(count),atom_label(count),xyz_molecule(count,3),local_axis(count,3,3))

    do i=1,count+1
       backspace(7)
    enddo

    ! fill in atom types
    do i=1,count
       read(7,*) atom_label(i),junkr,(xyz_molecule(i,j),j=1,3),junkc,atomtype(i)
    enddo

close(7)

end subroutine get_molecule_data

!*****************************************************
! this subroutine reads the positions and response values for the
! point to point response sites from the camcasp .p2p binary files
!******************************************************
  subroutine get_p2p(path,base)
    use variables
    character(*),intent(in)::base,path
  character(160) :: response_file
    character(1)  :: temp1
    character(2)  :: temp2
 integer*4 :: n_points
  integer:: i,j,k

! first open the static file
response_file=trim(absolute_path)//trim(path)//trim(base)//'/'//trim(base)//"_000.p2p"

open (unit=7,file=response_file,form='unformatted',access='sequential')

read(7) n_points

allocate(xyz_pts(n_points,3),p2p_response(n_freq+1,n_points,n_points))

do i=1,n_points
 read(7) ( xyz_pts(i,j),j=1,3)
enddo

do i=1,n_points
   do j=1, i
      read(7) p2p_response(n_freq+1,i,j)
   enddo
enddo

close(7)
! now read dynamic pol files

do k=1,n_freq

if(k < 10 ) then
write(temp1,'(I1)') k
response_file=trim(absolute_path)//trim(path)//trim(base)//'/'//trim(base)//"_00"//temp1//".p2p"
else
write(temp2,'(I2)') k
response_file=trim(absolute_path)//trim(path)//trim(base)//'/'//trim(base)//"_0"//temp2//".p2p"
endif

open (unit=7,file=response_file,form='unformatted',access='sequential')

read(7) n_points

do i=1,n_points
 read(7) ( xyz_pts(i,j),j=1,3)
enddo

do i=1,n_points
   do j=1, i
      read(7) p2p_response(k,i,j)
   enddo
enddo

close(7)

enddo

end subroutine get_p2p


!****************************************************
! this subroutine reads in frequency dependent polarizability
! values that will be hard constraints
!****************************************************
subroutine get_hard_constraints(ifile)
  use variables
  character(*),intent(in) :: ifile
  character(50) :: line
  integer :: count,i,j,rank_l


open(unit=7,file=ifile,status='old')

       read(7,'(A)') line
       read(7,'(A)') line
       read(7,'(A)') line
       
       count=1
       ! keep reading hard constraint parameters until end of file
    do
       read(7,*,Iostat=inputstatus) hard_cons_type(count)
       if(inputstatus < 0 ) exit
    !   do i=1,rank
       i=1
          read(7,*) rank_l
          read(7,*) hard_cons_values(count,n_freq+1,rank_l)
!!$          read(7,*) (hard_cons_values(count,j,rank_l),j=1,5 )
!!$          read(7,*) (hard_cons_values(count,j,rank_l),j=6,10 )
   !    enddo
       count=count+1
       if (count .ge. MAX_CONS_TYPE ) then
          stop "more atom types in hard constraint file than MAX_CONS_TYPE"
       endif
    enddo

close(7)

end subroutine get_hard_constraints


!******************************************************
! this subroutine fills in hard values for the distributed 
! polarizability tensors of atoms that match labels in
! the hard_constraint input file
!******************************************************
subroutine fill_in_hard_constraints
  use variables
  integer :: i,j,k

  do i=1,size(atomtype)
     do j=1,size(hard_cons_type)
        if (atomtype(i) .eq. hard_cons_type(j) ) then
           ! remember, dist_pol_tensor has an index for static frequency
           ! but hard_cons_values doesn't, so we have to loop here
           ! remember, parameters can be negative if we want to constrain them
           ! so take absolute value here

           do k=1,n_freq+1
           dist_pol_tensor(i,k,:) = abs(hard_cons_values(j,k,:))
           enddo
        endif
     enddo
  enddo

end subroutine fill_in_hard_constraints


!******************************************************
! this subroutine determines the number of independent atomtypes
! as well as the number of hydrogen atoms
!
! if parameter is initially set to zero or negative value, don't fit it
! if it's negative use absolute value of that number as a hard constraint
!******************************************************
SUBROUTINE numindatoms(atms)
  use variables
  INTEGER,INTENT(out)::atms
  CHARACTER(2)::value
  INTEGER::count1,count2,atoms,count,i,j,ind
  integer,dimension(MAX_CONS_TYPE):: hard_cons_zero
  CHARACTER(3),DIMENSION(:),ALLOCATABLE::typemon1

  atoms=SIZE(atomtype)
  hard_cons_zero=1

  ALLOCATE(typemon1(atoms))

  DO i=1,atoms
     WRITE(typemon1(i),'(A3)') atomtype(i)
  ENDDO

  count1=0
  count2=0
  DO i=1,atoms
     if ( hard_constraints .eq. "yes" ) then
!!!!! fit atom type if it is initially greater than zero
        do j=1,size(hard_cons_type)
           if(atomtype(i).eq.hard_cons_type(j)) then
              if ( hard_cons_values(j,n_freq+1,1) > 0d0 ) then
                 hard_cons_zero(j)=0
              endif
           endif
        enddo
     endif
  ENDDO


  ! number of atom types fit is number of non zero initial values
  if ( hard_constraints .eq. "yes" ) then
     do j=1,size(hard_cons_type)
        if ( hard_cons_zero(j) .eq. 0 ) then
           count2=count2 + 1
        endif
     enddo
  endif

  atms=count2


END SUBROUTINE numindatoms


 subroutine collapseatomtype(h_atms)
    use variables
    INTEGER,INTENT(out)::h_atms
    INTEGER::count,countt,count1,count2,atoms,atms,i,j,k,ind

    atoms=SIZE(atomtype)
    atms=SIZE(atmtype)

    count=0
    countt=0
    count1=0

    DO i=1,atms
       count=countt
       DO j=i+count,atoms
!!! check if we are fitting this atom type
          if ( hard_constraints .eq. "yes" ) then
             do k=1,size(hard_cons_type)
                if(atomtype(j).eq.hard_cons_type(k)) then
                   ! fitting if non zero initial value
                   if ( hard_cons_values(k,n_freq+1,1) > 0d0 ) then
                      ! nothing
                   else
                   countt=countt+1
                   go to 202 
                 endif
                endif
             enddo
          endif
          DO k=1,j-1
             IF(atomtype(k).EQ.atomtype(j)) THEN
                countt=countt+1
                go to 202
             ENDIF
          ENDDO
          atmtype(i)=atomtype(j)
          !see if this is a hydrogen
          ind=index(atomtype(j),"H")
          if (ind .ne. 0 ) count1=count1+1

          go to 203
202       CONTINUE
       ENDDO
203    CONTINUE
    ENDDO

    h_atms = count1 


  END SUBROUTINE collapseatomtype


!********************************************
! this subroutine calculates all distances, direction cosines which are needed
! for calculating all electrostatic interaction tensors of interest
! notation is consistent with The Theory of Intermolecular Forces, AJ Stone
!********************************************
  subroutine initialize_elec_int_data
    use variables
    integer :: n_points,n_atom,i,j,k
    real*8,dimension(3) :: rP , rA, Rap,ap_unit,ra_cos,rp_cos
    real*8 :: Rinv,Rinv2,Rinv3,Rinv4
    real*8,dimension(3,3) :: lP,lA

    n_points=size(p2p_response(1,:,1))
    n_atom=size(atomtype)

    allocate( R_inv_n(n_points,n_atom,3) , ra_cosine(n_points,n_atom,3) , rp_cosine(n_points,n_atom,3))

    ! points for response are always in global axis
    lP=0d0
    do i=1,3
       lP(i,i)=1d0
    enddo

    do i=1,n_points
       rP(:)=xyz_pts(i,:)
       do j=1,n_atom
          lA(:,:)=local_axis(j,:,:)
          rA(:)=xyz_molecule(j,:)

          ! keeping notation consistent with Stone, Rap is vector from a to p
          Rap=rP-rA
          Rinv2=1d0/dot_product(Rap,Rap)
          Rinv= sqrt(Rinv2)
          ap_unit=Rap*Rinv
          Rinv3 = Rinv*Rinv2
          Rinv4 = Rinv2*Rinv2

          ! compute direction cosines

          do k=1,3
             ra_cos(k) = dot_product(lA(k,:),ap_unit(:))
             rp_cos(k) = -dot_product(lP(k,:),ap_unit(:))
          enddo

          ! now store in global variables
          R_inv_n(i,j,1)=Rinv2;R_inv_n(i,j,2)=Rinv3;R_inv_n(i,j,3)=Rinv4
          ra_cosine(i,j,:) = ra_cos(:)
          rp_cosine(i,j,:) = rp_cos(:)
       enddo
    enddo


  end subroutine initialize_elec_int_data

  
 

!***************************************
!  this subroutine links_parameters to tensor components and fills in 
!  polarizability tensor
!****************************************

subroutine link_param(p0)
use variables
  real*8,dimension(:),intent(in)::p0
integer :: i,j,k,count,atms,atoms,ind

atoms=size(atomtype)
atms=size(atmtype)

count=1

! parameters are filled in order of atmtype, in order of rank
do i=1,atms
  ! fill in link_list which stores where atmtype rank 1 tensor element starts
   ! in the parameter vector
   link_list(i)=count

   ! check if this is a hydrogen atom
   ind=index(atmtype(i),"H")
  Select Case(ind)
     Case(0)
        ! non hydrogen, fill in all ranks
        do j=1,atoms
           if (atomtype(j) .eq. atmtype(i) ) then
              do k=1,rank
              dist_pol_tensor(j,i_freq,k)=p0(count+(k-1))
              enddo
           endif
        enddo
      count=count+rank
   case  default
      ! hydrogen, fill in first rank
      do j=1,atoms
         if (atomtype(j) .eq. atmtype(i) ) then
            dist_pol_tensor(j,i_freq,1)=p0(count)
            if (rank .gt. 1 ) then
               do k=2,rank
                  dist_pol_tensor(j,i_freq,k)=0d0
               enddo
            endif
         endif
      enddo
      count=count+1

end select

enddo

end subroutine link_param


subroutine fill_charges(param)
 use variables
 real*8,dimension(:),intent(in) :: param

integer :: i,j,atms,atoms

atms=size(atmtype)
atoms=size(atomtype)

do i=1, atoms
   do j=1, atms

if ( atmtype(j) .eq. atomtype(i) ) then
  shellcharge(i)= -abs(param(j))
endif
  enddo
enddo

end subroutine fill_charges



!**************************************************
! this subroutine prints all of the frequency dependent polarizabilities for each atom
! type in the format of a .casimir  output file so that we can use casimir (CAMCASP)
! to construct the C6...C10 coefficients from these polarizabilities
!**************************************************

subroutine print_casimir(base)
  use variables
character(*),intent(in) :: base
character(60) :: ofile
character(1) :: temp1
integer :: i,j

write(temp1,'(I1)') rank

ofile=trim(base)//"_ref_wt0_L"//temp1//".casimir"

open(unit=6,file=ofile,status="new")

write(6,*) "Title molecular frequencies"
write(6,*) "Frequencies   0.5    10"
write(6,*) "Skip  0"
write(6,*) "Print nonzero"
write(6,*)
write(6,*) "Molecule  ",base

do i=1,size(atomtype)
write(6,*) "Site  ",atom_label(i)," type  ",atomtype(i)
do j=1,rank
  Select Case(j)
Case(1)
write(6,*) "10 10"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "11c 11c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "11s 11s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
Case(2)
write(6,*) "20 20"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "21c 21c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "21s 21s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "22c 22c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "22s 22s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
Case(3)
write(6,*) "30 30"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "31c 31c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "31s 31s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "32c 32c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "32s 32s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "33c 33c"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
write(6,*) "33s 33s"
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=1,5)
write(6,'(3x,6F16.8)') (dist_pol_tensor(i,k,j),k=6,10)
end Select
enddo
write(6,*) "End"
enddo

write(6,*)
write(6,*) "CGdir  /home/jmcdaniel/apps/camcasp-5.5/data/realcg"
write(6,*) "Dispersion   12 ",base
write(6,*)
write(6,*) "Finish"

end subroutine print_casimir




end module routines
