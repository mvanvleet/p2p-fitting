
!************************** LINEAR FITTING ***************************************
! these subroutines construct matrices for linear fitting
! for linear fit we have Ga=f, where a is a vector with all the independent polarizabilities
! G is a matrix defined by a triple sum over point to point sites, independent components of the 
! particular tensor rank, and identical atom types over gij=TiTiTjTj where T are the interaction
! tensors, and i and j represent a rank of a particular atom type
! f is a vector defined with a similar sum over fi=TiTi
!********************************************************************************


!**************************************************
! this subroutine constructs the G_matrix
!**************************************************
subroutine construct_G_matrix(G_matrix)
  use variables
  use electrostatic_interaction_tensor
  real*8,dimension(:,:), intent(out) :: G_matrix

  integer :: i,j, k,l,m,n,o,p,atms,atoms,n_points,rankA,rankB,ind,param_A,param_B,mmA,mmB,nnA,nnB,limA,limB
  real*8 :: sum,Tpa,Taq,Tpb,Tbq
  real*8,dimension(:,:),allocatable:: R_AP,R_BP,rap_cos,rbp_cos,rpa_cos,rpb_cos
  real*8,dimension(3,0:3,0:1,2) :: e_int_tensor_pa,e_int_tensor_aq,e_int_tensor_pb,e_int_tensor_bq

  atms=size(atmtype)
  atoms=size(atomtype)
  n_points = size(xyz_pts(:,1))

  allocate( R_AP(n_points,3),rap_cos(n_points,3),rpa_cos(n_points,3) )
  allocate( R_BP(n_points,3),rbp_cos(n_points,3),rpb_cos(n_points,3) )

  G_matrix=0d0

  ! loop over all atomtypes
  do i=1, atms
     ! loop over all atoms of the same type
     do j=1, atoms
        if ( atomtype(j) .eq. atmtype(i) ) then
           ! locally store all the vectors and cosines for this atomtype
           R_AP(:,:)=R_inv_n(:,j,:)
           rap_cos(:,:)=ra_cosine(:,j,:)
           rpa_cos(:,:)=rp_cosine(:,j,:)
           do k=1, atms
              do l = 1, atoms
                 if (atomtype(l) .eq. atmtype(k) ) then
                    ! locally store all the vectors and cosines for this atomtype
                    R_BP(:,:)=R_inv_n(:,l,:)
                    rbp_cos(:,:)=ra_cosine(:,l,:) 
                    rpb_cos(:,:)=rp_cosine(:,l,:)
                    ! loop over all point to point responses
                    do m=1,n_points
                       do n=1,m

                          ! get interaction tensor for all ranks up to 3 (even if we aren't using them)
                          e_int_tensor_pa = elec_int_tensor( R_AP(m,:),rap_cos(m,:),rpa_cos(m,:),rank,1 )
                          e_int_tensor_aq = elec_int_tensor( R_AP(n,:),rap_cos(n,:),rpa_cos(n,:),rank,0 )
                          e_int_tensor_pb = elec_int_tensor( R_BP(m,:),rbp_cos(m,:),rpb_cos(m,:),rank,1 )
                          e_int_tensor_bq = elec_int_tensor( R_BP(n,:),rbp_cos(n,:),rpb_cos(n,:),rank,0 )

                          ! figure out what ranks we are fitting
                          rankA=rank;rankB=rank;
                          ! check for hydrogen
                          ind=index(atmtype(i),"H")
                          if (ind .ne. 0 ) rankA=1
                          ind=index(atmtype(k),"H")
                          if (ind .ne. 0 ) rankB=1

                          ! now loop over all ranks for these sites
                          do o=1,rankA
                             do p=1, rankB
                                ! determine which parameters which array element we're filling in
                                ! remember link_list stores where the rank 1 tensor component starts for
                                ! particular atom type
                                param_A=link_list(i)+(o-1)
                                param_B=link_list(k)+(p-1)

                                ! now loop over all components of the tensor on these sites
                                ! notice this corresponds to three indices each (since diagonal): rankA,mmA,nnA and rankB,mmB,nnB
                                ! the corresponding index for the tensor component for rank2 is given by
                                ! max(1, 2*m + n)
                                do mmA=0,o
                                   limA=min(1,mmA)
                                   do nnA=0,limA
                                      do mmB=0,p
                                         limB=min(1,mmB)
                                         do nnB=0,limB
                                            ! get interaction tensors
                                            Tpa= e_int_tensor_pa(o,mmA,nnA,2)
                                            Taq=e_int_tensor_aq(o,mmA,nnA,1)
                                            Tpb= e_int_tensor_pb(p,mmB,nnB,2)
                                            Tbq=e_int_tensor_bq(p,mmB,nnB,1)

                                            ! add this contribution to G_matrix
                                            G_matrix(param_A,param_B) = G_matrix(param_A,param_B) + Tpa*Taq*Tpb*Tbq

                                         enddo
                                      enddo
                                   enddo
                                enddo   ! end loop over mmA
                             enddo
                          enddo     ! end loop over o (rankA)
                       enddo
                    enddo    ! end loop over m (n_points)

                 endif
              enddo
           enddo ! end loop over k (atms)
        endif
     enddo
  enddo ! end loop over i (atms)

  deallocate( R_AP,rap_cos,rpa_cos )
  deallocate( R_BP,rbp_cos,rpb_cos )

end subroutine construct_G_matrix



!***************************************************
! this subroutine constructs the f_vector
!***************************************************
subroutine construct_f_vector(f_vector)
  use variables
  use electrostatic_interaction_tensor
  real*8,dimension(:),intent(out) :: f_vector
  integer :: i,j, k,l,m,n,o,p,atms,atoms,n_points,rankA,rankB,ind,param_A,param_B,mmA,mmB,nnA,nnB,limA,limB
  real*8 :: sum,Tpa,Taq,Tpb,Tbq
  real*8,dimension(3) :: R_P,R_Q,rap_cos,rpa_cos,raq_cos,rqa_cos
  real*8,dimension(3,0:3,0:1,2) :: e_int_tensor_pa,e_int_tensor_aq,e_int_tensor_pb,e_int_tensor_bq

  atms=size(atmtype)
  atoms=size(atomtype)
  n_points = size(xyz_pts(:,1))


  f_vector=0d0

  ! loop over all point to point responses
  do m=1,n_points
     do n=1,m
        ! first we need to subtract the contributions from the hard constraints from the p2p response
        sum=0d0
        if ( hard_constraints .eq. "yes" ) then
           do k=1, size(hard_cons_type)
              do l =1, atoms
                 if (atomtype(l) .eq. hard_cons_type(k) ) then
                    R_P(:) = R_inv_n(m,l,:)
                    R_Q(:) = R_inv_n(n,l,:)
                    rap_cos(:)=ra_cosine(m,l,:)
                    rpa_cos(:)=rp_cosine(m,l,:)
                    raq_cos(:)=ra_cosine(n,l,:)
                    rqa_cos(:)=rp_cosine(n,l,:)

                     ! get interaction tensor for all ranks
                          e_int_tensor_pb = elec_int_tensor( R_P,rap_cos,rpa_cos,rank,1 )
                          e_int_tensor_bq = elec_int_tensor( R_Q,raq_cos,rqa_cos,rank,0 )
                     ! figure out what ranks we're fitting
                    rankB=rank;
                    ind=index(hard_cons_type(k),"H")
                    if (ind .ne. 0 ) rankB=1
                    do p=1, rankB
                       do mmB=0,p
                          limB=min(1,mmB)
                          do nnB=0,limB
                             ! get interaction tensors
                             Tpb= e_int_tensor_pb(p,mmB,nnB,2)
                             Tbq= e_int_tensor_bq(p,mmB,nnB,1)
!!$                             sum=sum + dist_pol_tensor(l,i_freq,p) * Tpb * Tbq
                             sum=sum - dist_pol_tensor(l,i_freq,p) * Tpb * Tbq
                          enddo
                       enddo
                    enddo
                 endif
              enddo
           enddo
        endif

        ! now let add sum to p2p response
        sum = sum + p2p_response(i_freq,m,n)


        ! loop over all atomtypes
        do i=1, atms
           ! loop over all atoms of the same type
           do j=1, atoms
              if ( atomtype(j) .eq. atmtype(i) ) then
                    R_P(:) = R_inv_n(m,j,:)
                    R_Q(:) = R_inv_n(n,j,:)
                    rap_cos(:)=ra_cosine(m,j,:)
                    rpa_cos(:)=rp_cosine(m,j,:)
                    raq_cos(:)=ra_cosine(n,j,:)
                    rqa_cos(:)=rp_cosine(n,j,:)

                          ! get interaction tensor for all ranks
                          e_int_tensor_pa = elec_int_tensor( R_P,rap_cos,rpa_cos,rank,1 )
                          e_int_tensor_aq = elec_int_tensor( R_Q,raq_cos,rqa_cos,rank,0 )
                 ! figure out what ranks we are fitting
                 rankA=rank
                 ! check for hydrogen
                 ind=index(atmtype(i),"H")
                 if (ind .ne. 0 ) rankA=1

                 do o=1,rankA
                    ! determine which parameters which array element we're filling in
                    ! remember link_list stores where the rank 1 tensor component starts for
                    ! particular atom type
                    param_A=link_list(i)+(o-1)

                    ! now loop over all components of the tensor
                    ! this corresponds to three indices(since diagonal): rankA,mmA,nnA
                    do mmA=0,o
                       limA=min(1,mmA)
                       do nnA=0,limA
                          ! get interaction tensors
                          Tpa= e_int_tensor_pa(o,mmA,nnA,2)
                          Taq=e_int_tensor_aq(o,mmA,nnA,1)

                          ! add contribution to f_vector
                          f_vector(param_A) = f_vector(param_A) + sum * Tpa * Taq
                       enddo
                    enddo
                 enddo
              endif
           enddo
        enddo         ! end loop i over atomtypes
     enddo
  enddo  ! end loop over response points m

end subroutine construct_f_vector



! ************************** NON-LINEAR FITTING**************************************
! this is slow and probably shouldn't be used unless needed for constraint types
!**********************************************************************************

!!$function kaisq(p0)
!!$  use variables
!!$  real*8 :: kaisq
!!$  real*8,dimension(:),intent(in)::p0
!!$
!!$  interface
!!$     function elec_int_tensor(la,ma,csa,lb,mb,csb,rA,rB,locA,locB)
!!$       integer,intent(in) :: la,ma,csa,lb,mb,csb
!!$       real*8,dimension(3),intent(in)::rA,rB
!!$       real*8,dimension(3,3),intent(in)::locA,locB
!!$     end function elec_int_tensor
!!$     subroutine outparam(p0)
!!$       use variables
!!$       real*8,dimension(:),intent(in)::p0
!!$     end subroutine outparam
!!$  end interface
!!$
!!$  integer :: n_points,n_sites,i,j,k,l,m,n,ll,mm,nn,lim1,lim2,index1,index2
!!$  real*8 :: sum,rA(3),rP(3),rQ(3),lA(3,3),lP(3,3),lQ(3,3),Tpa,Taq
!!$
!!$  kaisq=0d0
!!$
!!$  n_points = size(xyz_pts(:,1))
!!$  n_sites=size(site_use)
!!$
!!$  call outparam(p0)
!!$
!!$  ! points for response are always in global axis
!!$  lP=0d0;lQ=0d0
!!$  do i=1,3
!!$     lP(i,i)=1d0
!!$     lQ(i,i)=1d0
!!$  enddo
!!$
!!$  ! loop over all point to point responses
!!$  do i=1,n_points
!!$     rP(:)=xyz_pts(i,:)
!!$     do j=1,i
!!$        rQ(:)=xyz_pts(j,:)
!!$        sum=0d0
!!$        ! loop over all sites where we are fitting polarizability tensors
!!$        do k=1,n_sites
!!$           lA(:,:)=local_axis(k,:,:)
!!$           rA(:)=xyz_molecule(k,:)
!!$           ! if this site is being used
!!$           if (site_use(k) .eq. 1) then
!!$              ! now loop over all components of the tensor on this site
!!$              ! notice this corresponds to six indices: l,m,n and ll,mm,nn
!!$              ! the corresponding index for the tensor component for rank2 is given by
!!$              ! max(1, 2*m + n)
!!$              l=rank
!!$              ll=rank
!!$              do m=0,l
!!$                 lim1=min(1,m)
!!$                 do n=0,lim1
!!$                    index1=max(1,2*m+n)
!!$                    do mm=0,ll
!!$                       lim2=min(1,mm)
!!$                       do nn=0,lim2
!!$                          index2=max(1,2*mm+nn)
!!$                          ! get interaction tensors
!!$                          Tpa= elec_int_tensor(0,0,0,l,m,n,rP,rA,lP,lA)
!!$                          Taq=elec_int_tensor(ll,mm,nn,0,0,0,rA,rQ,lA,lQ)
!!$                          ! now add contribution to the response from this tensor component
!!$                          sum=sum- Tpa* dist_pol_tensor(k,rank,index1,rank,index2)* Taq
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           endif
!!$        enddo
!!$
!!$        ! finish looping over all sites, "sum" now is the total response with contributions from
!!$        ! all the tensor components on all the sites
!!$        ! add to kaisq
!!$
!!$        kaisq=kaisq + ( sum - p2p_response(i_freq,i,j))**2
!!$     enddo
!!$  enddo
!!$
!!$end function kaisq
!!$
!!$
!!$
!!$
!!$function dkaisq(p0)
!!$  use variables
!!$  real*8,dimension(:),intent(in)::p0
!!$  real*8,dimension(size(p0)) :: dkaisq
!!$
!!$  interface
!!$     function elec_int_tensor(la,ma,csa,lb,mb,csb,rA,rB,locA,locB)
!!$       integer,intent(in) :: la,ma,csa,lb,mb,csb
!!$       real*8,dimension(3),intent(in)::rA,rB
!!$       real*8,dimension(3,3),intent(in)::locA,locB
!!$     end function elec_int_tensor
!!$  end interface
!!$
!!$  integer :: n_points,n_sites,i,j,k,l,m,n,ll,mm,nn,lim1,lim2,index1,index2,count
!!$  real*8 :: sum,rA(3),rP(3),rQ(3),lA(3,3),lP(3,3),lQ(3,3),Tpa,Taq,term1
!!$
!!$  dkaisq=0d0
!!$
!!$  n_points = size(xyz_pts(:,1))
!!$  n_sites=size(site_use)
!!$
!!$
!!$  ! points for response are always in global axis
!!$  lP=0d0;lQ=0d0
!!$  do i=1,3
!!$     lP(i,i)=1d0
!!$     lQ(i,i)=1d0
!!$  enddo
!!$
!!$  ! loop over all point to point responses
!!$  do i=1,n_points
!!$     rP(:)=xyz_pts(i,:)
!!$     do j=1,i
!!$        rQ(:)=xyz_pts(j,:)
!!$        sum=0d0
!!$        ! loop over all sites where we are fitting polarizability tensors
!!$        do k=1,n_sites
!!$           lA(:,:)=local_axis(k,:,:)
!!$           rA(:)=xyz_molecule(k,:)
!!$           ! if this site is being used
!!$           if (site_use(k) .eq. 1) then
!!$              ! now loop over all components of the tensor on this site
!!$              ! notice this corresponds to six indices: l,m,n and ll,mm,nn
!!$              ! the corresponding index for the tensor component for rank2 is given by
!!$              ! max(1, 2*m + n)
!!$              l=rank
!!$              ll=rank
!!$              do m=0,l
!!$                 lim1=min(1,m)
!!$                 do n=0,lim1
!!$                    index1=max(1,2*m+n)
!!$                    do mm=0,ll
!!$                       lim2=min(1,mm)
!!$                       do nn=0,lim2
!!$                          index2=max(1,2*mm+nn)
!!$                          ! get interaction tensors
!!$                          Tpa= elec_int_tensor(0,0,0,l,m,n,rP,rA,lP,lA)
!!$                          Taq=elec_int_tensor(ll,mm,nn,0,0,0,rA,rQ,lA,lQ)
!!$                          ! now add contribution to the response from this tensor component
!!$                          sum=sum- Tpa* dist_pol_tensor(k,rank,index1,rank,index2)* Taq
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$           endif
!!$        enddo
!!$
!!$        ! finish looping over all sites, "sum" now is the total response with contributions from
!!$        ! all the tensor components on all the sites
!!$
!!$
!!$        term1= 2d0* ( sum - p2p_response(i_freq,i,j))
!!$
!!$        ! now loop again over all polarizabilities to get derivatives
!!$
!!$        count=0
!!$        do k=1,n_sites
!!$           lA(:,:)=local_axis(k,:,:)
!!$           rA(:)=xyz_molecule(k,:)
!!$           ! if this site is being used
!!$           if (site_use(k) .eq. 1) then
!!$              ! now loop over all components of the tensor on this site
!!$              ! notice this corresponds to six indices: l,m,n and ll,mm,nn
!!$              ! the corresponding index for the tensor component for rank2 is given by
!!$              ! max(1, 2*m + n)
!!$              l=rank
!!$              ll=rank
!!$              do m=0,l
!!$                 lim1=min(1,m)
!!$                 do n=0,lim1
!!$                    index1=max(1,2*m+n)
!!$                    do mm=0,ll
!!$                       lim2=min(1,mm)
!!$                       do nn=0,lim2
!!$                          index2=max(1,2*mm+nn)
!!$                          ! get interaction tensors
!!$                          Tpa= elec_int_tensor(0,0,0,l,m,n,rP,rA,lP,lA)
!!$                          Taq=elec_int_tensor(ll,mm,nn,0,0,0,rA,rQ,lA,lQ)
!!$
!!$                          ! add to derivative 
!!$                          ! for rank 1 tensors only, each site has 9 parameters
!!$                          ! so corresponding index for derivative array is
!!$                          ! count*9 + (index1-1)*3 + index2
!!$
!!$                          dkaisq(count*9 +(index1-1)*3 +index2) = dkaisq(count*9 +(index1-1)*3 +index2)- Tpa* Taq * term1
!!$                       enddo
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$              count=count+1
!!$           endif
!!$        enddo
!!$
!!$
!!$     enddo
!!$  enddo
!!$
!!$end function dkaisq
