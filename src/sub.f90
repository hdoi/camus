subroutine remove_flow(natom, velocity)
  implicit none
  ! global
  integer natom
  double precision velocity(natom, 3)
  ! local
  integer ii
  double precision c_mass_velocity(3)

  do ii = 1, 3
    c_mass_velocity(ii) = sum( velocity(1:natom,ii)) / dble(natom)
    velocity(:,ii) = velocity(:,ii) - c_mass_velocity(ii)
  enddo
end subroutine

subroutine output_xyz(natom, part_type, coord, box_length, ene)
  use energy_module
  use filenames
  implicit none
  ! global
  integer natom, part_type(natom)
  double precision coord(natom,3), box_length(3)
  type ( energy_var ) ene
  ! local
  integer ii
  double precision half_length(1:3), temp_coord(3)
  character(3) atom_list(20)

  atom_list(1)  = "H  "
  atom_list(2)  = "He "
  atom_list(3)  = "Li "
  atom_list(4)  = "Be "
  atom_list(5)  = "B  "
  atom_list(6)  = "C  "
  atom_list(7)  = "N  "
  atom_list(8)  = "O  "
  atom_list(9)  = "F  "
  atom_list(10) = "Ne "
  atom_list(11) = "Na "
  atom_list(12) = "Mg "
  atom_list(13) = "Al "
  atom_list(14) = "Si "
  atom_list(15) = "P  "
  atom_list(16) = "S  "
  atom_list(17) = "Cl "
  atom_list(18) = "Ar "
  atom_list(19) = "K  "
  atom_list(20) = "Ca "


  half_length = box_length / 2.0d0

  open(30,file=xyz_file)

  write(30,"(i10)") natom+8
  write(30,"(3f15.5)") ene%pot(1), ene%kin, ene%pot(1) + ene%kin
  write(30,"(a3,3f15.3)") "Xe " ,   half_length(1) ,   half_length(2) ,  half_length(3)
  write(30,"(a3,3f15.3)") "Xe " ,   half_length(1) , - half_length(2) ,  half_length(3)
  write(30,"(a3,3f15.3)") "Xe " , - half_length(1) ,   half_length(2) ,  half_length(3)
  write(30,"(a3,3f15.3)") "Xe " , - half_length(1) , - half_length(2) ,  half_length(3)
  write(30,"(a3,3f15.3)") "Xe " ,   half_length(1) ,   half_length(2) ,- half_length(3)
  write(30,"(a3,3f15.3)") "Xe " ,   half_length(1) , - half_length(2) ,- half_length(3)
  write(30,"(a3,3f15.3)") "Xe " , - half_length(1) ,   half_length(2) ,- half_length(3)
  write(30,"(a3,3f15.3)") "Xe " , - half_length(1) , - half_length(2) ,- half_length(3)
  do ii = 1, natom
    temp_coord(:) = coord(ii,1:3)
    write(30,"(a3,3f15.3)") atom_list(part_type(ii)) , temp_coord(1) , temp_coord(2), temp_coord(3)
  enddo
end subroutine

subroutine check_bounds(natom, coord, box_length)
  implicit none
  ! global
  integer natom
  double precision coord(natom,3), box_length(3)
  ! local
  integer ii

  ! molecular moving
  do ii = 1, natom
    coord(ii,1:3) = coord(ii,1:3) - nint(coord(ii,1:3)/box_length(1:3))*box_length(1:3)
  enddo
end subroutine

subroutine divide_into_cell(natom, coord, box_length)
  use cell_index_param
  implicit none
  ! global
  integer natom
  double precision box_length(3)
  double precision coord(natom,3)
  ! local
  integer ii, ilocate,  max_counter
  integer i_index(3)
  double precision half_length(3), grid(3)

  half_length(1:3) =  box_length(1:3)/2.0d0
  grid(1:3) = ( box_length(1:3) / cell_size(1:3) )
  cell_counter(:,:,:) = 0
  ! calc max_counter
  do ii = 1, natom
    i_index(1:3) = int( (coord(ii,1:3) + half_length(1:3))/grid(1:3) ) + 1
    cell_counter(i_index(1), i_index(2), i_index(3)) = cell_counter(i_index(1), i_index(2), i_index(3)) + 1
  enddo

  max_counter = maxval( cell_counter )
  allocate( cell_list(0:cell_size(1)+1,0:cell_size(2)+1,0:cell_size(3)+1,max_counter+10) )

  ! calc cell_index
  cell_counter(:,:,:) = 0
  cell_list(:,:,:,:) = 0
  do ii = 1, natom
    i_index(1:3) = int( (coord(ii,1:3) + half_length(1:3))/grid(1:3) ) + 1
    cell_counter(i_index(1), i_index(2), i_index(3)) = cell_counter(i_index(1), i_index(2), i_index(3)) + 1
    ilocate = cell_counter(i_index(1), i_index(2), i_index(3))
    cell_list(i_index(1), i_index(2), i_index(3), ilocate ) = ii
  enddo

  ! x copy
  cell_counter(0,:,:)              = cell_counter(cell_size(1)  ,:,:)
  cell_counter(cell_size(1)+1,:,:) = cell_counter(1,:,:)
  cell_list(0,:,:,:)              = cell_list(cell_size(1)  ,:,:,:)
  cell_list(cell_size(1)+1,:,:,:) = cell_list(1,:,:,:)
  ! y copy
  cell_counter(:,0,:)              = cell_counter(:,cell_size(2)  ,:)
  cell_counter(:,cell_size(2)+1,:) = cell_counter(:,1,:)
  cell_list(:,0,:,:)              = cell_list(:,cell_size(2)  ,:,:)
  cell_list(:,cell_size(2)+1,:,:) = cell_list(:,1,:,:)
  ! z copy
  cell_counter(:,:,0)              = cell_counter(:,:,cell_size(3))
  cell_counter(:,:,cell_size(3)+1) = cell_counter(:,:,1)
  cell_list(:,:,0,:)              = cell_list(:,:,cell_size(3)  ,:)
  cell_list(:,:,cell_size(3)+1,:) = cell_list(:,:,1,:)

end subroutine
    
subroutine calc_force_cell_index(natom, part_type, coord, velocity, for_dpd, box_length, Epot)
  use cell_index_param
  use dpd_param, only : aij_list
  implicit none
  ! global
  integer natom
  double precision coord(natom,3), velocity(natom,3), box_length(3)
  integer part_type(natom)
  double precision for_dpd(natom,3)
  double precision Epot
  ! local
  integer ii, ij
  double precision dr(3), dv(3), offset(3)
  double precision f(3), half_length(3)
  double precision aij, aij_half_cutoff
  integer ix1,iy1,iz1,ix2,iy2,iz2, cs_x, cs_y, cs_z
  integer iatom1, iatom2, itype1, itype2
  integer myid, omp_get_thread_num
  double precision E_each
  integer linked_cell(13, 3), icell

  call check_bounds(natom, coord, box_length)
  call divide_into_cell(natom, coord, box_length)

  Epot = 0.0d0
  for_dpd     = 0.0d0
  half_length(:) = box_length(:)/2.0d0

  cs_x = cell_size(1)
  cs_y = cell_size(2)
  cs_z = cell_size(3)

  !linked_cell(1)  = ( 1  , 0 , 0)
  !linked_cell(2)  = ( 1  , 1 , 0)
  !linked_cell(3)  = ( 0  , 1 , 0)
  !linked_cell(4)  = ( -1 , 1 , 0)
  !linked_cell(5)  = ( 1  , 0 , -1)
  !linked_cell(6)  = ( 1  , 1 , -1)
  !linked_cell(7)  = ( 0  , 1 , -1)
  !linked_cell(8)  = ( -1 , 1 , -1)
  !linked_cell(9)  = ( 1  , 0 , 1)
  !linked_cell(10) = ( 1  , 1 , 1)
  !linked_cell(11) = ( 0  , 1 , 1)
  !linked_cell(12) = ( -1 , 1 , 1)
  !linked_cell(13) = ( 0  , 0 , 1)
  
  ! x
  linked_cell(1:4,3)  = 0
  linked_cell(5:8,3)  =-1
  linked_cell(9:13,3) = 1
  ! y
  linked_cell(1:13,2) = 1
  linked_cell(1,2) = 0
  linked_cell(5,2) = 0
  linked_cell(9,2) = 0
  linked_cell(13,2) = 0
  ! z
  linked_cell(1:13,1) = 1
  linked_cell(3,1) = 0
  linked_cell(7,1) = 0
  linked_cell(11,1) = 0
  linked_cell(13,1) = 0
  linked_cell(4,1) = -1
  linked_cell(8,1) = -1
  linked_cell(12,1) = -1
 
  !$OMP  parallel do &
  !$OMP& default(none) private(icell, myid, ix1, iy1, iz1, ix2, iy2, iz2, &
  !$OMP& ii,ij, iatom1, iatom2, itype1, itype2, dr, dv, f, E_each, aij, aij_half_cutoff, offset) &
  !$OMP& shared(natom, box_length, cell_counter, cell_list, velocity, coord, &
  !$OMP& cs_x, cs_y, cs_z, half_length, part_type) &
  !$OMP& firstprivate(aij_list, linked_cell) &
  !$OMP& reduction(+:for_dpd, Epot)
  do ix1 = 1, cs_x
  myid = omp_get_thread_num() + 1
  do iy1 = 1, cs_y
  do iz1 = 1, cs_z
    if( cell_counter(ix1,iy1,iz1) .eq. 0) cycle
    do icell = 1, 13
      ix2 = ix1 + linked_cell(icell,1)
      iy2 = iy1 + linked_cell(icell,2)
      iz2 = iz1 + linked_cell(icell,3)
      offset(:) = 0.0d0
      if(ix2 .eq. 0) offset(1) = -box_length(1)
      if(ix2 .eq. cs_x+1) offset(1) = box_length(1)
      if(iy2 .eq. 0) offset(2) = -box_length(2)
      if(iy2 .eq. cs_y+1) offset(2) = box_length(2)
      if(iz2 .eq. 0) offset(3) = -box_length(3)
      if(iz2 .eq. cs_z+1) offset(3) = box_length(3)
      if( cell_counter(ix2,iy2,iz2) .eq. 0) cycle
        
      do ii = 1, cell_counter(ix1,iy1,iz1)
        iatom1 = cell_list(ix1,iy1,iz1,ii)
        itype1 = part_type(iatom1)
        do ij = 1, cell_counter(ix2,iy2,iz2)
          iatom2 = cell_list(ix2,iy2,iz2,ij)
          itype2 = part_type(iatom2)
          ! calc parameter
          aij = aij_list(1, itype1, itype2)
          aij_half_cutoff = aij_list(2, itype1, itype2)
          ! calc force
          dr(1:3) = coord(iatom1,1:3) - (coord(iatom2,1:3) + offset(:))
          dv(1:3) = velocity(iatom1,1:3) - velocity(iatom2,1:3)
          call calc_force_body(dr, dv, f, E_each, myid, aij, aij_half_cutoff)
          for_dpd(iatom1,1:3) = for_dpd(iatom1,1:3) + f(1:3)
          for_dpd(iatom2,1:3) = for_dpd(iatom2,1:3) - f(1:3)
          ! calc potential energy
          Epot = Epot + E_each
        enddo
      enddo
    enddo
    ! self
    ix2 = ix1
    iy2 = iy1
    iz2 = iz1
    do ii = 1, cell_counter(ix1,iy1,iz1)
      iatom1 = cell_list(ix1,iy1,iz1,ii)
      itype1 = part_type(iatom1)
      do ij = 1, cell_counter(ix2,iy2,iz2)
        iatom2 = cell_list(ix2,iy2,iz2,ij)
        itype2 = part_type(iatom2)
        if( iatom1 .ge. iatom2) cycle
        ! calc parameter
        aij = aij_list(1, itype1, itype2)
        aij_half_cutoff = aij_list(2, itype1, itype2)
        ! calc force
        dr(1:3) = coord(iatom1,1:3) - coord(iatom2,1:3)
        dv(1:3) = velocity(iatom1,1:3) - velocity(iatom2,1:3)
        call calc_force_body(dr, dv, f, E_each, myid, aij, aij_half_cutoff)
        for_dpd(iatom1,1:3) = for_dpd(iatom1,1:3) + f(1:3)
        for_dpd(iatom2,1:3) = for_dpd(iatom2,1:3) - f(1:3)
        ! calc potential energy
        Epot = Epot + E_each
      enddo
    enddo
  enddo
  enddo
  enddo

  deallocate( cell_list )
end subroutine

subroutine calc_force_bond(natom, coord, force_bond, box_length, Ebond)
  use dpd_vars, only : num_bonds, bond_list
  ! global
  integer natom
  double precision coord(natom,3)
  double precision force_bond(natom, 3)
  double precision box_length(3)
  double precision Ebond
  ! local
  integer ii, now_bond, iatom1, iatom2
  double precision dr(3), dr0(3), f(3), length, length2, dx
  double precision exppart, expminusone, const_exp
  double precision half_length(3)
  
  half_length(1:3) = box_length(1:3)/2.0d0

  Ebond = 0.0
  force_bond = 0.0d0
  !$OMP parallel do &
  !$OMP& default(none) private(now_bond, ii, iatom1, iatom2, const_exp, f, &
  !$OMP& expminusone, exppart, dx, dr, dr0, length, length2) &
  !$OMP& shared(coord, bond_list, num_bonds) &
  !$OMP& firstprivate(box_length, half_length) &
  !$OMP& reduction(+:force_bond, Ebond)
  do now_bond = 1, num_bonds
    iatom1 = bond_list(now_bond)%part1
    iatom2 = bond_list(now_bond)%part2
    dr(1:3) = coord(iatom1,1:3) - coord(iatom2,1:3)
    do ii = 1, 3
      if( dr(ii) .gt. half_length(ii))then
        dr(ii) = dr(ii) - box_length(ii)
      elseif( dr(ii) .lt. -half_length(ii))then
        dr(ii) = dr(ii) + box_length(ii)
      endif
    enddo
    length2 = dot_product(dr,dr)
    length = sqrt(length2)
    dr0(1:3) = dr(1:3)/length
    dx = bond_list(now_bond)%length - length
    if( bond_list(now_bond)%bond_type .eq. 1 )then
      f(1:3) =  bond_list(now_bond)%const*dx*dr0
      Ebond = Ebond + 0.5*bond_list(now_bond)%const*dx*dx
    elseif( bond_list(now_bond)%bond_type .eq. 2)then
      !Morse bond
      !r=dr.length();
      !exppart = exp(-bconst*(r-r0));
      !expminusone = exppart - 1.0;
      !ftmp = -2*aconst*bconst*exppart*expminusone * dr/r;
      !ene = aconst*expminusone*expminusone;
      exppart = exp(-bond_list(now_bond)%const2*abs(dx))
      expminusone = exppart - 1.0
      const_exp = bond_list(now_bond)%const*expminusone
      f(1:3) =  2.0d0*const_exp*bond_list(now_bond)%const2*exppart*dr0
      !f(1:3) =  -2.0d0*const1*const2*exppart*expminusone*dr0
      Ebond = Ebond + const_exp*expminusone
    endif
    force_bond(iatom1,1:3) = force_bond(iatom1,1:3) + f(1:3)
    force_bond(iatom2,1:3) = force_bond(iatom2,1:3) - f(1:3)
  enddo
end subroutine

subroutine read_position_velocity(filename,natom,coord)
  implicit none
  ! global
  character(80) filename
  integer natom
  double precision coord(natom,3)
  ! local
  double precision dtemp(3)
  integer numter, iatom

  open(20,file=filename)

  do numter=1,natom
    read(20,*) iatom ,dtemp(1), dtemp(2), dtemp(3)
    coord(iatom,1:3) = dtemp(1:3)
  end do
  close(20)
end subroutine

subroutine random_initial_vel(natom,velocity)
  implicit none
  ! global
  integer natom
  double precision velocity(natom,3)
  ! local
  integer ii, ij
  double precision drand

  do ii = 1, natom
    do ij = 1, 3
      call genrand_wrapper(drand, 1)
      velocity(ii,ij) = drand -0.5d0
    enddo
  enddo
end subroutine

subroutine output_energy_dat(dt, nstep, natom,ene, times)
  use energy_module
  use filenames, only : ene_file
  implicit none
  ! global
  double precision dt
  integer nstep, natom, times
  type ( energy_var ) ene
  ! local

  ene%pot       = ene%pot    / dble(natom)
  ene%kin       = ene%kin    / dble(natom)
  ene%potave    = ene%potave / dble(natom)
  ene%kinave    = ene%kinave / dble(natom)

  if( times .eq. 1 )then
    open(31, file=ene_file)
    write(31,"(a82)") "#    1    2       3          4          5        6         7      8            9"
    write(31,"(a90)") "# time, total, pot_total, kin_total, total_ave, pot_ave, kin_ave, pot_nonbond, pot_bond"
  endif
  write(31,"(9f18.8)") dt*dble(nstep), ene%pot(1)+ene%kin, ene%pot(1), ene%kin, &
    ene%potave(1)+ene%kinave, ene%potave(1), ene%kinave, ene%pot(2), ene%pot(3)

end subroutine

subroutine output_checkpoint(filename, natom, coord)
  implicit none
  ! global
  character(80) filename
  integer natom
  double precision coord(natom,3)
  ! local
  integer iatom

  open(32, file= filename)
  do iatom = 1, natom
    write(32,"(i8,3f18.5)") iatom, coord(iatom,1), coord(iatom,2), coord(iatom,3)
  enddo

  close(32)
end subroutine

subroutine calc_force_body(dr, dv, f, E_each, myid, aij, aij_half_cutoff)
  use dpd_param, only : sigma_kdt, gamma, cutoff, cutoff2
  implicit none
  ! global
  double precision dr(3), dv(3), f(3), E_each
  double precision aij, aij_half_cutoff
  integer myid
  ! local
  double precision l2, l, dr0(3), ip, ll, fpair, fdpd, frand, zeta
  double precision drand1, drand2
  double precision , parameter :: pi_2 = 2.0d0*dacos(-1.0d0)

  ! kdt = sqrt(1.0/dt)
  ! sigma = sqrt(2.0*kbt*gamma)
  ! dr = r0 - r1
  ! dv = v0 - v1
  ! l = sqrt(sum(dr*dr))
  ! dr0 = dr/l
  ! ip = dr0*dv
  ! ll = 1.0 - l/cutoff
  ! wd = ll*ll
  ! wr = abs(ll)
  ! zeta = gaussian_random_number
  ! f = aij*cutoff*ll*dr/l
  ! fdpd = -gamma*wd*ip*dr0
  ! frand = sigma*wr*zeta*kdt*dr0
  ! f = f + fdpd + frand
  ! ene = 0.5*aij*cutoff*ll*ll

  f = 0.0d0
  E_each = 0.0d0
  l2 = dot_product(dr,dr)
  if (l2 .ge. cutoff2) return
  l = sqrt(l2)
  dr0 = dr/l
  ip = dot_product(dr0,dv)
  ll = 1.0d0 - l/cutoff
  call genrand_wrapper(drand1, myid)
  call genrand_wrapper(drand2, myid)
  zeta = sqrt( -2.0*log(drand1))*sin(pi_2*drand2)
  
  fpair = aij
  fdpd = -gamma*ll*ip
  frand = sigma_kdt*zeta
  f = (fpair + fdpd + frand)*ll*dr0
  E_each = aij_half_cutoff*ll*ll
  
end subroutine

subroutine read_parameter()
  use simulation_setting
  use dpd_param
  use cell_index_param
  use filenames
  implicit none
  ! local
  integer ii
  character (len=255) :: finp ! input file name

  namelist /input1/ natom, nstep_total, dump_step, dump_ene_step
  namelist /data_file/ init_pos_file, init_vel_file, aij_file, particle_file, bond_file, &
    pos_file, vel_file, xyz_file, ene_file, input_dump_file
  namelist /phys_param/ box_length, dt, kb, lambda, temperature, rand_seed
  namelist /calc_flags/ use_random_initial_coord
  namelist /dpd_parameter/ gamma, cutoff

  call get_command_argument(1, finp)
  open(99, file=finp, status='old')

  read(99, input1)

  read(99, data_file)
  
  read(99, phys_param )
  half_dt_mass = (dt/2.0d0)

  read(99, dpd_parameter)
  kdt = dsqrt(1.0d0/dt)
  sigma = dsqrt(2.0*kb*temperature*gamma)
  sigma_kdt = sigma*kdt
  cutoff2 = (cutoff)**2

  ! calc cell_index_param
  cell_size(1:3) = int(box_length(1:3)/(cutoff))
  do ii = 1, 3
    if( cell_size(ii) .le. 2)then
      write(*,*) "cell_size is too small. cell_size = ", cell_size(ii)
    endif
  enddo
  allocate( cell_counter(0:cell_size(1)+1,0:cell_size(2)+1,0:cell_size(3)+1))

  read(99,calc_flags)

  close(99)

  ! dump input data
  open(34, file=input_dump_file)
  call output_now_time(34, 'tot_start ')

  write(34,"(a30)")    "&input1"
  write(34,"(a30,i10)")"  natom = ", natom
  write(34,"(a30,i10)")"  nstep_total   = ", nstep_total
  write(34,"(a30,i10)")"  dump_step     = ", dump_step
  write(34,"(a30)")     "&end"
  write(34,"(a30)")""
  write(34,"(a30)")"&data_file"
  write(34,"(a30,a80)")" init_pos_file = " , init_pos_file
  write(34,"(a30,a80)")" init_vel_file = " , init_vel_file
  write(34,"(a30,a80)")" aij_file =      " , aij_file
  write(34,"(a30,a80)")" particle_file = " , particle_file
  write(34,"(a30,a80)")" bond_file     = " , bond_file
  write(34,"(a30,a80)")" pos_file      = " , pos_file
  write(34,"(a30,a80)")" vel_file      = " , vel_file
  write(34,"(a30,a80)")" xyz_file      = " , xyz_file
  write(34,"(a30,a80)")" ene_file      = " , ene_file
  write(34,"(a30,a80)")"input_dump_file= " , input_dump_file
  write(34,"(a30)")"&end"
  write(34,"(a30)")""
  write(34,"(a30)")"&phys_param"
  write(34,"(a30,f15.5)")"  box_length(1) = ", box_length(1)
  write(34,"(a30,f15.5)")"  box_length(2) = ", box_length(2)
  write(34,"(a30,f15.5)")"  box_length(3) = ", box_length(3)
  write(34,"(a30,f15.5)")"  dt          = ", dt
  write(34,"(a30,f15.5)")"  kb          = ", kb
  write(34,"(a30,f15.5)")"  temperature = ", temperature
  write(34,"(a30,f15.5)")"  lambda      = ", lambda
  write(34,"(a30,i15)")"  rand_seed   = ", rand_seed
  write(34,"(a30)")"&end"
  write(34,"(a30)")""
  write(34,"(a30)")"&dpd_parameter"
  write(34,"(a30,f15.5)")"  gamma         = ", gamma         
  write(34,"(a30,f15.5)")"  cutoff        = ", cutoff       
  write(34,"(a30)")"&end"
  write(34,"(a30,f15.5)")"  sigma         = ", sigma
  write(34,"(a30)")""
  write(34,"(a30)")"&calc_flags"
  write(34,"(a30,l3)")"  use_random_initial_coord   = ", use_random_initial_coord  
  write(34,"(a30)")"&end"

  write(34,"(a30,3i5)") 'cell_size = ', cell_size(1), cell_size(2), cell_size(3)

end subroutine

subroutine init_dpd_vars(natom)
  use simulation_setting, only : use_random_initial_coord, box_length, temperature
  use dpd_vars
  use filenames
  use dpd_param, only : aij_list
  implicit none
  ! global
  integer natom
  ! local
  integer part_num
  double precision now_temperature

  allocate(coord(natom,3), velocity(natom,3), force(natom,3))
  allocate(vel_dpd(natom,3), for_dpd(natom,3))
  allocate(force_bond(natom,3))
  allocate(part_type(natom), fix_flag(natom))

  coord = 0.0d0
  velocity = 0.0d0
  force = 0.0d0
  part_type = -1

  ! particle type reading
  call read_part_type(natom, part_type, fix_flag, particle_file)
  call read_part_max(aij_file, part_num)
  
  allocate(aij_list(2,part_num,part_num))
  call read_aij(aij_file, part_num, aij_list)
  ! bond reading
  call read_bond(bond_file)
  
  ! coord reading
  if( use_random_initial_coord )then
    call random_initial_coord(natom, coord, box_length)
    call random_initial_vel(natom, velocity)
    call clear_f_or_v_for_fixed(natom, velocity)
    
    ! relax structure
    call steepest_descent_method(natom, part_type, coord, box_length, 5)
  else
    call read_position_velocity(init_pos_file, natom,coord)
    call read_position_velocity(init_vel_file, natom,velocity)
  endif
  call check_bounds(natom, coord, box_length)
  
  ! init velocity
  call remove_flow( natom, velocity )
  now_temperature = ( sum(velocity(:,:)*velocity(:,:)))/( 3.0d0 * dble(sum(fix_flag)) )
  velocity(:,:) = dsqrt(temperature/now_temperature)*velocity(:,:)


end subroutine

subroutine read_bond(bond_file)
  use dpd_vars, only : bond_list, num_bonds
  implicit none
  ! global
  character(80) bond_file
  ! local
  integer bond_count, ii, part1, part2, itype
  double precision length, const, const2
  character(40) dummy
  
  bond_count = 0
  open(20, file = bond_file)
  do
    read(20,*, end=100) dummy
    bond_count = bond_count + 1
  enddo
  100 close(20)
  bond_count = bond_count - 1
  num_bonds = bond_count

  allocate(bond_list(num_bonds))

  open(20, file = bond_file)
  read(20,*) dummy
  do ii = 1, bond_count
  read(20,*) itype, part1, part2, length, const, const2
    bond_list(ii)%bond_type = itype
    bond_list(ii)%part1 = part1
    bond_list(ii)%part2 = part2
    bond_list(ii)%length = length
    bond_list(ii)%const = const
    if( bond_list(ii)%bond_type .eq. 1 )then
      bond_list(ii)%const2 = 0
    else
      bond_list(ii)%const2 = const2
    endif
  enddo
end subroutine

subroutine read_part_type(natom, part_type, fix_flag, particle_file)
  implicit none
  ! global
  integer natom, part_type(natom), fix_flag(natom)
  character(80) particle_file
  ! local
  integer ii, iatom, itemp1, itemp2
  character(40) dummy

  open(20, file = particle_file)
  read(20, *) dummy
  do ii = 1, natom
    read(20,*) iatom, itemp1, itemp2
    part_type(iatom) = itemp1
    fix_flag(iatom) = itemp2
  enddo
  close(20)
  ! check
  do ii = 1, natom
    if( part_type(ii) .eq. -1)then
      write(*,*) "not specified particle type. num = ", ii
    endif
  enddo

end subroutine

subroutine read_part_max(aij_file, part_num)
  implicit none
  ! global
  character(80) aij_file
  integer part_num
  ! local
  character(40) dummy
  integer i1, i2
  double precision dtemp
  
  part_num = 0
  open(20, file = aij_file)
  read(20,*) dummy
  do
    read(20,*, end=100) i1, i2, dtemp
    if( part_num .lt. i1 )then
      part_num = i1
    endif
    if( part_num .lt. i2 )then
      part_num = i2
    endif
  enddo
  100 close(20)
end subroutine

subroutine read_aij(aij_file, part_num, aij_list )
  use dpd_param, only : cutoff
  implicit none
  ! global
  character(80) aij_file
  integer part_num
  double precision aij_list(2, part_num, part_num)
  ! local
  character(40) dummy
  integer i1, i2, ii, ij
  double precision dtemp
  
  aij_list = -10000
  
  open(20, file = aij_file)
  read(20,*) dummy
  do
    read(20,*, end=100) i1, i2, dtemp
    aij_list(1, i1,i2) = dtemp
  enddo
  100 close(20)

  do ii = 1, part_num
    do ij = 1, part_num
      if( aij_list(1, ii,ij) .eq. -10000)then
        aij_list(1, ii,ij) = aij_list(1, ij,ii)
      endif
    enddo
  enddo
  
  do ii = 1, part_num
    do ij = 1, part_num
      if( aij_list(1, ii,ij) .eq. -10000)then
        write(*,*) 'aij parameter lack ', ii, ij
      endif
    enddo
  enddo

  ! calc aij_half_cutoff
  aij_list(2,:,:) = 0.5d0*cutoff*aij_list(1,:,:)
end subroutine

subroutine random_initial_coord(natom, coord, box_length)
  implicit none
  ! global
  integer natom
  double precision coord(natom,3), box_length(3)
  ! local
  integer iatom, ii
  double precision drand
  
  do iatom = 1, natom
    do ii = 1, 3
      call genrand_wrapper(drand,1)
      coord(iatom,ii) = box_length(ii)*(drand -0.5d0)
    enddo
  enddo
end subroutine

subroutine output_now_time(fnum,mode)
  use time_measure
  implicit none
  !global
  integer fnum
  character(10) mode
  ! output time
  character (len=10) :: stime1, stime2, stime3  ! for output date
  ! local
  integer(8) ii
  integer(8) jst(8)
  integer(4) t_end, diff, t_max, t_rate
  double precision delta_time

  if( mode .eq. 'tot_start ')then
    call date_and_time(stime1,stime2,stime3,jst)
    write(fnum,"(3i6,3x,3(i2,a1),i3)") &
      (jst(ii),ii=1,3),jst(5),":",jst(6),":", jst(7),".",jst(8)
    call system_clock(t_start_total)
  elseif( mode .eq. 'tot_end   ')then
    call system_clock(t_end, t_rate, t_max )
    if ( t_end < t_start_total ) then
      diff = t_end - t_start_total + t_max
    else
      diff = t_end - t_start_total
    endif
    write(fnum,"(a,f15.3,a)") " elapsed time :", diff/dble(t_rate), " ( seconds )"
    write(fnum,"(a,f15.3,a)") " elapsed time :", diff/(dble(t_rate)*3600), " ( hours )"
  endif

  if( mode .eq. 'step_start')then
    call system_clock(t_start_per_step)
  elseif( mode .eq. 'step_end  ')then
    call system_clock(t_end, t_rate, t_max )
      diff = t_end - t_start_per_step + t_max
    if ( t_end < t_start_per_step ) then
      diff = t_end - t_start_per_step + t_max
    else
      diff = t_end - t_start_per_step
    endif
    delta_time = (dble(diff)/dble(t_rate*n_step_for_time))
    write(fnum,"(a,f15.9,a)") " elapsed time per step :", delta_time, " ( seconds )"
  endif
end subroutine

! steepest descent method
subroutine steepest_descent_method(natom, part_type, coord, box_length, initialize_step)
  implicit none
  ! global
  integer natom, initialize_step
  double precision coord(natom,3)
  integer part_type(natom)
  double precision box_length(3)
  ! local
  double precision delta_E, dRMS, Epot
  double precision :: Epot_old = -1.0d100
  double precision :: delta_r = 1.0d0
  double precision , parameter :: delta_max = 20.0d0
  integer nstep
  double precision force_bond(natom,3), force(natom,3), Ebond
  double precision vel_dummy(natom, 3)

  vel_dummy = 0.0d0
  do nstep = 1, initialize_step
    call calc_init_bond(natom, coord, force_bond, box_length)
    call calc_force_bond(natom, coord, force_bond, box_length, Ebond)
    call calc_force_cell_index(natom, part_type, coord, vel_dummy, force, box_length, Epot)
    force = force + force_bond
    call clear_f_or_v_for_fixed(natom, force)
    Epot = Ebond + Epot

    write(*,*) "#", nstep, Epot

    delta_E = Epot - Epot_old
    dRMS = dsqrt( sum( force(:,:)*force(:,:) ))

    ! converged ?
    if( (dRMS .lt. 0.001d0 ).or.(dabs(delta_E) .lt. 1.0d0))then
      if( Epot .lt. 0.0d0)then
        exit
      endif
    endif

    if( delta_E .gt. 0.0d0)then
      delta_r = 0.5*delta_r
    else
      delta_r = min( delta_max, 1.2*delta_r)
    endif

    coord = coord + (delta_r/dRMS)*force
    call check_bounds(natom, coord, box_length)
    Epot_old = Epot
  enddo
end subroutine

subroutine calc_init_bond(natom, coord, force_bond, box_length)
  use dpd_vars, only : num_bonds, bond_list
  ! global
  integer natom
  double precision coord(natom,3)
  double precision force_bond(natom, 3)
  double precision box_length(3)
  ! local
  integer ii, now_bond, iatom1, iatom2
  double precision dr(3), dr0(3), length
  double precision half_length(3)
  
  half_length(1:3) = box_length(1:3)/2.0d0

  force_bond = 0.0d0
  do now_bond = 1, num_bonds
    iatom1 = bond_list(now_bond)%part1
    iatom2 = bond_list(now_bond)%part2
    dr(1:3) = coord(iatom2,1:3) - coord(iatom1,1:3)
    do ii = 1, 3
      if( dr(ii) .gt. half_length(ii))then
        dr(ii) = dr(ii) - box_length(ii)
      elseif( dr(ii) .lt. -half_length(ii))then
        dr(ii) = dr(ii) + box_length(ii)
      endif
    enddo
    length = sqrt(sum(dr(1:3)*dr(1:3)))
    dr0(1:3) = bond_list(now_bond)%length*(dr(1:3)/length)
    coord(iatom2,1:3) = coord(iatom1,1:3) + dr0
  enddo
end subroutine

subroutine clear_f_or_v_for_fixed(natom, velocity)
  use dpd_vars, only : fix_flag
  implicit none
  ! global
  integer natom
  double precision velocity(natom,3)
  ! local
  integer ii

  do ii = 1, natom
    velocity(ii,1:3) = velocity(ii,1:3)*fix_flag(ii)
  enddo
end subroutine
