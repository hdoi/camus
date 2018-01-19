module simulation_setting
  integer natom ! num of atoms
  integer nstep_total     ! num of steps
  integer dump_step       ! dump step for xyz structure
  integer dump_ene_step       ! dump step for energy
  
  double precision box_length(3)
  double precision  dt, lambda
  double precision  half_dt_mass
  double precision  kb
  double precision  temperature
  integer(4) rand_seed

  logical use_random_initial_coord

end module simulation_setting

module cell_index_param
  integer cell_size(3)
  integer , allocatable :: cell_counter(:,:,:)
  integer , allocatable :: cell_list(:,:,:,:)
end module cell_index_param

module dpd_param
  double precision gamma, sigma
  double precision kdt, sigma_kdt
  double precision cutoff, cutoff2
  double precision , allocatable :: aij_list(:,:,:) ! (1:aij 2:aij_half_cutoff, part type, part type )
end module dpd_param

module energy_module
  type energy_var
    double precision pot(3) ! potential (1:total, 2:nonbond, 3:bond)
    double precision kin ! kinetic
    double precision potave(3)
    double precision kinave
  end type energy_var
end module energy_module

module dpd_vars
  double precision , allocatable :: coord(:,:), velocity(:,:), force(:,:)
  double precision , allocatable :: vel_dpd(:,:), for_dpd(:,:)
  double precision , allocatable :: force_bond(:,:)
  integer , allocatable :: part_type(:), fix_flag(:)
  
  integer num_bonds
  type bond_param
    integer part1
    integer part2
    integer bond_type
    double precision length
    double precision const
    double precision const2
  end type bond_param
  type (bond_param)  , allocatable :: bond_list(:)
end module dpd_vars

module filenames
  character(80) init_pos_file, init_vel_file
  character(80) aij_file, particle_file, bond_file
  character(80) pos_file, vel_file, xyz_file, ene_file, input_dump_file
end module filenames

module time_measure
  integer t_start_total
  integer t_start_per_step
  integer n_step_for_time
end module time_measure

program main
  use simulation_setting
  use energy_module
  use dpd_vars
  use filenames
  use time_measure
  implicit none
  ! global
  type ( energy_var ) ene
  ! local
  integer nstep
  integer(4) ii
  
  ! initialize rand module
  do ii = 1, 32
    call genrand_init_wrapper(ii,rand_seed + ii)
  enddo

  ! initialize
  call read_parameter()
  call init_dpd_vars(natom)
  ene%pot(1:3) = 0.0d0
  ene%kin = 0.0d0
  ene%potave(1:3) = 0.0d0
  ene%kinave = 0.0d0


  call output_xyz(natom, part_type, coord, box_length, ene)
  
  call output_now_time(34, 'step_start' )
  n_step_for_time = nstep_total

  do nstep = 1, nstep_total
    ! velocity verlet method start
    ! 1, update coord
    coord = coord + dt*(velocity + half_dt_mass*force)

    ! 2, cal vel_dpd
    vel_dpd = velocity + lambda*dt*force

    ! 3, cal force
    call calc_force_cell_index(natom, part_type, coord, vel_dpd, for_dpd, box_length, ene%pot(2))
    call calc_force_bond(natom, coord, force_bond, box_length, ene%pot(3))
    ene%pot(1) = ene%pot(2) + ene%pot(3) ! total potential = nonbond + bond
    for_dpd = for_dpd + force_bond
    call clear_f_or_v_for_fixed(natom, for_dpd)
    ene%potave(1:3) = ene%potave(1:3) + ene%pot(1:3)

    ! 4, update vel
    velocity = velocity + half_dt_mass*(force + for_dpd)
    call clear_f_or_v_for_fixed(natom, velocity)
    ene%kin = 0.5d0*sum(velocity(:,:)*velocity(:,:))
    ene%kinave = ene%kinave + ene%kin
    
    force = for_dpd
    
    ! ***** analysis start *****
    if( mod( nstep, dump_ene_step ) .eq. 0 )then
      write(*,*) '# nstep ', nstep
      ene%potave(1:3) = ene%potave(1:3) / dble(dump_ene_step)
      ene%kinave = ene%kinave / dble(dump_ene_step)
      call output_energy_dat(dt, nstep, natom, ene, nstep / dump_ene_step)
    endif

    if( mod(nstep, dump_step).eq. 0)then
      call output_xyz(natom, part_type, coord, box_length, ene)
      call output_checkpoint(pos_file, natom, coord)
      call output_checkpoint(vel_file, natom, velocity)
    endif
  enddo

  call output_now_time(34, 'step_end  ')
  call output_now_time(34, 'tot_end   ' )
  close(34)
end program

