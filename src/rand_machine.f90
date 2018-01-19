module rand1
  include 'mt95.f90'
end module rand1

module rand2
  include 'mt95.f90'
end module rand2

module rand3
  include 'mt95.f90'
end module rand3

module rand4
  include 'mt95.f90'
end module rand4

module rand5
  include 'mt95.f90'
end module rand5

module rand6
  include 'mt95.f90'
end module rand6

module rand7
  include 'mt95.f90'
end module rand7

module rand8
  include 'mt95.f90'
end module rand8

module rand9
  include 'mt95.f90'
end module rand9

module rand10
  include 'mt95.f90'
end module rand10

module rand11
  include 'mt95.f90'
end module rand11

module rand12
  include 'mt95.f90'
end module rand12

module rand13
  include 'mt95.f90'
end module rand13

module rand14
  include 'mt95.f90'
end module rand14

module rand15
  include 'mt95.f90'
end module rand15

module rand16
  include 'mt95.f90'
end module rand16

module rand17
  include 'mt95.f90'
end module rand17

module rand18
  include 'mt95.f90'
end module rand18

module rand19
  include 'mt95.f90'
end module rand19

module rand20
  include 'mt95.f90'
end module rand20

module rand21
  include 'mt95.f90'
end module rand21

module rand22
  include 'mt95.f90'
end module rand22

module rand23
  include 'mt95.f90'
end module rand23

module rand24
  include 'mt95.f90'
end module rand24

module rand25
  include 'mt95.f90'
end module rand25

module rand26
  include 'mt95.f90'
end module rand26

module rand27
  include 'mt95.f90'
end module rand27

module rand28
  include 'mt95.f90'
end module rand28

module rand29
  include 'mt95.f90'
end module rand29

module rand30
  include 'mt95.f90'
end module rand30

module rand31
  include 'mt95.f90'
end module rand31

module rand32
  include 'mt95.f90'
end module rand32

subroutine rand1_init_wrapper(seed)
  use rand1
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand2_init_wrapper(seed)
  use rand2
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand3_init_wrapper(seed)
  use rand3
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand4_init_wrapper(seed)
  use rand4
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand5_init_wrapper(seed)
  use rand5
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand6_init_wrapper(seed)
  use rand6
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand7_init_wrapper(seed)
  use rand7
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand8_init_wrapper(seed)
  use rand8
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand9_init_wrapper(seed)
  use rand9
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand10_init_wrapper(seed)
  use rand10
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand11_init_wrapper(seed)
  use rand11
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand12_init_wrapper(seed)
  use rand12
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand13_init_wrapper(seed)
  use rand13
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand14_init_wrapper(seed)
  use rand14
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand15_init_wrapper(seed)
  use rand15
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand16_init_wrapper(seed)
  use rand16
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand17_init_wrapper(seed)
  use rand17
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand18_init_wrapper(seed)
  use rand18
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand19_init_wrapper(seed)
  use rand19
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand20_init_wrapper(seed)
  use rand20
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand21_init_wrapper(seed)
  use rand21
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand22_init_wrapper(seed)
  use rand22
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand23_init_wrapper(seed)
  use rand23
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand24_init_wrapper(seed)
  use rand24
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand25_init_wrapper(seed)
  use rand25
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand26_init_wrapper(seed)
  use rand26
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand27_init_wrapper(seed)
  use rand27
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand28_init_wrapper(seed)
  use rand28
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand29_init_wrapper(seed)
  use rand29
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand30_init_wrapper(seed)
  use rand30
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand31_init_wrapper(seed)
  use rand31
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand32_init_wrapper(seed)
  use rand32
  integer seed
  call genrand_init( put=seed )
end subroutine

subroutine rand1_wrapper(grnd)
  use rand1
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand2_wrapper(grnd)
  use rand2
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand3_wrapper(grnd)
  use rand3
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand4_wrapper(grnd)
  use rand4
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand5_wrapper(grnd)
  use rand5
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand6_wrapper(grnd)
  use rand6
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand7_wrapper(grnd)
  use rand7
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand8_wrapper(grnd)
  use rand8
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand9_wrapper(grnd)
  use rand9
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand10_wrapper(grnd)
  use rand10
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand11_wrapper(grnd)
  use rand11
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand12_wrapper(grnd)
  use rand12
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand13_wrapper(grnd)
  use rand13
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand14_wrapper(grnd)
  use rand14
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand15_wrapper(grnd)
  use rand15
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand16_wrapper(grnd)
  use rand16
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand17_wrapper(grnd)
  use rand17
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand18_wrapper(grnd)
  use rand18
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand19_wrapper(grnd)
  use rand19
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand20_wrapper(grnd)
  use rand20
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand21_wrapper(grnd)
  use rand21
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand22_wrapper(grnd)
  use rand22
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand23_wrapper(grnd)
  use rand23
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand24_wrapper(grnd)
  use rand24
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand25_wrapper(grnd)
  use rand25
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand26_wrapper(grnd)
  use rand26
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand27_wrapper(grnd)
  use rand27
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand28_wrapper(grnd)
  use rand28
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand29_wrapper(grnd)
  use rand29
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand30_wrapper(grnd)
  use rand30
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand31_wrapper(grnd)
  use rand31
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine rand32_wrapper(grnd)
  use rand32
  double precision grnd
  call genrand_real3(grnd)
end subroutine

subroutine genrand_init_wrapper(rep,seed)
  implicit none
  integer rep
  integer seed
  double precision drand
  integer ii
 
  select case (rep)
    case(1)
      call rand1_init_wrapper(seed) 
    case(2)
      call rand2_init_wrapper(seed) 
    case(3)
      call rand3_init_wrapper(seed) 
    case(4)
      call rand4_init_wrapper(seed) 
    case(5)
      call rand5_init_wrapper(seed) 
    case(6)
      call rand6_init_wrapper(seed) 
    case(7)
      call rand7_init_wrapper(seed) 
    case(8)
      call rand8_init_wrapper(seed) 
    case(9)
      call rand9_init_wrapper(seed) 
    case(10)
      call rand10_init_wrapper(seed) 
    case(11)
      call rand11_init_wrapper(seed) 
    case(12)
      call rand12_init_wrapper(seed) 
    case(13)
      call rand13_init_wrapper(seed) 
    case(14)
      call rand14_init_wrapper(seed) 
    case(15)
      call rand15_init_wrapper(seed) 
    case(16)
      call rand16_init_wrapper(seed) 
    case(17)
      call rand17_init_wrapper(seed) 
    case(18)
      call rand18_init_wrapper(seed) 
    case(19)
      call rand19_init_wrapper(seed) 
    case(20)
      call rand20_init_wrapper(seed) 
    case(21)
      call rand21_init_wrapper(seed) 
    case(22)
      call rand22_init_wrapper(seed) 
    case(23)
      call rand23_init_wrapper(seed) 
    case(24)
      call rand24_init_wrapper(seed) 
    case(25)
      call rand25_init_wrapper(seed) 
    case(26)
      call rand26_init_wrapper(seed) 
    case(27)
      call rand27_init_wrapper(seed) 
    case(28)
      call rand28_init_wrapper(seed) 
    case(29)
      call rand29_init_wrapper(seed) 
    case(30)
      call rand30_init_wrapper(seed) 
    case(31)
      call rand31_init_wrapper(seed) 
    case(32)
      call rand32_init_wrapper(seed) 
  end select
  do ii = 1, 1000
    call genrand_wrapper(drand,rep)
  enddo
end subroutine

subroutine genrand_wrapper(grnd, rep)
  implicit none
  integer rep
  double precision grnd

  select case ( rep )
    case(1)
      call rand1_wrapper(grnd)
    case(2)
      call rand2_wrapper(grnd)
    case(3)
      call rand3_wrapper(grnd)
    case(4)
      call rand4_wrapper(grnd)
    case(5)
      call rand5_wrapper(grnd)
    case(6)
      call rand6_wrapper(grnd)
    case(7)
      call rand7_wrapper(grnd)
    case(8)
      call rand8_wrapper(grnd)
    case(9)
      call rand9_wrapper(grnd)
    case(10)
      call rand10_wrapper(grnd)
    case(11)
      call rand11_wrapper(grnd)
    case(12)
      call rand12_wrapper(grnd)
    case(13)
      call rand13_wrapper(grnd)
    case(14)
      call rand14_wrapper(grnd)
    case(15)
      call rand15_wrapper(grnd)
    case(16)
      call rand16_wrapper(grnd)
    case(17)
      call rand17_wrapper(grnd)
    case(18)
      call rand18_wrapper(grnd)
    case(19)
      call rand19_wrapper(grnd)
    case(20)
      call rand20_wrapper(grnd)
    case(21)
      call rand21_wrapper(grnd)
    case(22)
      call rand22_wrapper(grnd)
    case(23)
      call rand23_wrapper(grnd)
    case(24)
      call rand24_wrapper(grnd)
    case(25)
      call rand25_wrapper(grnd)
    case(26)
      call rand26_wrapper(grnd)
    case(27)
      call rand27_wrapper(grnd)
    case(28)
      call rand28_wrapper(grnd)
    case(29)
      call rand29_wrapper(grnd)
    case(30)
      call rand30_wrapper(grnd)
    case(31)
      call rand31_wrapper(grnd)
    case(32)
      call rand32_wrapper(grnd)
  end select
end subroutine

