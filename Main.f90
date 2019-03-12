!==============================================================================!
  program Buga
!------------------------------------------------------------------------------!
  implicit none
!----------------------------------[Include]-----------------------------------!
  include 'mpif.h'
!---------------------------------[New Types]----------------------------------!
  type Particle_Type
    integer          :: cell
    double precision :: dval(3)
    double precision :: x, y, z
  end type
!------------------------------------------------------------------------------!
  integer, parameter             :: n_vars = 5     ! number of variables in type
  type(Particle_Type)            :: particle       ! struct to be exchanged
  integer                        :: blens(n_vars)  ! array of blocklengths
  integer                        :: types(n_vars)  ! array of types
  integer(kind=MPI_ADDRESS_KIND) :: disp(n_vars)   ! array of displacements
  integer                        :: status(MPI_STATUS_SIZE)
  integer                        :: des, len, sou, tag, i
  integer                        :: this_proc, n_proc, new_type, error
!==============================================================================!

  !---------------!
  !   Start MPI   !
  !---------------!
  call Mpi_Init(error)

  !-------------------------------------------------!
  !   Get number of processors and processor I.D.   !
  !-------------------------------------------------!
  call Mpi_Comm_Size(MPI_COMM_WORLD, n_proc, error) 
  call Mpi_Comm_Rank(MPI_COMM_WORLD, this_proc, error)
  this_proc = this_proc + 1

  !---------------------!
  !   Create new type   !
  !---------------------!

  ! Create array of blocklengths
  blens = (/1,  &  ! cell
            3,  &  ! dval(3)
            1,  &  ! x
            1,  &  ! y
            1/)    ! z

  ! Create array of types
  types = (/MPI_INTEGER,  &  ! cell
            MPI_DOUBLE,   &  ! dval(3)
            MPI_DOUBLE,   &  ! x
            MPI_DOUBLE,   &  ! y
            MPI_DOUBLE/)     ! z

  ! Create array of addresses, first take absolute addresses ...
  call Mpi_Get_Address(particle % cell, disp(1), error)
  call Mpi_Get_Address(particle % dval, disp(2), error)
  call Mpi_Get_Address(particle % x,    disp(3), error)
  call Mpi_Get_Address(particle % y,    disp(4), error)
  call Mpi_Get_Address(particle % z,    disp(5), error)

  ! ... then make them relative
  do i = n_vars, 1, -1
    disp(i) = disp(i) - disp(1)
  end do

  call Mpi_Type_Create_Struct(n_vars,    &
                              blens,     &
                              disp,      &
                              types,     &
                              new_type,  &
                              error)
  call Mpi_Type_Commit(new_type, error)

  !----------------------!
  !   Create some data   !
  !----------------------!
  particle % cell    = 0
  particle % dval(1) = 0.0
  particle % dval(2) = 0.0
  particle % dval(3) = 0.0

  if(this_proc .eq. 1) then
    particle % cell    = 128
    particle % dval(1) =  11.0
    particle % dval(2) =  22.0
    particle % dval(3) =  33.0
  end if

  !------------------------!
  !   Exchange some data   !
  !------------------------!
  if(this_proc .eq. 1) then
    tag = 7
    des = 1
    len = 1
    call Mpi_Send(particle,                         &
                  len,                              & ! length
                  new_type,                         & ! datatype
                  des,                              & ! dest,
                  tag,                              & ! sendtag,
                  MPI_COMM_WORLD,                   &
                  error)
  else if(this_proc .eq. 2) then
    tag = 7
    sou = 0
    len = 1
    call Mpi_Recv(particle,                         &
                  len,                              & ! length
                  new_type,                         & ! datatype
                  sou,                              & ! source,
                  tag,                              & ! recvtag,
                  MPI_COMM_WORLD,                   &
                  status,                           &
                  error)
  end if

  if(this_proc .eq. 2) then
    print *, particle % cell
    print *, particle % dval(1)
    print *, particle % dval(2)
    print *, particle % dval(3)
  end if

  !-------------!
  !   End MPI   !
  !-------------!
  call Mpi_Finalize(error)

  end program
