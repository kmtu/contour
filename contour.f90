PROGRAM contour
  IMPLICIT NONE
  INTEGER :: num_binx, num_biny
  REAL(KIND=8) :: box_sizex, box_sizey
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: contour_map
  INTEGER, PARAMETER :: hist_fileid = 11
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 6
  CHARACTER(LEN=80) :: hist_filename, input_filename
  CHARACTER(LEN=128) :: command, usage, arg, output_filename
  INTEGER :: num_frames, num_atoms, dummy_int
  INTEGER :: i, j, k, n
  REAL(KIND=8) :: dummy_real
  INTEGER :: ix, iy
  REAL(KIND=8) :: rx, ry
  INTEGER :: index_x, index_y
  INTEGER, PARAMETER :: num_halothane = 16
  REAL(KIND=8) :: ref_originx, ref_originy, ref_binwidthx, ref_binwidthy
  INTEGER, DIMENSION(2, 100) :: point
  REAL(KIND=8) :: middlepointx, middlepointy
  INTEGER, PARAMETER :: num_halothane_atoms = 8
  
  call initialize()

  open(unit=hist_fileid, file=hist_filename)
  read(hist_fileid,*) num_frames, dummy_int, dummy_int, num_atoms
  write(*,*) "number of frames = ", num_frames 
  write(*,*) "number of atoms = ", num_atoms
  write(*,*) "number of halothane = ", num_halothane
  write(*,*) "number of bins in x direction = ", num_binx
  write(*,*) "number of bins in y direction = ", num_biny

  !skip 9 rows
  do i=1, 9
     read(hist_fileid,*)
  end do

  do i = 1, num_frames  
     read(hist_fileid,*) box_sizex
     read(hist_fileid,*) dummy_real, box_sizey 
     !write(*,*) box_sizex, box_sizey
     read(hist_fileid,*) !skip z dimension

     !caculate the bin width of the reference frame (1st frame)
     if (i == 1) then
        ref_originx = (box_sizex / 2) - box_sizex
        ref_originy = (box_sizey / 2) - box_sizey     
        ref_binwidthx = box_sizex / num_binx
        ref_binwidthy = box_sizey / num_biny
     end if

     !read halothanes coorordinates
     do j = 1, num_halothane
        read(hist_fileid,*) !skip first atom
        read(hist_fileid,*) ix, iy
        !write(*,*) ix, iy
        call intre2(ix, iy, rx, ry)  !change the mode from integer to real 
        rx = rx + (box_sizex / 2)
        index_x = CEILING(rx / (box_sizex / num_binx))
        ry = ry + (box_sizey / 2)
        index_y = CEILING(ry / (box_sizey / num_biny))
        !write(*,*) rx, ry
        !write(*,*) index_x, index_y
        contour_map(index_x, index_y) = contour_map(index_x, index_y) + 1
        do k = 1, num_halothane_atoms - 2
           read(hist_fileid,*) !skip other six atoms
        end do
     end do
     
     !skip other atoms 
     do j = 1, num_atoms -(num_halothane * num_halothane_atoms)
        read(hist_fileid,*)
     end do
  end do

  call output()
     
  close(hist_fileid)
  STOP
  
CONTAINS
  SUBROUTINE initialize()
    IMPLICIT NONE
    INTEGER :: error, stat
    n = COMMAND_ARGUMENT_COUNT()
    call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
    usage = "Usage: " // TRIM(ADJUSTL(command)) // " -f <input_filename> &
         &-bx <num_binx> -by <num_biny>"
    
    if (n < REQUIRED_NUM_PAR) then
       write(*,"(A,I1,A)") "At least ", REQUIRED_NUM_PAR/2, "&
            & arguments are needed."
       write(*,*) usage
       call EXIT(1)
    end if
    
    i = 1
    do while (i <= n)
       call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
       select case (arg)
       case ('-f')
          i = i + 1
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=hist_filename,&
               & STATUS=stat)
          if (stat /= 0) then
             write(*,*) "Unable to read the value of argument -f"
             write(*,*) usage
             call EXIT(1)
          end if
          if (stat/=0) then
             write(*,*) "Unable to read the value of argument -f, a&
                  & input_filename is needed !"
             call EXIT(1)
          end if
       case ('-bx')
          i = i + 1
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
          if (stat /= 0) then
             write(*,*) "Unable to read the value of argument -bx"
             write(*,*) usage
             call EXIT(1)
          end if
          read(arg,*, IOSTAT=stat) num_binx
          if (stat/=0) then
             write(*,*) "Unable to read the value of argument -bx, a&
                  & integer is needed !"
             call EXIT(1)
          end if
       case ('-by')
          i = i + 1
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
          if (stat /= 0) then
             write(*,*) "Unable to read the value of argument -by"
             write(*,*) usage
             call EXIT(1)
          end if
          read(arg,*, IOSTAT=stat) num_biny
          if (stat/=0) then
             write(*,*) "Unable to read the value of argument -by, a&
                  & integer is needed !"
             call EXIT(1)
          end if
       case default
          write(*,*) "Undefined argument: ", arg
          write(*,*) usage
          call EXIT(1)
       end select
       i = i + 1
    end do
          
    !hist_filename = "/home/pc104/HISTORY.20100504"
    !num_binx = 40
    !num_biny = 40
    ALLOCATE(contour_map(num_binx, num_biny), STAT=error)
    if (error /= 0) then
       write(*,*) "Allocate contour_map error!"
       STOP
    end if
    contour_map = 0
  END SUBROUTINE initialize

  SUBROUTINE intre2(ixx, iyy, rxx, ryy)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ixx, iyy
    REAL(KIND=8), INTENT(OUT) :: rxx, ryy
    rxx = dble(ixx)*5.0d2/dble(2**30)
    ryy = dble(iyy)*5.0d2/dble(2**30)
  END SUBROUTINE intre2
  
  SUBROUTINE output()
    IMPLICIT NONE
    INTEGER, PARAMETER :: output_id = 12
    INTEGER :: num

    do i = 1, num_binx
       write(*,"(40I5)") contour_map(:, i)
    end do
    num = SCAN(hist_filename, "/", .TRUE.)
    num = num + 1
    !    output_filename(num:) = output_filename
    !corrected by TuTu
    output_filename = TRIM(hist_filename(num:))
    
    open(unit=output_id, file=TRIM(ADJUSTL(output_filename)) // ".contour")
    write(output_id,*) "#number of frames = ", num_frames 
    write(output_id,*) "#number of atoms = ", num_atoms
    write(output_id,*) "#number of halothane = ", num_halothane
    write(output_id,*) "#number of bins in x direction = ", num_binx
    write(output_id,*) "#number of bins in y direction = ", num_biny
    write(output_id,*) "# POINT(1)                POINT(2)                       map  "
    !write the result in file, and the form of the file is which we want
    do i = 1, num_binx
       do j = 1, num_biny
          middlepointx = ref_originx + (i * ref_binwidthx)- (ref_binwidthx / 2)
          middlepointy = ref_originy + (j * ref_binwidthy)- (ref_binwidthy / 2)
          write(output_id,*) middlepointx, middlepointy, contour_map(i, j)
       end do
       write(output_id,*)
    end do
    
    call get_middlepoint() !confirm the work is ok
  END SUBROUTINE output

  SUBROUTINE get_middlepoint()
    IMPLICIT NONE
    do i = 1, num_binx
       do j = 1, num_biny
          middlepointx = ref_originx + (i * ref_binwidthx)- (ref_binwidthx / 2)
          middlepointy = ref_originy + (j * ref_binwidthy)- (ref_binwidthy / 2)
          write(*,*) middlepointx, middlepointy, contour_map(i, j)
       end do
    end do
  END SUBROUTINE get_middlepoint
          
END PROGRAM contour
