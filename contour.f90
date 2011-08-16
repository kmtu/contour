PROGRAM contour
  IMPLICIT NONE
  INTEGER :: num_binx, num_biny
  REAL(KIND=8) :: box_sizex, box_sizey
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: contour_map
  INTEGER, PARAMETER :: REQUIRED_NUM_PAR = 8
  INTEGER, DIMENSION(:), ALLOCATABLE :: hist_fileid
  CHARACTER(LEN=128), DIMENSION(:), ALLOCATABLE :: hist_filename
  CHARACTER(LEN=128) :: command, usage, arg, output_filename
  INTEGER :: num_frames, num_atoms, dummy_int, num_hist_files
  INTEGER :: i, j, k, n, m, mark
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

  do m = 1, num_hist_files
     read(hist_fileid(m),*) num_frames, dummy_int, dummy_int, num_atoms
     write(*,*) "number of frames = ", num_frames 
     write(*,*) "number of atoms = ", num_atoms
     write(*,*) "number of halothane = ", num_halothane
     write(*,*) "number of bins in x direction = ", num_binx
     write(*,*) "number of bins in y direction = ", num_biny
     
     !skip 9 rows
     do i=1, 9
        read(hist_fileid(m),*)
     end do

     !read frame by frame
     mark = num_frames/100
     do i = 1, num_frames  
        !output current progress
        if (MOD(i, mark) == 0) then
           write(UNIT=*, FMT="(1X,'Reading %',I3,' of file ',I3)") &
                &i/mark, m
        end if

        read(hist_fileid(m),*) box_sizex
        read(hist_fileid(m),*) dummy_real, box_sizey 
        !write(*,*) box_sizex, box_sizey
        read(hist_fileid(m),*) !skip z dimension

        !caculate the bin width of the reference frame (only once, using 1st frame)
        if (m == 1 .and. i == 1) then
           ref_originx = (box_sizex / 2) - box_sizex
           ref_originy = (box_sizey / 2) - box_sizey     
           ref_binwidthx = box_sizex / num_binx
           ref_binwidthy = box_sizey / num_biny
        end if

        !read halothanes coorordinates
        do j = 1, num_halothane
           read(hist_fileid(m),*) !skip first atom
           read(hist_fileid(m),*) ix, iy
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
              read(hist_fileid(m),*) !skip other six atoms
           end do
        end do
     
        !skip other atoms 
        do j = 1, num_atoms -(num_halothane * num_halothane_atoms)
           read(hist_fileid(m),*)
        end do
     end do
  end do
  call output()
     
!  close(hist_fileid)
  STOP
  
CONTAINS
  SUBROUTINE initialize()
    IMPLICIT NONE
    INTEGER :: error, stat
    n = COMMAND_ARGUMENT_COUNT()
    call GET_COMMAND_ARGUMENT(NUMBER=0, VALUE=command)
    usage = "Usage: " // TRIM(ADJUSTL(command)) // " -f <input_filename(s)>... &
         &-bx <num_binx> -by <num_biny> -o <output_filename>"
    
    if (n < REQUIRED_NUM_PAR) then
       write(*,"(A,I1,A)") "At least ", REQUIRED_NUM_PAR/2, "&
            & arguments are needed."
       write(*,*) usage
       call EXIT(1)
    end if
    
    i = 1
    do while (i <= n)
       call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
       i = i + 1       
       select case (arg)
       case ('-f')
          !count number of data files 
          num_hist_files = 0
          j = i
          do while (j <= n)
             call GET_COMMAND_ARGUMENT(NUMBER=j, VALUE=arg, STATUS=stat)
             j = j + 1
             if (stat /= 0) then
                write(*,*) "Error: unable to count the number of arguments -f &
                     &<data file1> [<data file2> ...]"
                write(*,*) usage
                call EXIT(1)
             else if (arg(1:1) == '-') then !end of data file arguments
                EXIT
             end if
             num_hist_files = num_hist_files + 1
          end do
          if (num_hist_files == 0) then
             write(*,*) "Error: at least one data file must be provided!"
             write(*,*) usage
             call EXIT(1)
          end if

          ALLOCATE(hist_fileid(num_hist_files))
          ALLOCATE(hist_filename(num_hist_files))

          do j = 1, num_hist_files
             call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=hist_filename(j), STATUS=stat)
             i = i + 1
             if (stat /= 0) then
                write(*,*) "Error: unable to read the value of argument -f &
                     &<data file1> [<data file2> ...]"
                write(*,*) usage
                call EXIT(1)
             end if
          end do

       case ('-bx')
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
          i = i + 1
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
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=arg, STATUS=stat)
          i = i + 1
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
          
       case ('-o')
          call GET_COMMAND_ARGUMENT(NUMBER=i, VALUE=output_filename,&
               & STATUS=stat)
          i = i + 1
          if (stat /= 0) then
             write(*,*) "Unable to read the value of argument -o"
             write(*,*) usage
             call EXIT(1)
          end if
          
       case default
          write(*,*) "Undefined argument: ", arg
          write(*,*) usage
          call EXIT(1)
       end select
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
    
    !open every data file
    do i = 1, num_hist_files
       hist_fileid(i) = 100 + i
       open(UNIT=hist_fileid(i), FILE=hist_filename(i), IOSTAT=stat, &
            &STATUS="OLD", ACTION="READ")
       if (stat /=0) then
          write(*,*) "Error: unable to open file: ", TRIM(ADJUSTL(hist_filename(i)))
          call EXIT(1)
       end if
       write(*,*) "data file ",i ,":", TRIM(ADJUSTL(hist_filename(i)))
    end do
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

!    do i = 1, num_binx
!       write(*,"(40I5)") contour_map(:, i)
!    end do

!    num = SCAN(hist_filename, "/", .TRUE.)
!    num = num + 1
    !    output_filename(num:) = output_filename
    !corrected by TuTu
!   output_filename = TRIM(hist_filename(num:))
    
!    open(unit=output_id, file=TRIM(ADJUSTL(output_filename)) // ".contour")
    open(unit=output_id, file=TRIM(ADJUSTL(output_filename)))    
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
!          write(*,*) middlepointx, middlepointy, contour_map(i, j)
       end do
    end do
  END SUBROUTINE get_middlepoint
END PROGRAM contour
