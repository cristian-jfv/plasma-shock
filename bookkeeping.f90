! Created by cf on 01/08/23.

module bookkeeping
    implicit none
    public read_arguments, save_matrix, save_vector

contains

    subroutine read_arguments(mesh_configuration_file, initial_conditions_file)
        character(len=100), intent(out) :: mesh_configuration_file, initial_conditions_file

        if (command_argument_count() < 2 ) then
            write(*,*) "Usage ./lxf_kt <mesh configuration file> <initial conditions file>"
            stop 1
        end if

        call get_command_argument(1, mesh_configuration_file)
        call get_command_argument(2, initial_conditions_file)

        print *, "Mesh config. file: ", mesh_configuration_file
        print *, "Initial cond. file: ", initial_conditions_file
    end subroutine read_arguments

    subroutine save_matrix(m, filename)
        real, dimension(:,:), intent(in) :: m
        character(*), intent(in) :: filename
        integer :: fileunit, cols, rows, i, j
        rows = size(m, 1)
        cols = size(m, 2)

        print *, "Matrix shape: rows=", rows, ", cols=", cols

        open(newunit=fileunit, file=filename, action="write", status="replace")
        do i = 1,rows
            write(fileunit, *)(m(i,j), j=1,cols)
        end do

    end subroutine save_matrix

    subroutine save_vector(m, filename)
        real, dimension(:), intent(in) :: m
        character(*), intent(in) :: filename
        integer :: fileunit, rows, i
        rows = size(m, 1)

        print *, "Vector size: rows=", rows

        open(newunit=fileunit, file=filename, action="write", status="replace")
        do i = 1,rows
            write(fileunit, *)(m(i))
        end do

    end subroutine save_vector
end module bookkeeping