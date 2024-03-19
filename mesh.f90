! Created by cf on 31/07/23.

module mesh
    use bookkeeping, only: read_arguments
    use model, only: calc_primitives, calc_cons_quantities, init_cons_quantities, gamma

    implicit none
    private
    public T, M, D, dx, u, &
            initialize_mesh

    integer :: T ! Number of time steps
    integer :: M ! Number of cells for the range [xl, xr]
    integer, parameter :: D = 8 ! Length of the individual vectors u

    real :: dx, xl, xr
    namelist /mesh_parameters_nml/ T, M, xl, xr, gamma

    real :: pl0(D), pr0(D)
    namelist /initial_conditions_nml/ pl0, pr0

    real, dimension(:,:,:), allocatable :: u

    private xl, xr, mesh_parameters_nml
    public read_mesh_parameters

    contains

        subroutine read_mesh_parameters(mesh_configuration_file)
            character(len=100), intent(in) :: mesh_configuration_file
            integer :: fileunit

            open(newunit=fileunit, file=mesh_configuration_file, action="read")
            read(fileunit, nml=mesh_parameters_nml)
            print mesh_parameters_nml

            1 close(fileunit)

        end subroutine read_mesh_parameters

        subroutine read_initial_conditions(initial_conditions_file)
            character(len=100), intent(in) :: initial_conditions_file
            integer :: fileunit, status

            open(newunit=fileunit, file=initial_conditions_file, action="read")
            read(fileunit, nml=initial_conditions_nml)
            print initial_conditions_nml

            1 close(fileunit)
        end subroutine read_initial_conditions

        subroutine initialize_mesh(nr_ghost_cells)
            integer, intent(in) :: nr_ghost_cells
            character(len=100):: mesh_configuration_file, initial_conditions_file

            call read_arguments(mesh_configuration_file, initial_conditions_file)
            call read_mesh_parameters(mesh_configuration_file)
            call read_initial_conditions(initial_conditions_file)

            call create_mesh(pl0, pr0, nr_ghost_cells)

        end subroutine initialize_mesh

        subroutine create_mesh(pl0, pr0, nr_ghost_cells)
            integer, intent(in) :: nr_ghost_cells
            integer :: rows, cols, i
            real, intent(in) :: pl0(:), pr0(:)

            !if (.not. present(offset)) offset = 2

            dx = (xr - xl)/M

            allocate(u(D, M+nr_ghost_cells, T))
            call init_cons_quantities(pl0,u(:,:(M+nr_ghost_cells)/2,1))
            call init_cons_quantities(pr0,u(:,(M+nr_ghost_cells)/2+1:,1))

            cols = size(u,2)
            ! only print if the matrix is small
            if (cols < 9) then
                print *, "Initial values of the conserved quantities"
                rows = size(u,1)
                do i = 1, rows
                    write(*,*) u(i,:, 1)
                end do
            end if
        end subroutine create_mesh

end module mesh

