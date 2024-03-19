module scheme
    use model, only: gamma, calc_primitives, calc_cons_quantities, &
                    alfven_speed, fast_speed, slow_speed, sound_speed, f
    implicit none
    private
    public solve, print_matrix


    interface
    module function flux(u, dx, dt) result(ans)
            real, dimension(:,:), intent(in) :: u
            real, intent(in):: dx, dt
            real, dimension(size(u,1), size(u,2)-1) :: ans

    end function flux
    end interface

    contains

        subroutine solve(u, T, dx, times, whole_b)
            real, target, dimension(:,:,:), intent(inout) :: u
            integer, intent(in) :: T, whole_b
            real, intent(in) :: dx
            real, dimension(:), intent(inout) :: times
            real, dimension(size(u,2)) :: a, va, vf
            real, pointer :: rho(:), Bx(:), pg(:), vx(:)
            real, pointer, dimension(:,:) :: b
            integer :: i, j, k
            real :: dt, total_t

            total_t = 0.0
            times(1) = total_t

            do i = 1,T-1
                ! Determine stable timestep
                rho => u(1,:,i)
                Bx => u(5,:,i)
                b => u(5:7,:,i)
                pg => u(8,:,i)
                vx => u(2,:,i)

                !print *, "after pointer assignation"

                a = sound_speed(pg=pg, rho=rho)
                va = alfven_speed(Bx=Bx, rho=rho)
                vf = fast_speed(a=a, va=va, b=b, rho=rho, whole_b=whole_b)

                dt = courant_condition(vx=vx/rho, vf=vf, dx=dx)
                !print *, "after courant condition"
                total_t = total_t + dt

                print *, "============================================================================================"
                print *, "step: ", i, "; dt= ", dt, "; t= ", total_t
                call print_matrix(u(:,:,i), "u")
                
                u(:,:,i+1) = step(u=u(:,:,i), dx=dx, dt=dt)
                times(i+1) = total_t

            end do

        end subroutine solve

        function courant_condition(vx, vf, dx) result(ans)
            real, dimension(:), intent(in) :: vx, vf
            real, intent(in) :: dx
            real :: ans

            !print *, "dx=", dx
            !print *, "vx\n", vx
            !print *, "vf\n", vf

            ans = 0.3*dx/maxval(abs(vx) + vf)

        end function courant_condition

        function step(u, dx, dt) result(ans)
            real, dimension(:,:), intent(in) :: u
            real, intent(in) :: dx, dt
            real, dimension(size(u, 1), size(u, 2)) :: ans
            real, dimension(size(u, 1), size(u, 2)-1) :: flux_val
            integer :: flux_cols, u_cols

            flux_cols = size(flux_val, 2)
            u_cols = size(u, 2)

            flux_val = flux(u, dx=dx, dt=dt)
            !print *, "after flux calculation"
            call print_matrix(flux_val, "flux")

            ! ans(:,2:u_cols-1) = u(:,2:u_cols-1) + (dt/dx)*(flux_val(:,1:flux_cols-1) - flux_val(:,2:flux_cols))
            ans(:,2:u_cols-1) = u(:,2:u_cols-1) + (dt/dx)*(flux_val(:,1:flux_cols-1) - flux_val(:,2:flux_cols))
            ! Values for the ghost cells
            ans(:,1) = ans(:,2)
            ans(:,u_cols) = ans(:,u_cols-1)

        end function step

        subroutine print_matrix(m, matrix_name)
            real, dimension(:,:), intent(in) :: m
            character(*), intent(in) :: matrix_name
            integer :: r, c, i, j

            r = size(m,1)
            c = size(m,2)

            if (c < 9) then
                print *, matrix_name
                do i = 1,r
                    write(*,*) m(i,:)
                end do
            end if
        end subroutine print_matrix
            
end module scheme 
