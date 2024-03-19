! Created by cf on 02/08/23.

module avg_grad_scheme
    use model, only: f, sound_speed, alfven_speed, fast_speed, calc_primitives, calc_cons_quantities
    use lxf_scheme, only: lxf_flux => flux, courant_condition
    implicit none
    private
    public solve

    contains

        subroutine solve(u, T, dx, times)
            real, target, dimension(:,:,:), intent(inout) :: u
            integer, intent(in) :: T
            real, intent(in) :: dx
            real, dimension(:), intent(inout) :: times
            real, dimension(size(u,2)) :: cs, va, vf
            real, pointer :: rho(:), Bx(:), pg(:), vx(:)
            integer :: i, k
            real :: dt, total_t

            total_t = 0.0
            times(1) = total_t

            do i = 1,T-1
                ! Determine stable timestep
                rho => u(1,:,i)
                Bx => u(5,:,i)
                pg => u(8,:,i)
                vx => u(2,:,i)

                !print *, "after pointer assignation"

                cs = sound_speed(pg=pg, rho=rho)
                va = alfven_speed(Bx=Bx, rho=rho)
                vf = fast_speed(cs=cs, va=va)

                dt = courant_condition(vx=vx/rho, vf=vf, dx=dx)
                if (isnan(dt)) stop 1
                !print *, "after courant condition"
                total_t = total_t + dt
                print *, "step: ", i, "; dt= ", dt, "; t= ", total_t

                if (size(u, 2) < 9) then
                    do k = 1, 8
                        write(*,*) u(k,:,i)
                    end do
                end if

                u(:,:,i+1) = step(u=u(:,:,i), dx=dx, dt=dt)
                times(i+1) = total_t

            end do

        end subroutine solve

        function step(u, dx, dt) result(ans)
            real, intent(in) :: u(:,:)
            real, intent(in) :: dx, dt
            real, dimension(size(u,1), size(u,2)) :: ans
            real, dimension(size(u,1), size(u,2) - 5) :: avg_flux
            integer :: cols, avg_flux_cols, i
            cols = size(u,2)
            avg_flux_cols = size(avg_flux, 2)

            avg_flux = aflux(u=u, dx=dx, dt=dt)

            ans(:,4:cols-3) = u(:,4:cols-3) + (dt/dx)*(avg_flux(:,1:avg_flux_cols-1) - avg_flux(:,2:avg_flux_cols))

            ! Ghost cells
            do i = 1,3
                ans(:,i) = ans(:,4)
            end do

            do i = cols-2,cols
                ans(:,i) = ans(:,cols-3)
            end do

        end function step

        function aflux(u, dx, dt) result(avg_flux)
            real, intent(in) :: u(:,:)
            real, intent(in) :: dx, dt
            real, dimension(size(u,1), size(u,2) - 5) :: pl, pr, ul, ur, fl, fr, avg_flux

            call avg_gradient(u=u, dx=dx, dt=dt, pl=pl, pr=pr)
            call calc_cons_quantities(pl, ul)
            call calc_cons_quantities(pr, ur)
            call f(ur, fr)
            call f(ul, fl)
            avg_flux = 0.5*(fr + fl) - (dt/dx)*(ur - ul)

        end function aflux

        subroutine avg_gradient(u, dx, dt, pl, pr)
            real, intent(in) :: u(:,:)
            real, intent(in) :: dx, dt
            real, dimension(size(u,1), size(u,2) - 5), intent(out) :: pl, pr
            real, dimension(size(u,1), size(u,2) - 4) :: p, avg_grad !, diff_p
            integer :: cols, flux_cols, p_half_cols, diff_p_cols, p_cols
            real, dimension(size(u,1), size(u,2) - 3) :: diff_p
            real, dimension(size(u,1), size(u,2) - 2) :: u_half, p_half
            real, dimension(size(u,1), size(u,2) - 1) :: lxf_flux_val

            cols = size(u,2)
            diff_p_cols = size(diff_p, 2)
            flux_cols = size(lxf_flux_val, 2)
            p_half_cols = size(p_half, 2)
            p_cols = size(p, 2)

            !print *, "shape(u_half)", shape(u_half)

            lxf_flux_val = lxf_flux(u, dx=dx, dt=0.5*dt)

            u_half = u(:,2:cols-1) + 0.5*(dt/dx)*(lxf_flux_val(:,1:flux_cols-1) - lxf_flux_val(:,2:flux_cols))
            call calc_primitives(u_half, p_half)

            diff_p = p_half(:,1:p_half_cols-1) - p_half(:,2:p_half_cols)

            avg_grad = avg(diff_p(:,2:diff_p_cols), diff_p(:,1:diff_p_cols-1))/dx

            p = p_half(:,2:p_half_cols-1)

            pl = p(:,1:p_cols-1) + 0.5*dx*avg_grad(:,1:p_cols-1)
            pr = P(:,2:p_cols) - 0.5*dx*avg_grad(:,2:p_cols)

        end subroutine avg_gradient

        function avg(a, b) result(ans)
            real, dimension(:,:), intent(in) :: a, b
            real, dimension(size(a,1), size(a,2)) :: ans
            integer :: cols, rows, i, j
            rows = size(a, 1)
            cols = size(a, 2)

            do concurrent(i = 1:rows)
                do j = 1, cols
                    ans(i,j) = scalar_avg(a(i, j), b(i, j))
                end do
            end do

        end function avg

        pure function scalar_avg(a, b) result(ans)
            real, intent(in) :: a, b
            real :: ans
            if (a**2 + b**2 /= 0 .and. a*b > 0) then
                ans = (b*a**2 + a*b**2)/(a**2 + b**2)
            else
                ans = 0.0
            end if
        end function scalar_avg

end module avg_grad_scheme