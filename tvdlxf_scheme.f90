submodule (scheme) tvd_lxf_scheme
    use model, only: f, sound_speed, alfven_speed, fast_speed
    implicit none

    contains

        module function flux(u, dx, dt) result(ans)
            real, dimension(:,:), intent(in) :: u
            real, intent(in) :: dx, dt
            real, dimension(size(u,1), size(u,2)-1) :: ans, alpha_val
            real, dimension(size(u,1), size(u,2)) :: flux_val
            integer :: cols

            cols = size(u,2)

            call f(u, flux=flux_val)

            alpha_val = alpha(u=u, flux_val=flux_val, dx=dx, dt=dt)

            ans = 0.5*(flux_val(:,1:cols-1) + flux_val(:,2:cols)) &
                    + 0.5*abs(alpha_val)*(u(:,1:cols-1) - u(:,2:cols))

        end function flux

        function alpha(u, flux_val, dx, dt) result(ans)
            real, dimension(:,:), intent(in) :: u, flux_val
            real, intent(in) :: dx, dt
            real, dimension(size(u,1), size(u,2)-1) :: ans
            integer :: i, j, rows, cols
            real :: s
            s = dx/dt
            cols = size(u,2)

            ans = (flux_val(:,2:cols) - flux_val(:,1:cols-1))/(u(:,2:cols) - u(:,1:cols-1))
            
            do j = 1,cols-1
                do i = 1,8
                    if (u(i,j+1) - u(i,j) < 1e-10) then
                        ans(i,j) = s
                    end if
                end do
            end do

            call print_matrix(ans, "alpha")

        end function alpha

end submodule tvd_lxf_scheme
