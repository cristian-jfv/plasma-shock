submodule (scheme) lxf_scheme
    use model, only: f, sound_speed, alfven_speed, fast_speed
    implicit none

    contains

        module function flux(u, dx, dt) result(ans)
            real, intent(in) :: u(:,:)
            real, intent(in) :: dx, dt
            real, dimension(size(u, 1), size(u, 2)-1) :: ans
            real, dimension(size(u, 1), size(u, 2)) :: flux_val
            integer :: cols
            cols = size(u, 2)

            call f(u, flux_val)
            !print *, "after f calculation"

            ans = 0.5*(flux_val(:,1:cols-1) + flux_val(:,2:cols)) &
                    + 0.5*(dx/dt)*(u(:,1:cols-1) - u(:,2:cols))

        end function flux

end submodule lxf_scheme
