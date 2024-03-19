module conserved_foamm_scheme
    use model, only: f
    implicit none

    contains

        function flux(u, dx, dt) result(ans)
            real, target, dimension(:,:), intent(in) :: u
            real, intent(in) :: dx, dt
            real, dimension(size(u, 1), size(u, 2)-1) :: ans, sumatory
            real, dimension(size(u, 1), size(u, 2)) :: model_flux

            call f(u, flux_val)

            sumatory = calc_sumatory(u, dx, dt)

            ans = 0.5*(flux_val(:,1:cols-1) + flux_val(:,2:cols)) &
                    - 0.5*sumatory
            
        end function flux

        function calc_sumatory(u, dx, dt) result(ans)
            real, target, dimension(:,:), intent(in) :: u
            real, intent(in) :: dx, dt

            real, dimension(7, size(u,1)) :: R
            real, dimension(7) :: lambda, eta
            real, dimension(size(u,1), size(u,2)-1) :: ans
            integer :: c, u_cols, k
            
            u_cols = size(u,2)
            
            ! Proceed column-wise
            do c = 1, u_cols
                ! Calculate eigenvalues
                ! Calculate characteristics
                ! Calculate right eigenvectors
            end do
            
        end function calc_sumatory

        function calc_eigenvalues(p) result(ans)
            ! Calculate the eigenvalues at every position given by the matrix p
            ! Asumes p is the averaged matrix at the interfaces of the cells
            real, target, dimension(:,:), intent(in) :: p
            real, dimension(size(p,2)) :: v, cs, ca, cf, a
            real, pointer, dimension(:) :: pg, rho, bx
            real, dimension(7,size(p,2)) :: ans
            !       1    2   3   4   5   6   7   8
            ! P = [rho, vx, vy, vz, Bx, By, Bz, pg]

            rho => p(1,:)
            bx => p(5,:)
            pg => p(8,:)

            v = sqrt(sum(p(2:4,:)**2, dim=1))

            ! a: sound speed
            a = sound_speed(pg, rho)
            ! ca: Alfven speed
            ca = alfven_speed(bx=bx, rho=rho)
            ! cf: Fast magnetosonic speed
            cf = fast_speed(a=a, va=ca)
            ! cs: Slow magnetosonic speed
            cs = slow_speed(a=a, va=ca)

            ! lambda_1 = v - cf
            ans(1,:) = v - cf

            ! lambda 2: v - ca
            ans(2,:) = v - ca

            ! lambda 3: v - cs
            ans(3,:) = v - cs

            ! lambda 4: v
            ans(4,:) = v
            
            ! lambda 5: v + cs
            ans(5,:) = v + cs

            ! lambda 6: v + ca
            ans(6,:) = v + ca
        
            ! lambda 7: v + cf
            ans(7,:) = v + cf
    
        end function calc_eigenvalues

        function calc_eigenvalues(u) result(ans)
            real, dimension(:,:), intent(in) :: u
            real, dimension(size(u,1), size(u,2) - 1) :: u_avg
            real, dimension(7) :: ans
            integer :: u_cols
            real :: cf, ca, cs

            u_cols = size(u,2)
            
            u_avg = 0.5*(u(:,1:u_cols-1) - u(:,2:u_cols))

            ! 1: u - cf
            ! 2: u - ca
            ! 3: u - cs
            ! 4: u
            ! 5: u + cs
            ! 6: u + ca
            ! 7: u + cf
            
        end function calc_eigenvalues

        function calc_characteristics(ul, ur) result(ans)
            real, dimension(8), intent(in) :: ul, ur
            real, dimension(8) :: u_avg
            real, dimension(7) :: ans
            
        end function calc_characteristics

end module conserved_foamm_scheme

