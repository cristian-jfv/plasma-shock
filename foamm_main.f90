program lxf_kt

    use mesh, only: initialize_mesh, u, T, dx
    !use lxf_scheme, only: solve
    !use avg_grad_scheme, only: solve
    !use tvd_lxf_scheme, only: solve
    use scheme, only: solve
    use bookkeeping, only: save_matrix, save_vector
    use model, only: pg

    real, dimension(:), allocatable :: times
    real, dimension(:,:), allocatable :: total_p
    integer :: i

    call initialize_mesh(nr_ghost_cells=2)
    allocate(total_p(size(u, 2), size(u,3)))
    allocate(times(T))

    call solve(u=u, T=T, dx=dx, times=times, whole_b=1)

    call save_matrix(transpose(u(1,:,:)), "rho.txt")
    call save_matrix(transpose(u(2,:,:)), "rho_vx.txt")
    call save_matrix(transpose(u(3,:,:)), "rho_vy.txt")
    call save_matrix(transpose(u(4,:,:)), "rho_vz.txt")
    call save_matrix(transpose(u(5,:,:)), "Bx.txt")
    call save_matrix(transpose(u(6,:,:)), "By.txt")
    call save_matrix(transpose(u(7,:,:)), "Bz.txt")
    call save_matrix(transpose(u(8,:,:)), "e.txt")

    call save_vector(times, "times.txt")

    do i = 1, T
        total_p(:,i) = pg(e=u(8,:,i), B=u(5:7,:,i), vx=u(2,:,i)/u(1,:,i), &
                vy=u(3,:,i)/u(1,:,i), vz=u(4,:,i)/u(1,:,i), rho=u(1,:,i))
    end do

    call save_matrix(transpose(total_p), "p.txt")

    deallocate(u)
    deallocate(times)
    deallocate(total_p)

end program
