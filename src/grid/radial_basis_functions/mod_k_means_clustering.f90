module mod_k_means_clustering
    use mod_kinds,              only: ik, rk
    use mod_constants,          only: ZERO, ONE, TWO
    use mod_rbf_tools
    implicit none


contains

    ! Given an array of node coordinates and integer k, this function returns
    ! an array of the indices of the k nodes that are closest to the k-means cluster centers.
    function get_k_means_node_indices(nodes, k) result(center_indices)
        real(rk),       intent(in)  :: nodes(:,:)
        integer(ik),    intent(in)  :: k

        integer(ik) :: center_indices(k)


        real(rk) :: cluster_center(k,3), cluster_center_old(k,3)
        integer(ik) :: nnodes, inode
        integer(ik), allocatable :: cluster_ID(:)
        integer(ik) :: cluster_count(k)

        integer(ik) :: icluster, icenter, num_centers

        real(rk) :: cc_diff, cc_diff_max, dist_test, dist_min, dist_current, tol, dist_max

        logical :: clustering_incomplete = .true.

        nnodes = size(nodes,1)

        tol = 1.0e-5_rk

        allocate(cluster_ID(nnodes))


        !coord_min(1) = minval(nodes(:,1))
        !coord_min(2) = minval(nodes(:,2))
        !coord_min(3) = minval(nodes(:,3))

        !coord_max(1) = maxval(nodes(:,1))
        !coord_max(2) = maxval(nodes(:,2))
        !coord_max(3) = maxval(nodes(:,3))
        ! Randomly assign  initial cluster centers

        !
        ! Perform k-means iteration to find cluster indices
        !

        cluster_center(1,:) = nodes(1,:)
        cluster_ID(1) = 1
        
        num_centers = 1
        do while (num_centers<k)

            !Loop over the nodes and find the one with the maximal sum of distances to pre-existing cluster centers
            dist_max = 0.0_rk
            do inode = 1, nnodes

                dist_current = 0.0_rk
                do icluster = 1, num_centers

                    dist_current = dist_current + node_dist(nodes(inode,:), cluster_center(icluster,:))**TWO
                    
                end do

                if ((dist_current > dist_max) .and. (.not. any(cluster_ID==inode)) ) then

                    dist_max = dist_current
                    icenter = inode

                end if

            end do

            
            cluster_ID(num_centers) = icenter
            num_centers = num_centers+1
            cluster_center(num_centers,:) = nodes(icenter,:)
        end do

        !do icluster = 1, k

        !    cluster_center(icluster,:) = nodes(icluster,:)
        !end do
        cluster_center_old = cluster_center
        do while (clustering_incomplete)
            ! Assign nodes to clusters
            do inode = 1, nnodes
                
                
                dist_current = 1.0e15_rk
                do icluster = 1, k

                    dist_test = node_dist(nodes(inode,:), cluster_center(icluster,:))



                    if (dist_test < dist_current) then

                        cluster_ID(inode) = icluster
                        dist_current = dist_test
                        
                    end if


                end do


            end do

            ! Recompute cluster centers


            cluster_count = 0
            cluster_center = 0
            do inode = 1, nnodes


                icluster = cluster_ID(inode)
                cluster_count(icluster) = cluster_count(icluster)+1
                cluster_center(icluster,:) = cluster_center(icluster,:) + nodes(inode,:)


            end do


            cc_diff_max = 0
            do icluster = 1,k

                cluster_center(icluster,:) = &
                    (cluster_center(icluster,:)+cluster_center_old(icluster,:))/real(cluster_count(icluster)+1, rk)



                cc_diff = node_dist(cluster_center(icluster,:),cluster_center_old(icluster,:))

                if (cc_diff > cc_diff_max) cc_diff_max = cc_diff
            end do
            cluster_center_old = cluster_center
            
            clustering_incomplete = (cc_diff_max > tol)

        end do

        ! Loop over each cluster and see which input node is closest to the cluster center, and return this node

        do icluster = 1, k


            dist_min = 1.0e30_rk
            do inode = 1, nnodes

                if (cluster_ID(inode) == icluster) then

                    
                    dist_test = node_dist(nodes(inode,:), cluster_center(icluster,:))

                    if (dist_test < dist_min) then

                        dist_min = dist_test
                        center_indices(icluster) = inode



                    end if

                end if

            end do

        end do

    end function get_k_means_node_indices




end module mod_k_means_clustering
