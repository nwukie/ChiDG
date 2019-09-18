for name in bc_state_k_lnw_*
do
    newname=bc_state_sst_"$(echo "$name" | cut -c16-)"
    mv "$name" "$newname"
done
