for name in k_lnw_*
do
    newname=sst_"$(echo "$name" | cut -c7-)"
    mv "$name" "$newname"
done
