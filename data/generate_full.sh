awk < ELG_vista_des_target_jan13_1line5sig_i23.5_zprob.cat '{if($7==1)print "1 ",$0}' > new_coord_elg_TAG
awk < LRG_vista_des_target_jan13_sn2_i22.5_zphot_sup0.5.cat '{if($7==1)print "2 ",$0}'> new_coord_lrg_TAG
cat new_coord_elg_TAG new_coord_lrg_TAG > 20130129_gal_cat.dat


