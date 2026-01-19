(1) prepare the target
    Generate param files for the ncaa and the ligand. The ncaa param file is used for rosetta sequence design and the ligand param file is used to generate rif.

(2) run rifgen for the fluorophore only
    ~/Rosetta/rifine/build/apps/rosetta/rifgen @rifgen_input/rifgen_req.flag

(3) dump the rif residues
    ---- for visualization
        ---- /home/longxing/Rosetta/rifine/build/apps/rosetta/devel_rifgen @rifgen_input/dump_rifres_req.flag -dump_mode RANDOM_DUMP -dump_fraction 0.00001 -dump_fname random_dump.pdb -dump_score_cutoff 0.0
    ---- for wefold
        ---- /home/longxing/Rosetta/rifine/build/apps/rosetta/devel_rifgen @rifgen_input/dump_rifres_req.flag -dump_mode DUMP_PHMAP -dump_fraction 1.0 -dump_fname BLA_rif_cart1.0_ang16.0_bound512.phmap -dump_cart_resl 1.0 -dump_ori_resl 16.0 -dump_max_cart 512 -dump_score_cutoff 0.0

(4) generate the protein backbone
    ./covalent_ligand_binder -len 120 -num_repeats 1 -symmetry C1 -frag_path ~/database/awesome_7HL_frags_20220427.bin -rpx_db ~/database/rpx_cart2.0_ang26.0_ss_ALLAA.phmap -rpx_cart_resl 2.0 -rpx_ang_resl 26.0 -gzip -use_ss -hallucinate -nstruct 200 -total_score_cutoff -13.0 -sidechain_neighbor_cutoff 2.40 -num_helix 5 -context_pdb ~/wuRFP/wefold_input/BLA_fluorophore_Pfr.pdb -rif_table ~/wuRFP/wefold_input/BLA_Pfr_rif_hbond-3.0_clipped_cart1.0_ang16.0_bound512.phmap -rif_cart_resl 1.0 -rif_ang_resl 16.0 -rif_score_cutoff -750 -rif_weight 10 -context_clash_weight 5.0 -context_clash_radius 1.4 -CB_swelling_factor 1.1 -output_dir ./2253/ -mcmc_outer_cycles 10 -mcmc_inner_cycles 20000 -ncaa ~/wuRFP/wefold_input/BLA_Pfr_helix.pdb -output_prefix wuRFPPfr -motif_ss HHHH -ligand_pdb ~/wuRFP/wefold_input/BLA_fluorophore_Pfr.pdb -ligand_neighbors_cutoff 1.4

(5) run fast design
    /home/caolongxingLab/caolongxing/Rosetta/rosetta_devel/source/cmake/build_hdf5/rosetta_scripts -s input.pdb -parser:protocol /home/caolongxingLab/caolongxing/wuRFP/scripts/fastdesign.xml @/home/caolongxingLab/caolongxing/wuRFP/scripts/fastdesign.flags -out:file:scorefile rosetta_scores.sc -out:path:all ./ -precompute_ig

(6) run mpnn
    ~/anaconda3/envs/pytorch/bin/python ~/wuRFP/scripts/mpnn_design.py -pdb_list input_pdbs.list -num_trials 10 -relax_cycles 3 -omit_AAs CX -temperature 0.05 -output_path ./ -scorefile mpnn.sc -prefix MPNN_
