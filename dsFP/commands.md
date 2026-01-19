(1) run rifgen to generate interacting residues with FDC
    /home/longxing/Rosetta/rifine/build/apps/rosetta/rifgen @rifgen_input/rifgen.flag


(2) dump rifres
    # dump the rif table as phmap map
    $ /home/longxing/Rosetta/rifine/build/apps/rosetta/devel_rifgen @rifgen_input/dump_rifres.flag -dump_mode DUMP_PHMAP -dump_fraction 1.0 -dump_fname FDC_rif_cart1.0_ang16.0_bound512.phmap -dump_cart_resl 1.0 -dump_ori_resl 16.0 -dump_max_cart 512 -dump_score_cutoff 0.0

    # dump the rif res for visualization
    $ /home/longxing/Rosetta/rifine/build/apps/rosetta/devel_rifgen @rifgen_input/dump_rifres.flag -dump_mode RANDOM_DUMP -dump_fraction 0.001 -dump_fname random_dump.pdb -dump_score_cutoff 0.0

(3) protein backbone generation via covalent_ligand_binder
    ./covalent_ligand_binder -len 88 -num_repeats 1 -symmetry C1 -frag_path ~/WeFold/database/awesome_7HL_frags_20220427.bin -rpx_db ~/WeFold/database/rpx_cart2.0_ang26.0_ss_ALLAA.phmap -rpx_cart_resl 2.0 -rpx_ang_resl 26.0 -gzip -dump_polyG -use_ss -hallucinate -nstruct 1 -total_score_cutoff -10.0 -sidechain_neighbor_cutoff 2.0 -num_helix 4 -context_pdb ~/Projects/AzoF/params/ABT.pdb -rif_table ~/Projects/AzoF/AFT_rif_cart1.0_ang16.0_bound512.phmap -rif_cart_resl 1.0 -rif_ang_resl 16.0 -rif_score_cutoff -10 -rif_weight 20 -context_clash_weight 5.0 -context_clash_radius 1.5 -output_dir ./ -mcmc_outer_cycles 10 -mcmc_inner_cycles 20000 -ncaa ~/Projects/AzoF/AFT_rotamers/AFT_27.pdb -ncaa_pos 40

(4) run Rosetta fastdesign
    /home/caolongxingLab/caolongxing/Rosetta/rosetta_devel/source/cmake/build_release/rosetta_scripts -s input.pdb -parser:protocol /home/caolongxingLab/caolongxing/IRFP/scripts/fastdesign.xml @/home/caolongxingLab/caolongxing/IRFP/scripts/fastdesign.flags -out:file:scorefile test.sc -out:path:all ./

(5) run ProteinMPNN to further optimize the protein sequence
    /home/caolongxingLab/caolongxing/anaconda3/envs/pytorch/bin/python /home/caolongxingLab/caolongxing/IRFP/scripts/mpnn_design.py -pdb_list input_pdbs.list -num_trials 5 -relax_cycles 2 -temperature 0.02 -omit_AAs CX -output_path ./ -scorefile 0688_mpnn1.sc -prefix MPNN_
