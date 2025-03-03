#!/usr/bin/python3
import argparse
import sys
from coralme.builder.main import MEBuilder
import os, shutil
from pathlib import Path
#df_reaction_keff_consts
def parse_arguments():
    parser = argparse.ArgumentParser(description='coralME: COmprehensive Reconstruction ALgorithm for ME-models') 
    # Mandatory inputs
    parser.add_argument('--organism-json', help='Path to organism.json file')
    parser.add_argument('--m_model_path', help='Path to M-model file (.json or .xml)')
    parser.add_argument('--genbank_path', help='Path to GenBank file (.gb or .gbff)')
    
    # Optional parameters
    parser.add_argument('--run-blastp', action='store_true', help='Run BLASTp')
    parser.add_argument('--e-value', type=float, default=0.001, help='E-value cutoff')
    parser.add_argument('--locus-tag', default='locus_tag', help='Locus tag format (locus_tag or old_locus_tag)')
    parser.add_argument('--cores', type=int, default=1, help='Number of cores to use')
    parser.add_argument('--reference', help='Path to reference file')
    parser.add_argument('--include-pseudogenes', action='store_true', help='Include pseudogenes')
    parser.add_argument('--estimate-keffs', action='store_true', help='Estimate Keffs')
    parser.add_argument('--add-lipoproteins', action='store_true', help='Add lipoproteins')
    
    # Directory paths
    parser.add_argument('--log-directory', help='Path to logging directory')
    parser.add_argument('--out-directory', help='Path to output directory')
    
    # Optional file inputs
    parser.add_argument('--organism-matrix', help='Path to organism-specific matrix file')
    parser.add_argument('--tu-file', help='Path to Transcription Units file')
    parser.add_argument('--reaction-file', help='Path to Reaction file')
    parser.add_argument('--subreactions-file', help='Path to Subreactions file')
    parser.add_argument('--reactions-metadata', help='Path to Reactions metadata file')
    parser.add_argument('--metabolites-metadata', help='Path to Metabolites metadata file')
    
    # BioCyc related inputs
    parser.add_argument('--biocyc-genes', help='Path to BioCyc genes file')
    parser.add_argument('--biocyc-proteins', help='Path to BioCyc proteins file')
    parser.add_argument('--biocyc-tu', help='Path to BioCyc TU file')
    parser.add_argument('--biocyc-rna', help='Path to BioCyc RNA file')
    parser.add_argument('--biocyc-sequences', help='Path to BioCyc sequences file')
    
    return parser.parse_args()

if __name__ == "__main__":
    try:
        import coralme
        print(coralme.__file__)
        
        args = parse_arguments()
        print(f"Arguments: {args}")
        
        # Build configuration dictionary from args
        config = {
            "m-model-path": args.m_model_path,
            "genbank-path": args.genbank_path,
            "e_value_cutoff": args.e_value,
            "locus_tag": args.locus_tag,  # Now accepts either "locus_tag" or "old_locus_tag"
            "blast_threads": args.cores,
            "run_bbh_blast": args.run_blastp,
            "include_pseudo_genes": args.include_pseudogenes,
            "estimate_keffs": args.estimate_keffs,
            "add_lipoproteins": args.add_lipoproteins,
            "dev_reference": True
        }
        print(args.genbank_path)
        
        # Add optional parameters if provided
        if args.reference:
            config["reference-path"] = args.reference
        if args.log_directory:
            config["log_directory"] = args.log_directory
        if args.out_directory:
            config["out_directory"] = args.out_directory
        if args.organism_matrix:
            config["df_gene_cplxs_mods_rxns"] = args.organism_matrix
            
        # Add optional file inputs if provided
        if args.tu_file:
            config["transcription-units-file"] = args.tu_file
        if args.reaction_file:
            config["reaction-file"] = args.reaction_file
        if args.subreactions_file:
            config["subreactions-file"] = args.subreactions_file
        if args.reactions_metadata:
            config["reactions-metadata-file"] = args.reactions_metadata
        if args.metabolites_metadata:
            config["metabolites-metadata-file"] = args.metabolites_metadata
            
        # Add BioCyc inputs if provided
        if args.biocyc_genes:
            config["biocyc-genes-file"] = args.biocyc_genes
        if args.biocyc_proteins:
            config["biocyc-proteins-file"] = args.biocyc_proteins
        if args.biocyc_tu:
            config["biocyc-tu-file"] = args.biocyc_tu
        if args.biocyc_rna:
            config["biocyc-rna-file"] = args.biocyc_rna
        if args.biocyc_sequences:
            config["biocyc-sequences-file"] = args.biocyc_sequences
        
        # Print configuration for debugging
        print("Configuration:", config)
        
        # Initialize builder with configuration, including organism.json as the first parameter
        builder = MEBuilder(*[args.organism_json], **config)
        print(builder.configuration)
        
        # You can uncomment these lines to perform actual operations
        builder.generate_files(overwrite=True)
        shutil.copy2(args.organism_matrix, args.out_directory)
        builder.build_me_model(overwrite=False)
        # builder.troubleshoot(growth_key_and_value = { builder.me_model.mu : 0.001 })
        
        print("Script executed successfully.")
        sys.exit(0)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
