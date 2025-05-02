#!/usr/bin/python3
import argparse
from coralme.builder.main import MEBuilder
import sys
from pathlib import Path
import coralme

def parse_arguments():
    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    parser = argparse.ArgumentParser(description='coralME: COmprehensive Reconstruction ALgorithm for ME-models')
    # Mandatory inputs
    parser.add_argument('--m-model-path', help='Path to M-model file (.json or .xml)',required=True)
    parser.add_argument('--genbank-path', help='Path to GenBank file (.gb or .gbff)', required=True)

    # What modules to run
    parser.add_argument('--run-synchronize', help='Synchronize and Complement database files', default=True, type=str2bool)
    parser.add_argument('--run-build', help='Build the ME-model', default=True, type=str2bool)
    parser.add_argument('--run-troubleshoot', help='Troubleshoot the reconstructed ME-model', default=True, type=str2bool)

    # Optional parameters
    parser.add_argument('--organism-json', help='Path to organism.json configuration file', default=None)
    parser.add_argument('--run-bbh-blast', help='Run bi-directional BLASTp hit search', default=True)
    parser.add_argument('--e-value-cutoff', type=float, default=1e-10, help='E-value cutoff')
    parser.add_argument('--locus-tag', default='locus_tag', help='Locus tag format (locus_tag or old_locus_tag)')
    parser.add_argument('--blast-threads', type=int, default=1, help='Number of cores to use for BLASTp')
    parser.add_argument('--user-reference', help='Path to reference file', default=None)
    parser.add_argument('--include-pseudo-genes', help='Include pseudogenes', default=True)
    parser.add_argument('--estimate-keffs', help='Estimate Keffs', default=True)
    parser.add_argument('--add-lipoproteins', help='Add lipoproteins', default=True)
    parser.add_argument('--add-translocases', help='Add translocases associated to CPLX_dummy', default=False)

    # Directory paths
    parser.add_argument('--log-directory', help='Path to logging directory', default="./")
    parser.add_argument('--out-directory', help='Path to output directory', default="./")

    # Optional file inputs
    parser.add_argument('--df-gene-cplxs-mods-rxns', help='Path to organism-specific matrix file',default="./automated-org-with-refs.xlsx")
    parser.add_argument('--df-TranscriptionalUnits', help='Path to Transcription Units file', default=None)
    parser.add_argument('--df_matrix_stoichiometry', help='Path to Reaction file', default=None)
    parser.add_argument('--df_matrix_subrxn_stoich', help='Path to Subreactions file', default=None)
    parser.add_argument('--df_metadata_orphan_rxns', help='Path to Orphan Reactions file', default=None)
    parser.add_argument('--df_metadata_metabolites', help='Path to Metabolites mappings to E-matrix components file', default=None)

    # BioCyc related inputs
    parser.add_argument('--biocyc-genes', help='Path to BioCyc genes file', default=None)
    parser.add_argument('--biocyc-proteins', help='Path to BioCyc proteins file', default=None)
    parser.add_argument('--biocyc-tu', help='Path to BioCyc TU file', default=None)
    parser.add_argument('--biocyc-rna', help='Path to BioCyc RNA file', default=None)
    parser.add_argument('--biocyc-sequences', help='Path to BioCyc sequences file', default=None)

    return parser.parse_args()

def main():
    try:
        # import coralme
        # Only parse arguments if running from CLI
        config = {}
        args = parse_arguments()
        print(f"Arguments: {args}\n")
        for key, value in vars(args).items():
            if value is not None:
                if "biocyc" in key:
                    key = key.replace("biocyc_", "biocyc.")
                if "run_" in key:
                    continue
                config[key] = value
        # Print configuration for debugging
        print("Configuration:", config)

        # Initialize builder with configuration, including organism.json as the first parameter
        if args.organism_json:
            builder_args = [args.organism_json]
        else:
            builder_args = []
            config["ME-Model-ID"] = "coralME"

        builder = MEBuilder(*builder_args, **config)
        print(builder.configuration)

        # Run coralME modules based on command line arguments
        model_loaded = False
        if args.run_synchronize:
            builder.generate_files(overwrite=True)
            # builder.save_builder_info()

        if args.run_build:
            builder.build_me_model(overwrite=False)
            model_loaded = True

        if not model_loaded:
            builder.me_model = coralme.io.pickle.load_pickle_me_model(
                builder.configuration["out_directory"] + "MEModel-step2-{}.pkl".format(builder.configuration["ME-Model-ID"])
            )
            model_loaded = True

        if args.run_troubleshoot:
            builder.troubleshoot(growth_key_and_value={builder.me_model.mu: 0.001})

        print("Script executed successfully.")
        sys.exit(0)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
