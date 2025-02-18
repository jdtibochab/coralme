#!/usr/bin/python3

if __name__ == "__main__":
    import coralme
    print(coralme.__file__)

    #if len(sys.argv) < 4:
        #print("Usage: python script.py <file1_path> <file2_path> <output_file_path>")
        #sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    #output_file_path = sys.argv[3]

    #print(f"File 1 Path: {file1_path}")
    #print(f"File 2 Path: {file2_path}")
    #print(f"Output File Path: {output_file_path}")

    from coralme.builder.main import MEBuilder
    organism = './organism.json'
    inputs = './input.json'

    builder = MEBuilder(*[organism, inputs], **{ "m-model-path" : file1_path, "genbank-path" : file2_path })
    builder.generate_files(overwrite=True)
    builder.build_me_model(overwrite=False)
    builder.troubleshoot(growth_key_and_value = { builder.me_model.mu : 0.001 })

    # You can now proceed to read and process the files
    # using the provided paths
