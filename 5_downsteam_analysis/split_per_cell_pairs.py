import os
import argparse

#def split_file_by_first_column(input_file, output_dir, prefix):
def split_file_by_first_column(input_file, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    file_data = {}
    with open(input_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            columns = line.split('\t')
            first_col = columns[0]
            remaining_cols = '\t'.join(columns[1:])
            if first_col not in file_data:
                file_data[first_col] = []
            file_data[first_col].append(remaining_cols)
#    for first_col, lines in file_data.items():
    for first_col, lines in file_data.items():
        #output_file = os.path.join(output_dir, f"{prefix}_{first_col}.txt")
        output_file = os.path.join(output_dir, f"{first_col}.txt")
        with open(output_file, 'w') as out_f:
            out_f.write('\n'.join(lines) + '\n')
    #print(f"{prefix} finished")
    print("split finished")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split a file by its first column into multiple files.")
    parser.add_argument("--input_file", type=str, required=True, help="Path to the input file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the output files.")
    #parser.add_argument("--prefix", type=str, default="", help="Prefix to add to each output file name.")
    args = parser.parse_args()
    #split_file_by_first_column(args.input_file, args.output_dir, args.prefix)
    split_file_by_first_column(args.input_file, args.output_dir)

