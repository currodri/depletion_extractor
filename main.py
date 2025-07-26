import argparse
from extractor import main as extract_main
from plotting import plot_depletion_comparison

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract depletion data into HDF5 format.")
    parser.add_argument("path_list_file", help="Text file containing simulation folder paths (one per line)")
    parser.add_argument("-f", "--filename", default="depletion_data.h5", help="Output HDF5 filename (default: depletion_data.h5)")
    parser.add_argument("-o", "--output_number", type=int, default=81, help="Output number for plotting (default: 81)")
    args = parser.parse_args()

    extract_main(args.path_list_file, args.filename)

    plot_depletion_comparison(args.filename, output_number=args.output_number)
