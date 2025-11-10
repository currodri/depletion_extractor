import argparse
from extractor import main as extract_main
from plotting import plot_depletion_comparison, plot_dtm_dtg_vs_metallicity


def parse_outputs(spec):
    """
    Parse output specification which can be:
      - single integer: "81"
      - comma separated list: "10,20,30"
      - range: "10-20" (inclusive)
    Returns a sorted list of unique ints.
    """
    if isinstance(spec, (int,)):
        return [spec]
    spec = str(spec).strip()
    if not spec:
        return []
    parts = []
    for part in spec.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            try:
                a, b = part.split('-', 1)
                a_i = int(a); b_i = int(b)
                if a_i <= b_i:
                    parts.extend(range(a_i, b_i+1))
                else:
                    parts.extend(range(b_i, a_i+1))
            except Exception:
                continue
        else:
            try:
                parts.append(int(part))
            except Exception:
                continue
    # unique and sorted
    return sorted(set(parts))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract depletion data into HDF5 format.")
    parser.add_argument("path_list_file", help="Text file containing simulation folder paths (one per line)")
    parser.add_argument("-f", "--filename", default="depletion_data.h5", help="Output HDF5 filename (default: depletion_data.h5)")
    parser.add_argument("-o", "--outputs", default="81", help="Output number(s) to extract. Single int, comma-separated list, or range like 10-20. Default: 81")
    args = parser.parse_args()

    outputs = parse_outputs(args.outputs)
    if not outputs:
        print("No valid outputs specified; exiting.")
        raise SystemExit(1)

    # Pass list of outputs to extractor (it accepts a single int or iterable)
    extract_main(args.path_list_file, args.filename, outputs)

    # Only call plotting helpers if a single output was requested (backwards compatible behavior)
    if len(outputs) == 1:
        out = outputs[0]
        plot_depletion_comparison(args.filename, output_number=out)
        plot_dtm_dtg_vs_metallicity(args.filename, output_number=out)
    else:
        print(f"Extracted outputs {outputs} into {args.filename}. Skipping automatic plotting for multiple outputs.")
