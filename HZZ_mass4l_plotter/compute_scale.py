import re
from collections import defaultdict


def parse_file(filename):
    """
    Parse input file into:
    data[bin][channel][process] = yield (float)
    """
    data = defaultdict(lambda: defaultdict(dict))
    current_bin = None

    with open(filename, "r") as f:
        for line in f:
            line = line.strip()

            # Detect bin names (all caps lines)
            if line.isupper():
                current_bin = line
                continue

            if not line:
                continue

            # Match data lines
            match = re.match(r"(\w+)\s+(\S+)\s+\$(.*?)\^", line)
            if not match:
                continue

            process = match.group(1)
            label = match.group(2)
            value = float(match.group(3))

            # Extract channel (last part after "_")
            channel = label.split("_")[-1]

            data[current_bin][channel][process] = value

    return data


def compute_scales(prefit, postfit):
    """
    Compute:
    - per-bin scales
    - total (inclusive) scales
    """
    scales = defaultdict(lambda: defaultdict(dict))

    # Flat totals
    totals_prefit = defaultdict(float)
    totals_postfit = defaultdict(float)

    for bin_name in prefit:
        for channel in prefit[bin_name]:
            for process in prefit[bin_name][channel]:

                if process in ["total", "background"]:
                    continue

                pre = prefit[bin_name][channel].get(process, 0.0)
                post = postfit[bin_name][channel].get(process, 0.0)

                # Per-bin scale
                scale = post / pre if pre > 0 else 0.0
                scales[bin_name][channel][process] = scale

                # Accumulate totals
                key = (channel, process)
                totals_prefit[key] += pre
                totals_postfit[key] += post

    # Compute total scales
    total_scales = {}
    for key in totals_prefit:
        pre = totals_prefit[key]
        post = totals_postfit[key]
        total_scales[key] = post / pre if pre > 0 else 0.0

    return scales, total_scales


def print_results(scales, total_scales, output_file="scales_output.txt"):
    with open(output_file, "w") as f:
        f.write("=== Per-bin scales ===\n\n")
        for bin_name in scales:
            f.write(f"[{bin_name}]\n")
            for channel in scales[bin_name]:
                for process, scale in scales[bin_name][channel].items():
                    f.write(f"{channel:6s} {process:10s} : {scale:.3f}\n")
            f.write("\n")

        f.write("=== Total (inclusive) scales ===\n\n")
        for (channel, process), scale in total_scales.items():
            f.write(f"{channel:6s} {process:10s} : {scale:.3f}\n")

    print(f"Results saved to {output_file}")


if __name__ == "__main__":
    prefit_file = "rates_prefit_zxestimate_new.log"
    postfit_file = "rates_postfit_zxestimate_new.log"
    output_file = "scales_output_final_new_dphijj.txt"

    print("Reading files...")
    prefit = parse_file(prefit_file)
    postfit = parse_file(postfit_file)

    print("Computing scales...")
    scales, total_scales = compute_scales(prefit, postfit)

    print("Saving results...")
    print_results(scales, total_scales, output_file)