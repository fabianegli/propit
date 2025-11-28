# © Fabian Egli, FGCZ, ETHZ
# GPLv3+
"""Console script for propit."""

from pathlib import Path

from . import convert


def main():
    import argparse

    def existing_directory(directory):
        directory = Path(directory)
        if not directory.is_dir():
            raise argparse.ArgumentTypeError(f"'{directory.absolute}' does not exist.")
        return directory

    def existing_pin_file(pin_file):
        pin_file = Path(pin_file)
        if not pin_file.is_file():
            raise argparse.ArgumentTypeError(f"'{pin_file.absolute}' does not exist.")
        return pin_file

    def float_between_0_and_1(strimput):
        try:
            value = float(strimput)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"'{strimput}' is not a valid float.") from exc
        if not 0 <= value <= 1:
            raise argparse.ArgumentTypeError(f"'{value}' is not in the expected range of [0,1].")
        return value

    parser = argparse.ArgumentParser(description="Convert Proteomics files.")

    subparsers = parser.add_subparsers(required=True)

    # percolator -> FlashLFQ
    percolator2flashlfq = subparsers.add_parser("percolator2flashlfq", help="Make generic input files for FlashLFQ.")
    percolator2flashlfq.add_argument(
        "--percolator-output-dir",
        dest="pout",
        action="store",
        type=existing_directory,
        required=True,
    )
    percolator2flashlfq.add_argument(
        "--max-q-value",
        dest="mq",
        action="store",
        type=float_between_0_and_1,
        default=0.01,
    )
    percolator2flashlfq.add_argument("--remove-contaminants", dest="rc", action="store_true", default=True)
    percolator2flashlfq.add_argument("--leave-contaminants", dest="rc", action="store_false")
    percolator2flashlfq.add_argument("--generic-flashlfq-file-name", dest="gf", action="store", default=None)
    percolator2flashlfq.add_argument("--comet_perc", dest="comet_perc", action="store_true", default=False)

    def p2f(args):
        convert.Percolator2FlashLFQ(
            percolator_output_dir=args.pout,
            max_q_value=args.mq,
            remove_contaminants=args.rc,
            comet_perc=args.comet_perc,
        ).write(generic_flashlfq_input_file=args.gf)

    percolator2flashlfq.set_defaults(func=p2f)

    # Sage -> FlashLFQ
    sage2flashlfq = subparsers.add_parser("sage2flashlfq", help="Make generic input files for FlashLFQ.")
    sage2flashlfq.add_argument(
        "--sage-output-dir",
        dest="sout",
        action="store",
        type=existing_directory,
        required=True,
    )
    sage2flashlfq.add_argument("--flash-wd", dest="fwd", action="store", type=existing_directory, required=True)

    def s2f(args):
        convert.sage_results_2_generic_flashlfq_input(args.sout, flashlfq_wd=args.fwd)

    sage2flashlfq.set_defaults(func=s2f)

    # FlashLFQ -> ProteoBench
    flfq_pep2pb = subparsers.add_parser(
        "flashlfq2pb",
        help="Make generic input files for ProteoBench from FlashLFQ peptides.",
    )
    flfq_pep2pb.add_argument(
        "--flashlfq-output-dir",
        dest="fout",
        action="store",
        type=existing_directory,
        required=True,
    )

    def fpep2pb(args):
        convert.flashlfq_peptides_to_proteobench(args.fout)

    flfq_pep2pb.set_defaults(func=fpep2pb)

    # Comet -> percolator
    comet2perc = subparsers.add_parser(
        "comet2perc",
        help="This tool removes the deltaCn column from pin files.",
        description="Some comet versions sporadically write nans into the deltCn column which tripps off percolator.",
    )
    comet2perc.add_argument("--pin-file", dest="pin", action="store", type=existing_pin_file, required=True)

    def com2per(args):
        convert.strip_deltacn_col(args.pin)

    comet2perc.set_defaults(func=com2per)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
