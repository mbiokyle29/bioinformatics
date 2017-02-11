#!/usr/bin/env python
"""
listDiff.py
Kyle McChesney

Compare two sorted lists

"""

import logging, argparse, os, time

def main():
    # args
    parser = argparse.ArgumentParser(
        description = (" Compare two txt file lists "),
    )

    time_stamp = str(time.time()).replace(".","")
    log_file = os.path.join("./","{}.{}.{}".format("listDiff",time_stamp,"log"))

    parser.add_argument("--A", help="First list file", required=True)
    parser.add_argument("--B", help="Second list file", required=True)
    parser.add_argument("--output", help="Name of output file", default=log_file)
    args = parser.parse_args()

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log_formatter = logging.Formatter('%(asctime)s {%(levelname)s}: %(message)s')

    file_handler = logging.FileHandler(args.output)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(log_formatter)

    # console log
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(log_formatter)

    # set it all up
    log.addHandler(file_handler)
    log.addHandler(stream_handler)
    log.info("######### ListDiff ##########")
    
    a_list = [ line.rstrip('\n') for line in open(args.A)]
    b_list = [ line.rstrip('\n') for line in open(args.B)]

    log.info("Comparing %s to %s", args.A, args.B)
    log.info("%s has %i entries", args.A, len(a_list))
    log.info("%s has %i entries", args.B, len(b_list))
    log.info("The lists differ by %i entries", abs(len(a_list)-len(b_list)))

    a_set = set(a_list)
    b_set = set(b_list)

    log.info("The intersection of the lists has %i values", len(a_set.intersection(b_set)))
    log.info("The union of the lists has %i values", len(a_set.union(b_set)))

    if a_set.issubset(b_set):
        log.info("%s is a strict subset of %s", args.A, args.B)
    else:
        log.info("%s is NOT a strict subset of %s", args.A, args.B)

    if b_set.issubset(a_set):
        log.info("%s is a strict subset of %s", args.B, args.A)
    else:
        log.info("%s is NOT a strict subset of %s", args.B, args.A)

if __name__ == "__main__":
    main()